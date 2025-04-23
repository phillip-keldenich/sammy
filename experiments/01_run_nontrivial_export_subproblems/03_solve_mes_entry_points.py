import os
import subprocess
import sys
from pathlib import Path
import slurminade
import traceback
import random
import lzma
import json
from filelock import FileLock
from tempfile import TemporaryDirectory
import tqdm


IMPORTANT_SOLVERS = [
    "fixed_sat[kissat]",
    "fixed_sat[cadical]",
    "fixed_sat[cryptominisat]",
    "fixed_sat[lingeling]",
    "incremental_sat[cadical]",
    "incremental_sat[cryptominisat]",
    "incremental_sat[lingeling]",
    "satdsatur[cadical]",
    "satdsatur[cryptominisat]",
    "satdsatur[lingeling]",
]
RUNS_PER_SOLVER = 3
EXECUTABLE = "../../build/Release/src/solve_subproblems_mes_entry_point"
TIME_LIMIT = 600


slurminade.update_default_configuration(
    partition="alg", exclusive=True, constraint="alggen05"
)
slurminade.set_dispatch_limit(100_000)


def mes_entry_point_files():
    for d in os.listdir("02_output"):
        path = Path(os.path.join("02_output", d))
        if not os.path.isdir(path):
            continue
        yield from path.glob("mes_entry_points-*.json.xz")


def mes_sizes(mes_entry_data):
    yield len(mes_entry_data["initial_mes"])
    if "mes_history" in mes_entry_data and mes_entry_data["mes_history"]:
        for mhist in mes_entry_data["mes_history"]:
            yield len(mhist["mes"])


def unresolved_mes_solver_pairs():
    for f in mes_entry_point_files():
        with lzma.open(f, "rt") as fin:
            mes_entry_data = json.load(fin)
        print(f)
        solutions_data = None
        solutions_file = f.parent / f.name.replace("mes_entry_points-", "solutions-")
        if solutions_file.exists():
            with lzma.open(solutions_file, "rt") as fin:
                try:
                    solutions_data = json.load(fin)
                except Exception as e:
                    print("ERROR READING:", solutions_file)
                    raise e
        if not solutions_data or "solver_runs" not in solutions_data or \
           not solutions_data["solver_runs"]:
            for mes_size in mes_sizes(mes_entry_data):
                for solver in IMPORTANT_SOLVERS:
                    yield (f, solutions_file, mes_size, solver, RUNS_PER_SOLVER)
        else:
            solver_and_size_counts = {}
            for solver_run in solutions_data["solver_runs"]:
                for solver_name, ind_runs in solver_run.items():
                    for ind_run in ind_runs:
                        mes_size = ind_run["mes_size"]
                        key = (solver_name, mes_size)
                        if ind_run['outcome'] == 'OOM':
                            solver_and_size_counts[key] = RUNS_PER_SOLVER + 100
                        else:
                            v = solver_and_size_counts.get(key, 0) + 1
                            solver_and_size_counts[key] = v
            for mes_size in mes_sizes(mes_entry_data):
                for solver in IMPORTANT_SOLVERS:
                    key = (solver, mes_size)
                    existing = solver_and_size_counts.get(key, 0)
                    if existing < RUNS_PER_SOLVER:
                        yield (f, solutions_file, mes_size, solver, 
                               RUNS_PER_SOLVER - existing)


def get_run_mes_size(solver_run):
    for _, ind_runs in solver_run.items():
        return ind_runs[0]["mes_size"]
    raise RuntimeError("No run found in solver run!")


def create_solutions_file(solutions_file):
    sfname = solutions_file.name
    parent = solutions_file.parent
    entry_point_file = parent
    entry_point_file /= sfname.replace("solutions-", "mes_entry_points-")
    formula_file = parent
    formula_file /= "universe_and_clauses.json.xz"
    subproblem_file = parent
    subproblem_file /= sfname.replace("solutions-", "subproblem-")
    with lzma.open(subproblem_file, "rt") as f:
        removal_size = json.load(f)['num_removed_configs']
    toplevel = {
        "entry_point_file": str(entry_point_file),
        "formula_file": str(formula_file),
        "solve_time_limit": TIME_LIMIT,
        "solver_runs": [],
        "subproblem_file": str(subproblem_file),
        "removal_size": removal_size
    }
    with lzma.open(solutions_file, "wt") as f:
        json.dump(toplevel, f, indent=2)


def add_given_outcome(solutions_file, mes_size, solver, outcome):
    with FileLock(str(solutions_file) + ".lock"):
        if not solutions_file.exists():
            create_solutions_file(solutions_file)
        with lzma.open(solutions_file, "rt") as f:
            solutions_data = json.load(f)
        solver_runs = solutions_data["solver_runs"]
        found_run = None
        for solver_run in solver_runs:
            if get_run_mes_size(solver_run) == mes_size:
                found_run = solver_run
                break
        if found_run is None:
            solver_runs.append(dict())
            found_run = solver_runs[-1]
        if not isinstance(outcome, list):
            found_run.setdefault(solver, []).append(outcome)
        else:
            found_run.setdefault(solver, []).extend(outcome)
        with lzma.open(solutions_file, "wt") as f:
            json.dump(solutions_data, f, indent=2)


def add_oom_outcome(solutions_file, mes_size, solver):
    oom_outcome = {"mes_size": mes_size, "outcome": "OOM", 
        "build_time": TIME_LIMIT, "solve_time": TIME_LIMIT,
        "configurations": None, "events": []}
    add_given_outcome(solutions_file, mes_size, solver, oom_outcome)


def add_outcome(solutions_file, mes_size, solver, output):
    with lzma.open(output, "rt") as f:
        solver_run = json.load(f)
    for ind_run in solver_run:
        ind_run["mes_size"] = mes_size
    add_given_outcome(solutions_file, mes_size, solver, solver_run)


def run_mes_solver_pair(temp_dir, mes_entry_file, solutions_file, 
                        mes_size, solver, num_runs):
    output = os.path.join(temp_dir, f"solutions.json.xz")
    completed = subprocess.run([EXECUTABLE, 
                                "--entry-point-file", str(mes_entry_file),
                                "--mes-size", str(mes_size),
                                "--solver", solver,
                                "--solve-time-limit", str(TIME_LIMIT),
                                "--num-runs", str(num_runs),
                                "-o", str(output)],
                                capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        if completed.returncode == -9:
            # OOM outcome
            add_oom_outcome(solutions_file, mes_size, solver)
        else:
            print("Error in called process!")
            print(completed.stdout)
            print(completed.stderr)
            raise RuntimeError(f"Solver {solver} failed with return " +
                               f"code {completed.returncode}")
    else:
        # regular outcome
        add_outcome(solutions_file, mes_size, solver, Path(output))
        os.remove(output)


@slurminade.slurmify()
def run_mes_solver_pairs(entries):
    with TemporaryDirectory() as temp_dir:
        for mes_entry_file, sol_file, mes_size, solver, num_runs in entries:
            mes_entry_file = Path(mes_entry_file)
            sol_file = Path(sol_file)
            run_mes_solver_pair(temp_dir, mes_entry_file, sol_file, mes_size,
                                solver, num_runs)


if __name__ == "__main__":
    ur = []
    for mef, sf, ms, sol, nr in unresolved_mes_solver_pairs():
        ur.append((str(mef), str(sf), ms, sol, nr))
        if len(ur) >= 200:
            run_mes_solver_pairs.distribute(ur)
            ur = []
    if ur:
        run_mes_solver_pairs.distribute(ur)
