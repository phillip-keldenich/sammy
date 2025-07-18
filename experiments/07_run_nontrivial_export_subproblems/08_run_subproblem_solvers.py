import json
import slurminade
from pathlib import Path
import random
import subprocess


# ensure right systems and exclusive use
slurminade.update_default_configuration(
    partition="alg", exclusive=True, constraint="alggen05"
)
slurminade.set_dispatch_limit(10000)


@slurminade.slurmify()
def run_command(command):
    output_file = command[command.index("-o") + 1]
    try:
        subprocess.run(command, capture_output=True, check=True)
        subprocess.run(["xz", "-9", "-e", "-T", "8", output_file], check=True, 
                       capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command {command}: {e}")
        print(f"Output: {e.output.decode()}")
        print(f"Error output: {e.stderr.decode()}")


if __name__ == "__main__":
    SUBPROBLEM_SOLVERS = [
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
        "newsatdsatur[cadical]",
        "newsatdsatur[cryptominisat]",
        "newsatdsatur[lingeling]"
    ]
    NUM_SUBPROBLEMS = 500
    EXECUTABLE = Path(__file__).absolute().parent.parent.parent /\
                 "build" / "Release" / "src" / "solve_amended_subproblem"
    if not EXECUTABLE.is_file():
        raise RuntimeError(f"Could not find solver executable {EXECUTABLE}!")
    with open("07_shuffled_nonoptimal.json", "r") as file:
        subproblem_files = json.load(file)
        subproblems = subproblem_files[:NUM_SUBPROBLEMS]
        commands = []
        for subproblem in subproblems:       
            subfile = Path(subproblem)
            if not subfile.is_file():
                raise RuntimeError("Missing subproblem file!")
            for solver in SUBPROBLEM_SOLVERS:
                solvername = solver.replace("[", "_").replace("]", "")
                output_file_name = subfile.name.replace("amended-", f"solved-{solvername}-")
                raw_output_name = output_file_name.replace(".xz", "")
                output_file = subfile.parent / output_file_name
                if output_file.is_file():
                    continue
                raw_output_file = subfile.parent / raw_output_name
                command = [str(EXECUTABLE), '--solver-name', solver, '-o', str(raw_output_file), '--amended-subproblem-file', str(subfile)]
                commands.append(command)
        random.shuffle(commands)
        with slurminade.JobBundling(max_size=50):
            for command in commands:
                run_command.distribute(command)

