import lzma
import json
from pathlib import Path
import subprocess
import tqdm
import re

this_dir = Path(__file__).absolute().parent
output = this_dir / "02_output"
summary_file = this_dir / "05_amended_outcomes.json.xz"

out_of_iterations = []
mes_precludes_improvement = []
mes_optimal_does_not_preclude = []
mes_sizes = {}
rem_sizes = {}
pre_handled = set()
mes_re = re.compile(r'"mes_size": +([0-9]+)[^0-9]+$')
rem_re = re.compile(r'"num_removed_configs": +([0-9]+)[^0-9]+$')


def extract_mes_size(filename):
    res = subprocess.run(f'xzcat "{str(filename)}" | grep -Eo \'("mes_size"|"num_removed_configs"): [0-9]+[^0-9]\'',
                         shell=True, capture_output=True, text=True)
    max_mes_size = None
    rem_size = None
    for line in res.stdout.splitlines():
        match = mes_re.match(line)
        if match:
            max_mes_size = max(max_mes_size, int(match[1])) if max_mes_size else int(match[1])
        match = rem_re.match(line)
        if match:
            rem_size = int(match[1])
    if max_mes_size is None or rem_size is None:
        raise RuntimeError(str(filename) + ": " + res.stdout)
    return max_mes_size, rem_size


try:
    if summary_file.is_file():
        with lzma.open(summary_file, "rt") as f:
            data = json.load(f)
            out_of_iterations = data['out_of_iterations']
            mes_precludes_improvement = data['mes_precludes']
            mes_optimal_does_not_preclude = data['mes_optimal']
            mes_sizes = data['mes_sizes']
            rem_sizes = data['removal_sizes']
            for s in [out_of_iterations, mes_precludes_improvement, mes_optimal_does_not_preclude]:
                for f in s:
                    pre_handled.add(str(f))

    for subproblem in tqdm.tqdm(list(output.glob("*/amended-*.xz"))):
        if str(subproblem) in pre_handled:
            continue

        # avoid huge memory usage via grep and xzcat (most of the time)
        mes_sizes[str(subproblem)], rem_sizes[str(subproblem)] = extract_mes_size(subproblem)
        res = subprocess.run(f'xzcat "{str(subproblem)}" | grep -o LNS_MES_BY_CNP_PRECLUDES_IMPROVEMENT',
                             shell=True, capture_output=True)
        out = res.stdout.strip()
        if out:
            mes_precludes_improvement.append(subproblem)
            continue

        res2 = subprocess.run(f'xzcat "{str(subproblem)}" | grep -o LNS_MES_STOPPING_DUE_TO_MAX_ITERATIONS',
                              shell=True, capture_output=True)
        out2 = res2.stdout.strip()
        if out2:
            out_of_iterations.append(subproblem)
            continue

        try:
            with lzma.open(subproblem, "rt") as f:
                data = json.load(f)
                num_removed = data['num_removed_configs']
                mes_size = data['amendment']['mes_size']
                if num_removed <= mes_size:
                    raise RuntimeError(f"MES stopped for weird reason ({subproblem})")
                try:
                    if data['amendment']['events'][-2]['type'] == "DONE_SOLVE_FULL_RELAXATION":
                        if 'status' not in data['amendment']['events'][-2]:
                            out_of_iterations.append(subproblem)
                            continue
                        if data['amendment']['events'][-2]['status'] == "optimum on subgraph":
                            mes_optimal_does_not_preclude.append(subproblem)
                            continue
                    if data['amendment']['events'][-3]['type'] == "DONE_SOLVE_FULL_RELAXATION":
                        if 'status' not in data['amendment']['events'][-3]:
                            out_of_iterations.append(subproblem)
                            continue
                        if data['amendment']['events'][-3]['status'] == "optimum on subgraph":
                            mes_optimal_does_not_preclude.append(subproblem)
                            continue
                    if data['amendment']['events'][-2]['type'] == "ADDED_VERTICES" and \
                        data['amendment']['events'][-2]['new_count'] == 0:
                        mes_optimal_does_not_preclude.append(subproblem)
                        continue
                except KeyError as k:
                    raise RuntimeError(f"KEY ERROR: {subproblem}") from k
                raise RuntimeError(f"UNEXPECTED EVENT: {subproblem}")
        except Exception as e:
                raise RuntimeError(f"UNKNOWN ERROR: {subproblem}, {str(e)}") from e
finally:
    with lzma.open(summary_file, "wt") as f:
        json.dump({
            "out_of_iterations": [str(f) for f in out_of_iterations],
            "mes_precludes": [str(f) for f in mes_precludes_improvement],
            "mes_optimal": [str(f) for f in mes_optimal_does_not_preclude],
            "mes_sizes": mes_sizes,
            "removal_sizes": rem_sizes
        }, f)
