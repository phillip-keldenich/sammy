import slurminade
import lzma
import json
import subprocess
from pathlib import Path
import sys
import random
import tqdm


slurminade.update_default_configuration(
    partition="alg", exclusive=True, constraint="alggen05"
)


@slurminade.slurmify()
def run_commands(commands, output_file_flag):
    for command in commands:
        try:
            output_file = command[command.index(output_file_flag) + 1]
            subprocess.run(command, capture_output=True, check=True, text=True)
            subprocess.run(['xz', '-e', '-9', '-T', '24', str(output_file)],
                           capture_output=True, check=True, text=True)
        except CalledProcessError as e:
            with lzma.open(output_file, "wt") as f:
                error_data = {
                    'command': command,
                    'outcome': 'error',
                    'stdout': str(e.output),
                    'stderr': str(e.stderr),
                    'exitcode': e.returncode
                }
                json.dump(error_data, f)


if __name__ == "__main__":
    output_directory = Path(__file__).parent / "02_output"
    if not output_directory.is_dir():
        raise RuntimeError("Output directory does not exist - run step 1 experiment first")

    amender_file = Path(__file__).absolute().parent.parent.parent / "build" /\
                   "Release" / "src" / "amend_subproblem_with_default_mes"
    if not amender_file.is_file():
        raise RuntimeError("Could not find amend_subproblem_with_default_mes executable!")


    commands = []
    result_files = list(output_directory.glob("*.json.xz"))
    for result_file in tqdm.tqdm(result_files):
        name = result_file.name[:-8]
        subproblem_dir = output_directory / name
        if not subproblem_dir.is_dir():
            print(f"Subproblem directory does not exist for instance {name}")
            continue
        universe_and_clauses = subproblem_dir / 'universe_and_clauses.json.xz'
        for subproblem_file in subproblem_dir.glob("subproblem-*.json.xz"):
            output_file = subproblem_file.parent / subproblem_file.name.replace("subproblem-", "amended-").replace(".json.xz", ".json")
            output_xz = output_file.with_suffix(".json.xz")
            if output_xz.is_file():
                continue
            command = [str(amender_file), '--universe-and-clauses', str(universe_and_clauses),
                       '--subproblem', str(subproblem_file), '--output', str(output_file)]
            commands.append(command)
    random.shuffle(commands)
    with slurminade.JobBundling(max_size=4):
        while commands:
            current_batch = []
            while len(current_batch) < 1000 and commands:
                current_batch.append(commands.pop())
            run_commands.distribute(current_batch, '--output')

