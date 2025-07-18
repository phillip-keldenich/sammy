import os
import subprocess
from subprocess import CalledProcessError
import lzma
import json
import os
import random
import slurminade
import sys
from pathlib import Path


slurminade.update_default_configuration(
    partition="alg", exclusive=True, constraint="alggen05"
)


@slurminade.slurmify()
def run_commands(commands, output_file_flag, output_dir_flag):
    for command in commands:
        try:
            output_file = command[command.index(output_file_flag) + 1]
            output_dir = command[command.index(output_dir_flag) + 1]
            subprocess.run(command, capture_output=True, check=True, text=True)
            subprocess.run(['xz', '-e', '-9', '-T', '24', str(output_file)], capture_output=True, check=True, text=True)
            for file in Path(output_dir).glob("*.json"):
                subprocess.run(['xz', '-e', '-9', '-T', '24', str(file)], capture_output=True, check=True, text=True)
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
    instances_dir = "../06_identify_nontrivial/06_nontrivial_instances"
    output_dir = "02_output"
    command_path = "../../build/Release/src/sammy_solve"
    os.makedirs(output_dir, exist_ok=True)

    # avoid OOM conditions
    lns_worker_limits = {
        "freetz": 2,
        "FreeBSD-8_0_0": 2,
        "ocelot": 3,
        "uClinux": 7,
    }

    commands = []
    for f in os.listdir(instances_dir):
        if not f.endswith(".scm.json.xz"):
            continue
        instance_name = f[:-12]
        instance_path = os.path.join(instances_dir, f)
        subproblem_dir = os.path.join(output_dir, instance_name)
        output_file = subproblem_dir + ".json"
        if os.path.isfile(output_file + ".xz"):
            print("Skipping:", instance_name)
            continue
        extra = []
        if instance_name in lns_worker_limits:
            extra = ["--max-lns-workers", str(lns_worker_limits[instance_name])]
        commands.append([command_path, "--print-events", "--print-global-stats", "-o", output_file,
                        "--report-subproblems-to", subproblem_dir, *extra, instance_path])
        os.makedirs(subproblem_dir, exist_ok=True)

    with slurminade.JobBundling(max_size=1):
        while commands:
            next_batch = []
            while commands and len(next_batch) < 1:
                next_batch.append(commands.pop())
            run_commands.distribute(next_batch, '-o', '--report-subproblems-to')

