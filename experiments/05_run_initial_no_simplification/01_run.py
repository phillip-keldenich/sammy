import os
import subprocess
import sys
import slurminade
import traceback
import random
import socket
import json


# number of runs per instance
RUN_REPEATS = 5


# ensure right systems and exclusive use
slurminade.update_default_configuration(
    partition="alg", exclusive=True, constraint="alggen05"
)


# run a command, but allow distribution
@slurminade.slurmify()
def run_command(command):
    output_file = command[command.index("--dump-initial-phase") + 1]
    try:
        subprocess.run(command, capture_output=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command {command}: {e}")
        print(f"Output: {e.output.decode()}")
        print(f"Error output: {e.stderr.decode()}")
        with open(output_file, "w") as f:
            json.dump({
                "outcome": "error",
                "command": command,
                "returncode": e.returncode,
                "error": str(e),
                "output": e.output.decode(),
                "stderr": e.stderr.decode(),
            }, f)
    subprocess.run(["xz", "-9", "-e", output_file], check=True, 
                   capture_output=True)


def generate_commands(instances_dir, output_dir, command_path,
                      instance_filter):
    commands = []
    for f in os.listdir(instances_dir):
        if not f.endswith(".scm.json.xz"):
            continue
        if not instance_filter(f):
            continue
        instance_name = f[:-12]
        instance_path = os.path.join(instances_dir, f)
        for repeat in range(1, RUN_REPEATS + 1):
            output_file = os.path.join(output_dir, f"{instance_name}-nosimp_repeat{repeat}.out.json")
            if os.path.exists(output_file + ".xz"):
                print(f"Skipping {output_file} as it already exists.")
                continue
            commands.append([command_path, instance_path, "--dont-simplify",
                             "--dump-initial-phase", output_file])
    return commands


if __name__ == "__main__":
    instances_dir = "../../sammy_benchmark_instances"
    output_dir = "02_output"
    command_path = "../../build/Release/src/sammy_solve"
    os.makedirs(output_dir, exist_ok=True)
    # small instances that can be run regularly via slurm 
    # and don't risk OOM killing at 88 GiB by slurm
    commands = generate_commands(
        instances_dir=instances_dir,
        output_dir=output_dir,
        command_path=command_path,
        instance_filter=lambda f: ("Automotive02" not in f)
    )
    # these run into memory issues; it may be necessary to
    # run them manually, allowing swap or at least the full
    # 96 GiB instead of only 88 GiB that our SLURM allows
    commands += generate_commands(
        instances_dir=instances_dir,
        output_dir=output_dir,
        command_path=command_path,
        instance_filter=lambda f: ("Automotive02" in f),
    )
    random.shuffle(commands)
    with slurminade.JobBundling(max_size=5):
        for command in commands:
            run_command.distribute(command)

