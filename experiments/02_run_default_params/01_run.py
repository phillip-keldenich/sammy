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
    output_file = command[command.index("-o") + 1]
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
                      instance_filter, max_workers=None):
    commands = []
    for f in os.listdir(instances_dir):
        if not f.endswith(".scm.json.xz"):
            continue
        if not instance_filter(f):
            continue
        instance_name = f[:-12]
        instance_path = os.path.join(instances_dir, f)
        for repeat in range(1, RUN_REPEATS + 1):
            output_file = os.path.join(output_dir, f"{instance_name}_repeat{repeat}.out.json")
            if os.path.exists(output_file + ".xz"):
                print(f"Skipping {output_file} as it already exists.")
                continue
            commands.append([command_path, instance_path, "-o", output_file])
            if max_workers is not None:
                commands[-1] += ['--max-lns-workers', str(max_workers)]
    return commands


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Run either with command 'small_instances' (via slurm) or 'huge_instances' (directly on an algra)")
        sys.exit(0)
    
    instances_dir = "../../sammy_benchmark_instances"
    output_dir = "02_output"
    command_path = "../../build/Release/src/sammy_solve"
    os.makedirs(output_dir, exist_ok=True)

    if sys.argv[1] == 'small_instances':
        # small instances that can be run regularly via slurm 
        # and don't risk OOM killing at 88 GiB by slurm
        commands = generate_commands(
            instances_dir=instances_dir,
            output_dir=output_dir,
            command_path=command_path,
            instance_filter=lambda f: ("Automotive02" not in f)
        )
        random.shuffle(commands)
        with slurminade.JobBundling(max_size=3):
            for command in commands:
                run_command.distribute(command)
    elif sys.argv[1] == 'huge_instances':
        # huge instances that may risk OOM killing
        if not socket.gethostname().startswith("algra"):
            print("Run 'huge_instances' directly on an algra after reserving it with sleep")
            sys.exit(1)
        commands = generate_commands(
            instances_dir=instances_dir,
            output_dir=output_dir,
            command_path=command_path,
            instance_filter=lambda f: ("Automotive02" in f),
            max_workers=1  # limit to 1 worker to avoid OOM (some more would prob work)
        )
        for command in commands:
            print(f"Running command: {command}")
            # run directly on algra, no slurmification needed
            run_command(command)
    else:
        print("Run either with command 'small_instances' (via slurm) or 'huge instances' (directly on an algra)")
        sys.exit(1)
