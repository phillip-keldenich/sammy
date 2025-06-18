import os
import subprocess
import sys
import slurminade
import traceback
import random
import socket
import json


# ensure right systems and exclusive use
slurminade.update_default_configuration(
    partition="alg", exclusive=True, constraint="alggen05"
)

# run a command, but allow distribution
@slurminade.slurmify()
def run_command(command):
    try:
        subprocess.run(command, capture_output=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {command}")
        print(f"Return code: {e.returncode}")
        print(f"Output: {e.stdout.decode()}")
        print(f"Error: {e.stderr.decode()}")
        file = command[command.index("--dump-initial-phase") + 1]
        with open(file, "w") as f:
            data = {
                "error": "Command failed",
                "command": command,
                "return_code": e.returncode,
                "stdout": e.stdout.decode(),
                "stderr": e.stderr.decode(),
            }
            json.dump(data, f, indent=2)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Run either with command 'small_instances' (via slurm) or 'huge_instances' (directly on an algra)")
        sys.exit(0)

    instances_dir = "../../sammy_benchmark_instances"
    output_dir = "02_output"
    command_path = "../../build/Release/src/sammy_solve"
    os.makedirs(output_dir, exist_ok=True)
    commands = []

    if sys.argv[1] == 'small_instances':
        # small instances that can be run regularly via slurm 
        # and don't risk OOM killing
        for f in os.listdir(instances_dir):
            if not f.endswith(".scm.json"):
                continue
            if "Automotive02" in f:
                continue
            instance_name = f[:-9]
            instance_path = os.path.join(instances_dir, f)
            output_file = os.path.join(output_dir, f"{instance_name}-initial-nosimp.out.json")
            if os.path.exists(output_file):
                continue
            commands.append([command_path, instance_path, 
                             "--dump-initial-phase", output_file,
                             "--dont-simplify"])
        with slurminade.JobBundling(max_size=1):
            for command in commands:
                run_command.distribute(command)
    elif sys.argv[1] == 'huge_instances':
        # huge instances that may risk OOM killing without simplification
        if not socket.gethostname().startswith("algra"):
            print("Run 'huge_instances' directly on an algra after reserving it with sleep")
            sys.exit(1)
        for f in os.listdir(instances_dir):
            if not f.endswith(".scm.json"):
                continue
            if "Automotive02" not in f:
                continue
            instance_name = f[:-9]
            instance_path = os.path.join(instances_dir, f)
            output_file = os.path.join(output_dir, f"{instance_name}-initial-nosimp.out.json")
            if os.path.exists(output_file):
                continue
            commands.append([command_path, instance_path, 
                             "--dump-initial-phase", output_file,
                             "--dont-simplify", "--print-events"])
            for command in commands:
                run_command(command)
 