import os
import subprocess
import sys
import slurminade
import traceback
import random


slurminade.update_default_configuration(
    partition="alg", exclusive=True, constraint="alggen05"
)

@slurminade.slurmify()
def run_command(command):
    subprocess.run(command, capture_output=True, check=True)


LIMITS = {
    "Automotive02": 0,
}

if __name__ == "__main__":
    instances_dir = "../../sammy_benchmark_instances"
    output_dir = "02_output"
    command_path = "../../build/Release/src/sammy_solve"
    os.makedirs(output_dir, exist_ok=True)
    commands = []

    for f in os.listdir(instances_dir):
        if not f.endswith(".scm.json"):
            continue
        instance_name = f[:-9]
        instance_path = os.path.join(instances_dir, f)
        output_file = os.path.join(output_dir, f"{instance_name}.out.json")
        if os.path.exists(output_file):
            continue
        extra = []
        for key, limit in LIMITS.items():
            if key in f:
                if limit > 0:
                    extra += ["--max-lns-workers", str(limit)]
                else:
                    extra += ["--only-initial-phase"]
                break
        commands.append([command_path, instance_path, "-o", output_file, *extra])

    random.shuffle(commands)
    with slurminade.JobBundling(max_size=1):
        for command in commands:
            run_command.distribute(command)

