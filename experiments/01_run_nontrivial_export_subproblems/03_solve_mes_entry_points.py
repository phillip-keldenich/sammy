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
def run_commands(commands):
    exc = None
    for command in commands:
        try:
            subprocess.run(command, capture_output=True, check=True)
        except subprocess.CalledProcessError as e:
            print("Error in called process!")
            traceback.print_exc()
            exc = e
    if exc is not None:
        raise exc


if __name__ == "__main__":
    output_dir = "02_output"
    solve_binary = "../../build/Release/src/solve_subproblems_mes_entry_point"
    time_limit = 10 * 60

    commands = []

    for f in os.listdir(output_dir):
        path = os.path.join(output_dir, f)
        if not os.path.isdir(path):
            continue
        instance_name = f
        for entry_point in os.listdir(path):
            if not entry_point.startswith("mes_entry_points-"):
                continue
            mep_name = entry_point[17:]
            mep_path = os.path.join(path, entry_point)
            out_path = os.path.join(path, f"solutions-{mep_name}")
            if not os.path.isfile(out_path):
                commands.append([solve_binary, "--entry-point-file", mep_path, "-o", out_path, "--solve-time-limit", str(time_limit)])

    random.shuffle(commands)
    with slurminade.JobBundling(max_size=10):
        current_subbatch = []
        for command in commands:
            current_subbatch.append(command)
            if len(current_subbatch) >= 100:
                run_commands.distribute(current_subbatch)
                current_subbatch = []
        if current_subbatch:
            run_commands.distribute(current_subbatch)

