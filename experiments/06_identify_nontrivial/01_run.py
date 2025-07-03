from pathlib import Path
import slurminade
import subprocess
from subprocess import CalledProcessError
import lzma
import json
import os
import random


slurminade.update_default_configuration(
    partition="alg",
    exclusive=True,
    constraint="alggen05",
    mail_type="FAIL"
)  # global options for slurm

slurminade.set_dispatch_limit(300)

@slurminade.slurmify()
def run_commands(commands, output_file_key):
    for command in commands:
        output_file = command[command.index(output_file_key)+1]
        try:
            subprocess.run(command, capture_output=True, check=True, text=True)
            subprocess.run(['xz', '-9', '-e', output_file, '-T', '24'])
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
    instances_dir = Path(__file__).parent.parent.parent / 'full_instances'
    benchmark_dir = Path(__file__).parent.parent.parent / 'sammy_benchmark_instances'
    executable_path = Path(__file__).parent.parent.parent / 'build' / 'Release' / 'src' /  'sammy_solve'
    executable_path = executable_path.absolute()
    if not executable_path.exists():
        raise RuntimeError("Could not find executable!")
    instances = []
    benchmark_names = {f.name for f in benchmark_dir.glob('*.scm.json.xz')}
    for file in instances_dir.glob('*.scm.json.xz'):
        if file.name not in benchmark_names:
            instances.append(file.absolute())
    output_dir = Path(__file__).parent / "02_output"
    os.makedirs(output_dir, exist_ok=True)
    commands_list = []
    for instance in instances:
        name = instance.name
        output = (output_dir / name.replace(".scm.json.xz", ".out.json")).absolute()
        if os.path.exists(str(output) + ".xz"):
            print("SKIPPING (existing output file):", instance.name)
            continue
        command = [str(executable_path), str(instance), '--print-events',
                   '--max-lns-workers=2', '-o', str(output)]
        commands_list.append(command)
    random.shuffle(commands_list)
    with slurminade.JobBundling(max_size=5):
        while commands_list:
            next_batch = []
            while commands_list and len(next_batch) < 5:
                next_batch.append(commands_list.pop())
            run_commands.distribute(next_batch, '-o')

