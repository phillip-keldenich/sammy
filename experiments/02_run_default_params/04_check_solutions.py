if __name__ != "__main__":
    raise ImportError("This script is not meant to be imported!")


import subprocess
from subprocess import CalledProcessError, TimeoutExpired
from pathlib import Path
import re
import tqdm
import json
import sys


instance_re = re.compile(r"^([^.]+)_repeat[0-9]+\.out\.json\.xz$")
checker_path = Path(__file__).parent.parent.parent / "build" / "Release" /\
               "src" / "check_solution_and_mes"
instances_dir = Path(__file__).parent.parent.parent /\
                "sammy_benchmark_instances"
output_path = Path(__file__).parent / "02_output"
if not output_path.is_dir():
    raise FileNotFoundError(f"Output directory {output_path} does not exist!")
if not instances_dir.is_dir():
    raise FileNotFoundError(f"Instances directory {instances_dir} not found!")
if not checker_path.is_file():
    raise FileNotFoundError(f"Checker executable {checker_path} not found!")


errors = []
for solution_file in tqdm.tqdm(list(output_path.glob("*.out.json.xz"))):
    match = instance_re.match(solution_file.name)
    if not match:
        continue
    instance_name = match.group(1)
    instance_file = instances_dir / f"{instance_name}.scm.json.xz"
    if not instance_file.is_file():
        raise FileNotFoundError(f"Instance file {instance_file} not found!")
    check_file = solution_file.parent / f"{solution_file.name[:-12]}.check.json"
    if check_file.is_file():
        continue  # skip if already checked
    try:
        command = [str(checker_path), str(solution_file), str(instance_file)]
        result = subprocess.run(command, check=True, capture_output=True,
                                text=True)
        with open(check_file, "w") as f:
            f.write(result.stdout)
    except CalledProcessError as e:
        errors.append({
            "instance": instance_name,
            "solution_file": str(solution_file),
            "error": str(e),
            "returncode": e.returncode,
            "output": str(e.output),
            "stderr": str(e.stderr)
        })


with open("errors.json", "w") as f:
    json.dump(errors, f, indent=4)

sys.exit(0 if not errors else 1)
