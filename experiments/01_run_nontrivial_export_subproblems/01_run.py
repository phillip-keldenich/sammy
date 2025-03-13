import os
import subprocess

instances_dir = "../../nontrivial_instances"
output_dir = "02_output"
command_path = "../../build/Release/src/sammy_solve"
os.makedirs(output_dir, exist_ok=True)

# avoid OOM conditions
lns_worker_limits = {
	"freetz": 3,
	"FreeBSD-8_0_0": 7,
	"uClinux": 7,
}

for f in os.listdir(instances_dir):
	if not f.endswith(".scm.json"):
		continue
	instance_name = f[:-9]
	instance_path = os.path.join(instances_dir, f)
	subproblem_dir = os.path.join(output_dir, instance_name)
	output_file = subproblem_dir + ".json"
	if os.path.isfile(output_file):
		print("Skipping:", instance_name)
		continue
	print("Running on", instance_name)
	extra = []
	if instance_name in lns_worker_limits:
		extra = ["--max-lns-workers", str(lns_worker_limits[instance_name])]
	subprocess.run([command_path, "--print-events", "--print-global-stats", "-o", output_file,
                    "--report-subproblems-to", subproblem_dir, *extra, instance_path], check=True)

