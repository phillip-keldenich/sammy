import os
import subprocess
import sys

output_dir = "02_output"
generate_binary = "../../build/Release/src/generate_subproblem_mes_entry_points"
time_limit = 5 * 60

for f in os.listdir(output_dir):
    path = os.path.join(output_dir, f)
    if not os.path.isdir(path):
        continue
    formula_file = os.path.join(path, "universe_and_clauses.json.xz")
    if not os.path.isfile(formula_file):
        continue
    for subproblem_file in os.listdir(path):
        if not subproblem_file.endswith(".json.xz") or \
           not subproblem_file.startswith("subproblem-"):
            continue
        output_file = "mes_entry_points" + subproblem_file[10:]
        subproblem_path = os.path.join(path, subproblem_file)
        output_path = os.path.join(path, output_file)
        if os.path.isfile(output_path):
            print("Skipping", subproblem_path)
            continue
        print("Generating for", subproblem_path)
        subprocess.run([generate_binary, "--formula", formula_file,
                        "--subproblem", subproblem_path,
                        "--output", output_path, "--time-limit", str(time_limit)], check=True)
