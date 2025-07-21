from pathlib import Path
import subprocess


if __name__ != "__main__":
    raise ImportError("Not meant to be imported!")

this_file = Path(__file__).absolute()
exp_dir = this_file.parent
summary_file = exp_dir / "author_experiment_summaries.zip"
sammy_data = exp_dir / "sammy-data"
raw_outputs = sammy_data / "raw_outputs"
default_params = raw_outputs / "02_run_default_params.zip"
nosimp = raw_outputs / "05_run_initial_no_simplification.zip"
identify_nontrivial = raw_outputs / "06_identify_nontrivial.zip"
run_subproblems = raw_outputs / "07_run_nontrivial_export_subproblems.zip"

if not summary_file.exists():
    raise FileNotFoundError(f"Expected {summary_file} to exist, but it does not.")

# check the raw files
for file in [default_params, nosimp, identify_nontrivial, run_subproblems]:
    if not file.exists():
        raise FileNotFoundError(f"Expected {file} to exist, but it does not.")

# extract the raw outputs
for file in [default_params, nosimp, identify_nontrivial, run_subproblems]:
    subprocess.run(["unzip", str(file)],
                   cwd=str(exp_dir), check=True)

# extract the summary files
subprocess.run(["unzip", "author_experiment_summaries.zip"], 
               cwd=str(exp_dir), check=True)
