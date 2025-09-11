import os
from pathlib import Path

this_path = Path(__file__).absolute().parent
exp_02_dir = this_path.parent.parent / "02_run_default_params"
exp_02_summary_archive = exp_02_dir / "05_sammy_summary_data.json.xz"

exp_03_baselines_archive = this_path.parent.parent / "sammy-data" / "baselines" / "02_baseline_summary_data.json.xz"
exp_03_samplns_archive = this_path.parent.parent / "sammy-data" / "samplns" / "02_samplns_summary_data.json.xz"

exp_05_dir = this_path.parent.parent / "05_run_initial_no_simplification"
exp_05_summary_archive = exp_05_dir / "04_nosimp_summary_data.json.xz"

for file in [exp_02_summary_archive, exp_03_baselines_archive, 
             exp_03_samplns_archive, exp_05_summary_archive]:
    if not file.exists():
        print(f"Expected summary archive {file} does not exist!")
        print("Please run the experiment or extract the authors' experiment data.")
        raise FileNotFoundError(f"Summary archive {file} not found!")

summary_02_target = this_path / "02_sammy_summary_data.json.xz"
if not summary_02_target.exists():
    os.symlink(exp_02_summary_archive, summary_02_target)

summary_03_samplns_target = this_path / "02_samplns_summary_data.json.xz"
if not summary_03_samplns_target.exists():
    os.symlink(exp_03_samplns_archive, summary_03_samplns_target)

summary_03_baselines_target = this_path / "02_baseline_summary_data.json.xz"
if not summary_03_baselines_target.exists():
    os.symlink(exp_03_baselines_archive, summary_03_baselines_target)

summary_05_target = this_path / "02_nosimp_summary_data.json.xz"
if not summary_05_target.exists():
    os.symlink(exp_05_summary_archive, summary_05_target)
