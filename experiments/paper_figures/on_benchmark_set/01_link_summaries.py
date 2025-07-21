import os
from pathlib import Path

this_path = Path(__file__).absolute().parent
exp_02_dir = this_path.parent.parent / "02_run_default_params"
exp_02_summary_archive = exp_02_dir / "05_sammy_summary_data.json.xz"
exp_05_dir = this_path.parent.parent / "05_run_initial_no_simplification"
exp_05_summary_archive = exp_05_dir / "04_nosimp_summary_data.json.xz"

if not exp_02_summary_archive.exists():
    print(f"Expected summary archive {exp_02_summary_archive} does not exist!")
    print("Please run the experiment or extract the authors' experiment data.")
    raise FileNotFoundError(f"Summary archive {exp_02_summary_archive} not found!")

if not exp_05_summary_archive.exists():
    print(f"Expected summary archive {exp_05_summary_archive} does not exist!")
    print("Please run the experiment or extract the authors' experiment data.")
    raise FileNotFoundError(f"Summary archive {exp_05_summary_archive} not found!")

summary_02_target = this_path / "02_sammy_summary_data.json.xz"
if not summary_02_target.exists():
    os.symlink(exp_02_summary_archive, summary_02_target)

summary_05_target = this_path / "02_nosimp_summary_data.json.xz"
if not summary_05_target.exists():
    os.symlink(exp_05_summary_archive, summary_05_target)
