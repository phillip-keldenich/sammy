#!/bin/sh

zip -9r sammy_package.zip custom_deps \
  experiments/02_run_default_params experiments/03_baselines experiments/05_run_initial_no_simplification \
  experiments/06_identify_nontrivial experiments/07_run_nontrivial_export_subproblems experiments/paper_figures experiments/extract_author_experiment_outputs.py \
  external include src test scripts CMakeLists.txt conanfile.txt \
  full_instances sammy_benchmark_instances \
  README.md experiments/author_experiment_summaries.zip \
  experiments/sammy-data/baselines \
  experiments/sammy-data/analyze \
  experiments/sammy-data/raw_outputs \
  experiments/sammy-data/sammy \
  experiments/sammy-data/samplns \
  experiments/sammy-data/instance_data \
  experiments/sammy-data/requirements.txt \
  -x "*/.git/*" -x "*/.idea/*" -x "*/.vscode/*" \
  -x '**/.*' -x '**/__MACOSX' -x '*/02_output/*' -x "*/__pycache__/*" \
  -x "*egg-info*" -x "*/.ipynb_checkpoints/*" \
  -x "*/05_sammy_summary_data.json.xz" \
  -x "*/04_nosimp_summary_data.json.xz" \
  -x "*/04_summary_data.json.xz" \
  -x "*/05_amended_outcomes.json.xz" \
  -x "*/07_shuffled_nonoptimal.json" \
  -x "*/10_summary_data.json.xz" \
  -x "*/on_benchmark_set/*.json.xz" \
  -x "experiments/sammy-data/analyze/_skbuild/*" \
  -x "experiments/sammy-data/analyze/src/sample_analyzer/verify/*.so" \
  -x '**/.*' \
  -x '**/CMakeUserPresets.json' \
  -x "*/sammy-data/README*"