#!/bin/sh

zip -9r sammy_package.zip custom_deps \
  experiments/02_run_default experiments/03_baselines experiments/05_run_initial_no_simplification experiments/06_identify_nontrivial experiments/07_run_nontrivial_export_subproblems experiments/paper_figures \
  external full_instances sammy_benchmark_instances CMakeLists.txt conanfile.txt README.md include src test scripts -x "*/.git/*" -x "*/.idea/*" -x "*/.vscode/*" -x '**/.*' -x '**/__MACOSX' -x '*/02_output/*' -x "*/__pycache__/*" -x "*egg-info*" -x "*/.ipynb_checkpoints/*"
