#!/bin/sh

set -euo pipefail

zip -0r sammy_evaluation_package.zip \
  custom_deps \
  experiments/02_run_default_params/*.py \
  experiments/03_baselines \
  experiments/05_run_initial_no_simplification/*.py \
  experiments/06_identify_nontrivial/*.py \
  experiments/06_identify_nontrivial/*.ipynb \
  experiments/07_run_nontrivial_export_subproblems/*.py \
  experiments/07_run_nontrivial_export_subproblems/*.ipynb \
  experiments/paper_figures/on_benchmark_set/*.py \
  experiments/paper_figures/on_benchmark_set/*.ipynb \
  experiments/sammy-data \
  external/ \
  full_instances/ \
  include/ \
  sammy_benchmark_instances/ \
  scripts/ \
  src/ \
  test/ \
  .clang-format \
  .dockerignore \
  .gitignore \
  .gitmodules \
  CMakeLists.txt \
  conanfile.docker.txt \
  conanfile.txt \
  Dockerfile \
  level1_generate_table.py \
  LICENSE \
  README.md \
  requirements.txt \
  run_level*.sh \
  solve_instances.sh \
  -x 'experiments/03_baselines/requirements.txt' \
  -x "*/.git/*" -x "*/.idea/*" -x "*/.vscode/*" \
  -x '**/__MACOSX' -x "*egg-info*" -x "*/.ipynb_checkpoints/*" \
  -x '**/CMakeUserPresets.json' \
  -x '**/*.so' -x '**/*.pyc' \
  -x '**/_skbuild/*' \
  -x "**/__pycache__/*"

rm -rf sammy_evaluation_package
mkdir -p sammy_evaluation_package
mv sammy_evaluation_package.zip sammy_evaluation_package/sammy_evaluation_package.zip
(cd sammy_evaluation_package && unzip sammy_evaluation_package.zip && rm sammy_evaluation_package.zip)
zip -9r sammy_evaluation_package.zip sammy_evaluation_package
rm -rf sammy_evaluation_package
