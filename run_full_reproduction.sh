#!/bin/sh

die() {
  echo "$*" 1>&2
  exit 1
}

python --version 1>/dev/null 2>/dev/null || die 'python not found'
jupyter --version 1>/dev/null 2>/dev/null || die 'jupyter not found'

# we assume gurobi license is setup; otherwise, this will detect that error
echo "Running Sammy on test instance to check for licenses and other issues"
./build/Release/src/sammy_solve sammy_benchmark_instances/soletta_2017-03-09_21-02-40.scm.json.xz || die 'failed to run sammy on a test instance'

# now the fun part: run the experiments; that will take quite long (expect weeks?);
# some instances require ~93 GiB of RAM
echo "Running Sammy on benchmark instances"
(cd experiments/02_run_default_params && python 01a_run_small.py) || die 'error running small instance experiments'
(cd experiments/02_run_default_params && python 01b_run_large.py) || die 'error running large instance experiments'
(cd experiments/02_run_default_params && python 03_check_solutions.py) || die 'error checking solutions'
(cd experiments/02_run_default_params && python 04_extract_summary.py) || die 'error extracting summary data'

echo "Running baseline algorithms"
(cd experiments/03_baselines && python 01_run_featjar_baselines.py 1>/dev/null) || die 'error running featjar baselines'
echo "Running YASA baselines"
(cd experiments/03_baselines && python 02_run_yasa.py 1>/dev/null) || die 'error running YASA baselines'
echo "Running SampLNS baseline"
(cd experiments/03_baselines && python 03_run_samplns.py 1>/dev/null) || die 'error running SampLNS baselines'
echo "Extracting benchmark summary data"
# TODO: these scripts are currently only in the IBR sammy-data repository...
(cd experiments/sammy-data/baselines && python 01_extract_baseline_summary_data.py) || die 'error extracting baseline data'
(cd experiments/sammy-data/samplns && python 01_extract_samplns_summary_data.py) || die 'error extracting SampLNS data'

echo "Running sammy without initial simplification"
(cd experiments/05_run_initial_no_simplification && python 01a_run_small.py) || die 'error running small instances without simplification'
(cd experiments/05_run_initial_no_simplification && python 01b_run_large.py) || die 'error running large instances without simplification'
(cd experiments/05_run_initial_no_simplification && python 03_extract_summary.py) || die 'error extracting summary data without simplification'

echo "Running sammy on full instance set"
(cd experiments/06_identify_nontrivial && python 01_run.py) || die 'error running identify nontrivial on full instance set'
(cd experiments/06_identify_nontrivial && 03_extract_summary.py) || die 'error extracting data for identify nontrivial on full instance set'
(
  cd experiments/06_identify_nontrivial
  jupyter nbconvert --to notebook --inplace --execute 05_look_at_results.ipynb --ExecutePreprocessor.timeout=3600
)

echo "Running different subproblem solvers"
(cd experiments/07_run_nontrivial_export_subproblems && python 01_run.py) || die 'error running nontrivial instances with subproblem export'
(cd experiments/07_run_nontrivial_export_subproblems && python 03_amend_subproblems.py) || die 'error amending nontrivial instances with subproblem export'
(cd experiments/07_run_nontrivial_export_subproblems && python 04_classify_amended_outcomes.py) || die 'error amending nontrivial instances with subproblem export'
(cd experiments/07_run_nontrivial_export_subproblems && python 06_sample_nonoptimal.py) || die 'error sampling non-optimal subproblems'
(cd experiments/07_run_nontrivial_export_subproblems && python 08_run_subproblem_solvers.py) || die 'error running subproblem solvers'
(cd experiments/07_run_nontrivial_export_subproblems && python 09_create_summary.py) || die 'error creating summary files'
(
  cd experiments/07_run_nontrivial_export_subproblems
  jupyter nbconvert --to notebook --inplace --execute 11_look_at_results.ipynb --ExecutePreprocessor.timeout=3600
)

echo "Running notebooks to produce plots"
python experiments/paper_figures/on_benchmark_set/01_link_summaries.py
(
  cd experiments/paper_figures/on_benchmark_set
  jupyter nbconvert --to notebook --inplace --execute 03_plots.ipynb --ExecutePreprocessor.timeout=3600
)
