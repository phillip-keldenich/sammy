#!/bin/sh

# Run the lowest-level, fastest reproduction using Docker,
# This runs our algorithm locally, on the machine executing this script,
# in a Docker container, on the benchmark instances from our paper;
# follow the instructions in the README to set up the required Docker container.
# This requires a Gurobi license that works in Docker (e.g., WLS license),
# and the license file should be placed in ./gurobi.lic.
#
# The aim is to reproduce part of Table E.1 in the paper,
# specifically the results for our algorithm 'sammy',
# not the other algorithms.
# Another part that is missing is the comparison to
# sammy without simplification/preprocessing.
#
# Expect errors due to out-of-memory conditions when running the
# largest instances (Automotive02_V[1-4]) on machines with less than 93 GiB
# of free RAM. The other instances should run fine on much weaker machines.

# check if Docker can be run
if ! docker info >/dev/null 2>&1; then
  echo "❌ Error: Docker does not seem to be running or is not installed."
  exit 1
fi

# check if the Gurobi license file is present
if [ ! -f ./gurobi.lic ]; then
    echo "❌ Gurobi license file not found at ./gurobi.lic."
    echo "   Please place your Docker-compatible Gurobi license file (WLS) in the current directory."
    exit 1
fi

# check if benchmark instances are in place (they are part of the repo)
if [ ! -d ./sammy_benchmark_instances ]; then
  echo "❌ Benchmark instance directory not found at ./sammy_benchmark_instances."
  exit 1
fi

# create results directory if it does not exist
mkdir -p ./level1_docker_results 1>/dev/null 2>/dev/null || true

# run each instance 5 times
for run_index in 1 2 3 4 5; do 
  # solve all instances, write results to ./level1_docker_results/tmp
  ./solve_instances.sh -i ./sammy_benchmark_instances -o ./sammy_benchmark_instances/tmp
  # move the results to ./level1_docker_results, renaming them to include the run index
  for output_file in ./sammy_benchmark_instances/tmp/*.json.xz; do
    target_name=$(echo $output_file | grep -Eo "[^/]+$" | sed s/.scm.json.xz//g)
    target_name="./level1_docker_results/${target_name}.scm.out${run_index}.json.xz"
    mv "$output_file" "$target_name"
  done
done

# remove ./instances if it is a symlink
unlink ./instances 1>/dev/null 2>/dev/null || true
python ./level1_generate_table.py ./level1_docker_results ./level1_table_E1.xlsx
echo "✅ Results written to ./level1_table_E1.xlsx"
