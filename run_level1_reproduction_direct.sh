#!/bin/sh

# Run the lowest-level, fastest reproduction directly,
# i.e., without Docker. This runs our algorithm
# locally, on the machine executing this script,
# on the set of benchmark instances used in the paper.
# The aim is to reproduce part of Table E.1 in the paper,
# specifically the results for our algorithm 'sammy',
# not the other algorithms.
# Another part that is missing is the comparison to
# sammy without simplification/preprocessing.
#
# Please follow the direct build instructions prior to running this script,
# including the setup of a Gurobi license (any working non-restricted Gurobi 
# license should work).
#
# Expect errors due to out-of-memory conditions when running the
# largest instances (Automotive02_V[1-4]) on machines with less than 93 GiB
# of free RAM. The other instances should run fine on much weaker machines.

# check if benchmark instances are in place (they are part of the repo)
if [ ! -d ./sammy_benchmark_instances ]; then
  echo "❌ Benchmark instance directory not found at ./sammy_benchmark_instances"
  exit 1
fi

# create results directory if it does not exist
mkdir -p ./level1_direct_results 1>/dev/null 2>/dev/null || true

# process all .json.xz files in ./sammy_benchmark_instances
# and write results to ./level1_direct_results with a time limit of 1h
for input_file in ./sammy_benchmark_instances/*.json.xz; do
  filename=$(basename "$input_file")
  nostem=$(echo "$filename" | sed 's/\.json$//;s/\.json\.xz$//')
  # run each instance 5 times
  for run_index in 1 2 3 4 5; do
    output_file="./level1_direct_results/${nostem}.out${run_index}.json.xz"

    if [ -f "$output_file" ]; then
        echo "✅ Skipping $filename run $run_index, output $output_file already exists."
        continue
    fi

    echo "Processing $filename -> $output_file ..."
    
    if ./build/Release/src/sammy_solve "$input_file" -o "$output_file" --print-events --time-limit 3600; then
        echo "✅ Completed: $filename -> $output_file"

    else
        echo "❌ Failed: $filename -> $output_file"
    fi
  done
  echo
done

# generate table from results; relies on pandas and openpyxl
python ./level1_generate_table.py ./level1_direct_results ./level1_table_E1.xlsx
echo "✅ Results written to ./level1_table_E1.xlsx"
