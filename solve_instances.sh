#!/bin/bash
# A script to simply solve a set of instances with sammy. You can easily use it to solve your own instances.
# This script processes all .json and .json.xz instance files in the ./instances/ directory
# and generates corresponding output files in the ./results/ directory using a Docker container.
# To build the docker container, run:
# docker build --platform linux/amd64 -t sammy .
#
# The input file need to be converted to a compatible format as the vanilla sammy does not support
# parsing the xml or dimacs instances but relies on other tools for that. A script to convert feature-ide XML-files
# or dimacs files can be found in scripts/instance_to_generic_json.py
#
# sammy requires a gurobi.lic file to be present in the current directory that is compatible with docker, e.g. a WLS license.


# Create results directory if it doesn't exist
mkdir -p ./results/

# check if gurobi.lic exists
if [ ! -f ./gurobi.lic ]; then
    echo "Error: gurobi.lic file not found!"
    echo "Please place your Gurobi license file (gurobi.lic) in the current directory. Note that this needs to be a docker compatible license file, such as WLS."
    exit 1
fi

instance_folder="./instances"
if [ ! -d "$instance_folder" ]; then
    echo "Error: instances directory '$instance_folder' not found!"
    exit 1
fi

sammy_args="--print-events --print-initial-progress --time-limit 60"

# Check if there are any .json or .json.xz files
instance_count=$(find "$instance_folder" -name "*.json" -o -name "*.json.xz" | wc -l)
if [ "$instance_count" -eq 0 ]; then
    echo "Warning: No .json or .json.xz files found in '$instance_folder'!"
    echo "Sammy requires its own json format but you can use the script scripts/instance_to_generic_json.py to convert xml or dimacs files into a compatible json format."
    exit 1
fi

echo "Found $instance_count instance(s) to process."

# Loop through all .json and .json.xz files in the instances directory
for instance_file in "$instance_folder"/*.json "$instance_folder"/*.json.xz; do
    # Skip if no files match the pattern
    [[ -f "$instance_file" ]] || continue
    
    # Get the filename without path
    filename=$(basename "$instance_file")
    
    # Create output filename (handle both .json and .json.xz extensions)
    if [[ "$filename" == *.json.xz ]]; then
        # For .json.xz files, remove .json.xz and add .json
        base_name="${filename%.json.xz}"
        output_file="/results/${base_name}.json"
    else
        # For .json files, just change the path
        output_file="/results/$filename"
    fi
    
    echo "Processing: $filename -> $output_file"
    
    docker run --platform linux/amd64 -it --rm \
        -v ./gurobi.lic:/opt/gurobi/gurobi.lic:ro \
        -v "$instance_folder":/instances/ \
        -v ./results/:/results/ \
        sammy \
        --input-file "/instances/$filename" \
        -o "$output_file" \
        $sammy_args
    
    echo "Completed: $filename"
    echo "---"
done

echo "All instances processed!"