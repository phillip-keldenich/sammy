#!/bin/bash
#
# SAMMY Instance Solver Script
# ============================
#
# This script processes all .json and .json.xz instance files in the ./instances/ directory
# and generates corresponding output files in the ./results/ directory using a Docker container.
#
# Prerequisites:
# - Docker container built with: docker build --platform linux/amd64 -t sammy .
# - Gurobi license file (gurobi.lic) compatible with Docker (e.g., WLS license)
# - Instance files in JSON format (use scripts/instance_to_generic_json.py to convert XML/DIMACS)
#
# Usage:
#   ./solve_instances.sh [OPTIONS]
#
# # Use defaults (original behavior)
# ./solve_instances.sh
# # Custom directories
# ./solve_instances.sh -i ./my_instances -o ./my_results
# # Different time limit
# ./solve_instances.sh --time-limit 7200
# # Custom Docker image
# ./solve_instances.sh --docker sammy:v2.0
# # Multiple options
# ./solve_instances.sh -i ./data -t 1800 -a "--verbose --seed 42"
# # Help
# ./solve_instances.sh --help

set -euo pipefail # Exit on error, undefined vars, pipe failures

# Default configuration
INSTANCE_FOLDER="./instances"
RESULTS_FOLDER="./results"
GUROBI_LICENSE="./gurobi.lic"
SAMMY_ARGS="--print-events --time-limit 3600"
DOCKER_IMAGE="sammy"
DOCKER_PLATFORM="linux/amd64"

# Help function
show_help() {
    cat << EOF
Usage: $0 [OPTIONS]

Process .json and .json.xz instance files using SAMMY solver in Docker.

OPTIONS:
    -i, --instances DIR     Instance directory (default: $INSTANCE_FOLDER)
    -o, --output DIR        Results directory (default: $RESULTS_FOLDER)
    -l, --license FILE      Gurobi license file (default: $GUROBI_LICENSE)
    -a, --args ARGS         Additional SAMMY arguments (default: "$SAMMY_ARGS")
    -d, --docker IMAGE      Docker image name (default: $DOCKER_IMAGE)
    -p, --platform PLATFORM Docker platform (default: $DOCKER_PLATFORM)
    -t, --time-limit SEC    Time limit in seconds (default: 3600)
    -h, --help             Show this help message

EXAMPLES:
    $0                                          # Use defaults
    $0 -i ./my_instances -o ./my_results        # Custom directories
    $0 -t 1800 -a "--verbose --seed 42"        # Custom time limit and args
    $0 --instances ./data --time-limit 7200     # Long option format
EOF
}

# Parse command line arguments
parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--instances)
                INSTANCE_FOLDER="$2"
                shift 2
                ;;
            -o|--output)
                RESULTS_FOLDER="$2"
                shift 2
                ;;
            -l|--license)
                GUROBI_LICENSE="$2"
                shift 2
                ;;
            -a|--args)
                SAMMY_ARGS="$2"
                shift 2
                ;;
            -d|--docker)
                DOCKER_IMAGE="$2"
                shift 2
                ;;
            -p|--platform)
                DOCKER_PLATFORM="$2"
                shift 2
                ;;
            -t|--time-limit)
                # Update time limit in SAMMY_ARGS
                SAMMY_ARGS=$(echo "$SAMMY_ARGS" | sed 's/--time-limit [0-9]\+//')
                SAMMY_ARGS="$SAMMY_ARGS --time-limit $2"
                shift 2
                ;;
            -h|--help)
                show_help
                exit 0
                ;;
            *)
                echo "❌ Unknown option: $1" >&2
                show_help >&2
                exit 1
                ;;
        esac
    done
}

# Validate prerequisites
validate_prerequisites() {
    echo "🔍 Validating prerequisites..."

    # Check if Docker can be run
    if ! docker info >/dev/null 2>&1; then
        echo "❌ Error: Docker does not seem to be running or is not installed."
        exit 1
    fi

    # Check if Docker image exists
    if ! docker image inspect "$DOCKER_IMAGE" >/dev/null 2>&1; then
        echo "❌ Error: Docker image '$DOCKER_IMAGE' not found!"
        echo "   Please build it with: docker build --platform $DOCKER_PLATFORM -t $DOCKER_IMAGE ."
        exit 1
    fi

    if [[ ! -f "$GUROBI_LICENSE" ]]; then
        echo "❌ Error: Gurobi license file not found at '$GUROBI_LICENSE'"
        echo "   Please place your Docker-compatible Gurobi license file (WLS) in the specified location."
        exit 1
    fi

    if [[ ! -d "$INSTANCE_FOLDER" ]]; then
        echo "❌ Error: Instance directory '$INSTANCE_FOLDER' not found! Create a folder named '$INSTANCE_FOLDER' and move all the instances you want to solve in it. For testing, we recommend to copy the benchmark instance 'full_instances/berkeleyDB1.scm.json.xz'."
        exit 1
    fi

    echo "✅ Prerequisites validated"
}

# Check for instances
check_instances() {
    local count=$(find -L "$INSTANCE_FOLDER" -name "*.json" -o -name "*.json.xz" 2>/dev/null | wc -l)

    if [[ "$count" -eq 0 ]]; then
        echo "❌ No .json or .json.xz files found in '$INSTANCE_FOLDER'"
        echo "   Use scripts/instance_to_generic_json.py to convert XML/DIMACS files to compatible JSON format."
        exit 1
    fi

    echo "📊 Found $count instance(s) to process"
}

# Process all instances
process_instances() {
    local processed=0
    local any_failed=0

    for instance_file in "$INSTANCE_FOLDER"/*.json "$INSTANCE_FOLDER"/*.json.xz; do
        [[ -f "$instance_file" ]] || continue

        local filename=$(basename "$instance_file")
        local output_file

        # Generate output filename
        if [[ "$filename" == *.json.xz ]]; then
            local base_name="${filename%.json.xz}"
            output_file="$RESULTS_FOLDER/${base_name}.json"
        else
            output_file="$RESULTS_FOLDER/$filename"
        fi

        ((processed++))
        echo "🔄 Processing [$processed]: $filename"
        echo "   Output: $(basename "$output_file")"

        # Run SAMMY in Docker
        if docker run --platform "$DOCKER_PLATFORM" --rm -it \
          -v "$GUROBI_LICENSE":/opt/gurobi/gurobi.lic:ro \
          -v "$INSTANCE_FOLDER":/instances/ \
          -v "$RESULTS_FOLDER":/results/ \
          "$DOCKER_IMAGE" \
          --input-file "/instances/$filename" \
          -o "/results/$(basename "$output_file").xz" \
          $SAMMY_ARGS; then
            echo "✅ Completed: $filename"
        else
            echo "❌ Failed: $filename"
            any_failed=1
        fi

        echo
    done

    echo "📈 Processed $processed instance(s)"
    return $any_failed
}

# Main execution
main() {
    parse_arguments "$@"
    
    echo "🚀 Starting SAMMY batch processing..."
    echo
    
    # Create results directory
    mkdir -p "$RESULTS_FOLDER"
    
    validate_prerequisites
    check_instances
    
    echo
    echo "⚙️  Configuration:"
    echo "   Instances: $INSTANCE_FOLDER"
    echo "   Results:   $RESULTS_FOLDER"
    echo "   License:   $GUROBI_LICENSE"
    echo "   Docker:    $DOCKER_IMAGE ($DOCKER_PLATFORM)"
    echo "   Args:      $SAMMY_ARGS"
    echo
    
    if process_instances; then
        echo "✅ All instances processed successfully!"
    else
        echo "❌ Some instances failed to process."
        exit 1
    fi
}

# Run main function with all arguments
main "$@"
