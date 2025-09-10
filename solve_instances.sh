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
#   ./solve_instances.sh
#

set -euo pipefail # Exit on error, undefined vars, pipe failures

# Configuration
INSTANCE_FOLDER="./instances"
RESULTS_FOLDER="./results"
GUROBI_LICENSE="./gurobi.lic"
SAMMY_ARGS="--print-events --time-limit 60"

# Create results directory
mkdir -p "$RESULTS_FOLDER"

# Validate prerequisites
validate_prerequisites() {
  echo "Validating prerequisites..."

  # check if Docker can be run
  if ! docker info >/dev/null 2>&1; then
    echo "❌ Error: Docker does not seem to be running or is not installed."
    exit 1
  fi

  if [[ ! -f "$GUROBI_LICENSE" ]]; then
    echo "❌ Error: Gurobi license file not found at '$GUROBI_LICENSE'"
    echo "   Please place your Docker-compatible Gurobi license file (WLS) in the current directory."
    exit 1
  fi

  if [[ ! -d "$INSTANCE_FOLDER" ]]; then
    echo "❌ Error: Instance directory '$INSTANCE_FOLDER' not found!"
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

# Main execution
main() {
  echo "🚀 Starting SAMMY batch processing..."
  echo

  validate_prerequisites
  check_instances

  echo
  echo "⚙️  Configuration:"
  echo "   Instances: $INSTANCE_FOLDER"
  echo "   Results:   $RESULTS_FOLDER"
  echo "   Args:      $SAMMY_ARGS"
  echo

  if process_instances; then
    echo "✅ All instances processed successfully!"
  else
    echo "❌ Some instances failed to process."
    exit 1
  fi
}

# Process all instances
process_instances() {
  local processed=0

  for instance_file in "$INSTANCE_FOLDER"/*.json "$INSTANCE_FOLDER"/*.json.xz; do
    [[ -f "$instance_file" ]] || continue

    local filename=$(basename "$instance_file")
    local output_file
    local any_failed=0

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
    if docker run --platform linux/amd64 --rm -it \
      -v "$GUROBI_LICENSE":/opt/gurobi/gurobi.lic:ro \
      -v "$INSTANCE_FOLDER":/instances/ \
      -v "$RESULTS_FOLDER":/results/ \
      sammy \
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

# Run main function
main "$@"
