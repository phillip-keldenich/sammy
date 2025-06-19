import json
import lzma
import bz2
import gzip
import sys
from pathlib import Path


def load(filename: str | Path):
    xz_endings = [".xz", ".lzma"]
    bz_endings = [".bz2", ".bzip2"]
    gz_endings = [".gz", ".gzip"]
    if isinstance(filename, Path):
        flower = filename.name.lower()
    else:
        flower = filename.lower()

    if any(flower.endswith(e) for e in xz_endings):
        with lzma.open(filename, "rt") as f:
            return json.load(f)
    elif any(flower.endswith(e) for e in bz_endings):
        with bz2.open(filename, "rt") as f:
            return json.load(f)
    elif any(flower.endswith(e) for e in gz_endings):
        with gzip.open(filename, "rt") as f:
            return json.load(f)
    else:
        with open(filename, "r") as f:
            return json.load(f)


def is_initial_phase_dump(data):
    return isinstance(data, dict) and "initial_phase_result" in data


def is_solution_output(data):
    return isinstance(data, dict) and "best_solution" in data and\
           isinstance(data["best_solution"], list)


def solution_to_mapping(instance, solution):
    if not solution or not isinstance(solution, list):
        raise ValueError("Best solution is not a list in the initial phase dump.")
    num_concrete = instance["num_concrete_features"]
    num_variables = instance["num_variables"]
    if any(len(sol) != num_variables for sol in solution):
        raise ValueError("Best solution does not match the number of variables in the instance.")
    if len(instance["labels"]) < num_concrete:
        raise ValueError("Number of labels in the instance is not at least the number of concrete features.")
    label_map = instance["labels"]
    mapping = []
    for sol in solution:
        sol_map = {}
        for name, cnf_lit in label_map.items():
            sol_index = cnf_lit - 1 if cnf_lit > 0 else -cnf_lit + 1
            sol_map[name] = sol[sol_index] if cnf_lit > 0 else not sol[sol_index]
        mapping.append(sol_map)
    return mapping


def initial_phase_dump_to_mapping(instance, initial_dump):
    best_solution = initial_dump["initial_phase_result"]["best_solution"]
    return solution_to_mapping(instance, best_solution)


def solution_output_to_mapping(instance, solution_output):
    best_solution = solution_output["best_solution"]
    return solution_to_mapping(instance, best_solution)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: solution_to_feature_name_mapping.py instance solution")
        sys.exit(1)

    instance = Path(sys.argv[1])
    solution = Path(sys.argv[2])

    instance = load(instance)
    solution = load(solution)
    if is_initial_phase_dump(solution):
        mapping = initial_phase_dump_to_mapping(instance, solution)
    elif is_solution_output(solution):
        mapping = solution_output_to_mapping(instance, solution)
    else:
        raise ValueError("Solution file is neither an initial phase dump nor a solution output.")
    print(json.dumps(mapping, indent=2))
