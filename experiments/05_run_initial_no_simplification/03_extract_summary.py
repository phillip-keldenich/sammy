from pathlib import Path
import re
import lzma
import json
import numpy as np
import tqdm
import sys


def scan_events(events):
    # we need: time to first solution,
    # time for initial phase,
    # number of initial phase iterations,
    # time to 10 iterations,
    # number of interactions
    extra_data = {}
    for e in events:
        t = e['type']
        if t == 'INITIAL_PHASE_DONE':
            extra_data["initial_phase_lb"] = e["lb"]
            extra_data["initial_phase_ub"] = e["ub"]
            extra_data["initial_phase_time"] = e["time"]
            extra_data["initial_num_interactions"] = e["num_interactions"]
        elif t == "FIRST_COLORING_DONE":
            extra_data["time_to_first_solution"] = e["time"]
    return extra_data


def extend_list_measure(dict, measure_name):
    measures = dict[measure_name + "s"]
    dict[measure_name + "_mean"] = np.mean(measures)
    dict[measure_name + "_med"] = np.median(measures)
    dict[measure_name + "_max"] = max(measures)
    dict[measure_name + "_min"] = min(measures)


if __name__ != "__main__":
    raise ImportError("This script is not meant to be imported.")


this_path = Path(__file__).absolute()
sammy_data_dir = this_path.parent / "02_output"
summary_data_archive = this_path.parent / "04_nosimp_summary_data.json.xz"

if summary_data_archive.exists():
    print("Summary data archive already exists - skipping!")
    sys.exit(0)

TOTAL_INSTANCES = 55
REPEATS_PER_INSTANCE = 5

instance_name_re = re.compile(r"^([^.]+)-nosimp_repeat([0-9]+).out.json.xz$")
instances_and_files = {}
for file in sammy_data_dir.glob("*.out.json.xz"):
    name = file.name
    m = instance_name_re.match(name)
    if not m:
        raise RuntimeError(f"Unexpected file name: {name}")
    instance_name = m[1]
    repeat = int(m[2])
    instances_and_files.setdefault(instance_name, []).append((file, repeat))

if len(instances_and_files) != TOTAL_INSTANCES:
    raise RuntimeError(f"Expected exactly {TOTAL_INSTANCES} different instances!")

instance_data = {}
for instance, repeats in tqdm.tqdm(instances_and_files.items()):
    if len(repeats) != REPEATS_PER_INSTANCE:
        raise RuntimeError(f"Instance {instance} has {len(repeats)} repeats; expected {REPEATS_PER_INSTANCE}")
    feature_count_nonreduced = None
    initial_lbs = []
    initial_ubs = []
    initial_times = []
    first_sol_times = []
    for file, repeat in repeats:
        with lzma.open(file, "rt") as input:
            js_data = json.load(input)
        extra_data = scan_events(js_data["events"])
        initial_lbs.append(extra_data["initial_phase_lb"])
        initial_ubs.append(extra_data["initial_phase_ub"])
        initial_times.append(extra_data["initial_phase_time"])
        first_sol_times.append(extra_data["time_to_first_solution"])
        if feature_count_nonreduced is None:
            feature_count_nonreduced = extra_data["initial_num_interactions"]
        else:
            if feature_count_nonreduced != extra_data["initial_num_interactions"]:
                raise RuntimeError(f"Mismatch between interaction counts: {file}")
    instance_data[instance] = {
        "initial_lbs": initial_lbs,
        "initial_ubs": initial_ubs,
        "initial_times": initial_times,
        "first_sol_times": first_sol_times,
        "feature_count_nonreduced": feature_count_nonreduced
    }
    extend_list_measure(instance_data[instance], "initial_lb")
    extend_list_measure(instance_data[instance], "initial_ub")
    extend_list_measure(instance_data[instance], "initial_time")
    extend_list_measure(instance_data[instance], "first_sol_time")

with lzma.open(summary_data_archive, "wt") as output:
    json.dump(instance_data, output, indent=2)
