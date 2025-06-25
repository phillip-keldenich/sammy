from pathlib import Path
import re
import lzma
import json
import numpy as np
import tqdm
import sys


def scan_events(events):
    extra_data = {"lb_history": [], "ub_history": [], "lb_is_optimal": False}
    for e in events:
        t = e["type"]
        if t == "INPUT_READ":
            for k in ("num_clauses", "num_concrete", "num_vars"):
                extra_data[k] = e[k]
        elif t == "INITIAL_PHASE_DONE":
            extra_data["initial_phase_lb"] = e["lb"]
            extra_data["initial_phase_ub"] = e["ub"]
            extra_data["initial_phase_time"] = e["time"]
            extra_data["initial_num_interactions"] = e["num_interactions"]
            extra_data["reduced_num_interactions"] = e["num_interactions"]
        elif t == "DONE_SIMPLIFICATION":
            for sstat, value in e["simplification"]["simplification_stats"].items():
                extra_data[f"simplification_{sstat}"] = value
        elif t == "FIRST_COLORING_DONE":
            extra_data["time_to_first_solution"] = e["time"]
        elif t == "DONE_IMPLIED_VERTEX_ELIMINATION":
            extra_data["reduced_num_interactions"] = e["reduced_size"]
        elif t in {"INITIAL_ITERATION_DONE","STRENGTHENED_CONSTRAINTS",} or "SOLVE_FULL_RELAX" in t:
            continue
        elif t == "LB_IMPROVED" or t == "IMPROVED_LB":
            source = None
            key = None
            if "lb" in e:
                key = "lb"
            else:
                key = "num_interactions"
            if "iteration" in e:
                source = "initial_phase"
            elif 'method' in e:
                source = e['method']
            elif 'source' in e:
                source = e['source']
            if source is None:
                raise RuntimeError("NO SOURCE FOR " + str(e))
            lb = e[key]
            if not extra_data["lb_history"] or extra_data["lb_history"][-1]["lb"] < lb:
                extra_data["lb_history"].append({"time": e["time"], "lb": lb, "source": source})
        elif t == "UB_IMPROVED" or t == "IMPROVED_UB" or t == "IMPROVED_SOLUTION":
            source = None
            key = "num_configs"
            if 'iteration' in e:
                source = "initial_phase"
            elif 'source' in e:
                source = e["source"]
            if key not in e:
                key = "ub"
            if key not in e:
                key = "size"
            if key not in e:
                raise RuntimeError("NO UB FOR "+str(e))
            if source is None:
                raise RuntimeError("NO SOURCE FOR " + str(e))
            extra_data["ub_history"].append({"time": e["time"], "ub": e[key], "source": source})
        elif t == "PRICING_FOUND_OPTIMALITY":
            extra_data["lb_is_optimal"] = True
    return extra_data


def store_list_measure(out_data, measure_name, measures):
    out_data[f"{measure_name}s"] = measures
    out_data[f"min_{measure_name}"] = min(measures)
    out_data[f"max_{measure_name}"] = max(measures)
    out_data[f"avg_{measure_name}"] = np.mean(measures)
    out_data[f"med_{measure_name}"] = np.median(measures)


if __name__ != "__main__":
    raise ImportError("This script is not meant to be imported.")


sammy_data_dir = Path(__file__).parent.parent.parent / \
                 "02_run_default_params" / "02_output"
summary_data_archive = Path(__file__).parent / "02_sammy_summary_data.json.xz"

if summary_data_archive.exists():
    print("Summary data archive already exists - skipping!")
    sys.exit(0)

TOTAL_INSTANCES = 55
REPEATS_PER_INSTANCE = 5

instance_name_re = re.compile(r"^([^.]+)_repeat([0-9]+).out.json.xz$")
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
    feature_count_reduced = []
    ubs = []
    lbs = []
    solve_times = []
    initial_lbs = []
    initial_ubs = []
    initial_times = []
    first_sol_times = []
    optimalities = []
    lb_optimalities = []
    simp_measures = {}
    instance_size_data = {}
    for file, repeat in repeats:
        with lzma.open(file, "rt") as input:
            js_data = json.load(input)
        extra_data = scan_events(js_data["events"])
        if js_data["lb"] == js_data["ub"]:
            extra_data["lb_is_optimal"] = True
        if feature_count_nonreduced is None:
            feature_count_nonreduced = extra_data["initial_num_interactions"]
        else:
            if feature_count_nonreduced != extra_data["initial_num_interactions"]:
                raise RuntimeError(f"Mismatch between interaction counts: {file}")
        ubs.append(js_data["ub"])
        lbs.append(js_data["lb"])
        initial_lbs.append(extra_data["initial_phase_lb"])
        initial_ubs.append(extra_data["initial_phase_ub"])
        initial_times.append(extra_data["initial_phase_time"])
        first_sol_times.append(extra_data["time_to_first_solution"])
        solve_times.append(js_data["events"][-1]["time"])
        optimalities.append(js_data["lb"] >= js_data["ub"])
        lb_optimalities.append(optimalities[-1] or extra_data["lb_is_optimal"])
        feature_count_reduced.append(extra_data["reduced_num_interactions"])
        if not simp_measures:
            for m_name, measure in extra_data.items():
                if m_name.startswith("simplification_"):
                    simp_measures[m_name] = measure
        if not instance_size_data:
            instance_size_data["num_clauses"] = extra_data["num_clauses"]
            instance_size_data["num_variables"] = extra_data["num_vars"]
            instance_size_data["num_concrete"] = extra_data["num_concrete"]
    best_lb = max(lbs)
    best_ub = min(ubs)
    if best_ub < best_lb:
        raise RuntimeError("LB/UB violation")
    out_data = {}
    instance_data[instance] = out_data
    out_data["one_optimal"] = (best_lb == best_ub)
    out_data["all_optimal"] = all(optimalities)
    out_data["optimalities"] = optimalities
    store_list_measure(out_data, "ub", ubs)
    store_list_measure(out_data, "lb", lbs)
    store_list_measure(out_data, "initial_lb", initial_lbs)
    store_list_measure(out_data, "initial_ub", initial_ubs)
    store_list_measure(out_data, "solve_time", solve_times)
    store_list_measure(out_data, "initial_time", initial_times)
    store_list_measure(out_data, "first_sol_time", first_sol_times)
    out_data["interactions_nonreduced"] = feature_count_nonreduced
    store_list_measure(out_data, "interactions_reduced", feature_count_reduced)
    out_data.update(simp_measures)
    out_data.update(instance_size_data)

with lzma.open(summary_data_archive, "wt") as f:
    json.dump(instance_data, f)
