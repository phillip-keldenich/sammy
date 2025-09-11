import lzma
import json
import algbench
import numpy as np
from pathlib import Path

instance_data = {}

def bound_history(iteration_info):
    prev_lb = None
    prev_ub = None
    lb_history = []
    ub_history = []
    for iteration_info in iteration_info:
        lb = iteration_info['lb']
        ub = iteration_info['ub']
        time = iteration_info['time']
        if lb != prev_lb:
            prev_lb = lb
            lb_history.append({'lb': lb, 'time': time})
        if ub != prev_ub:
            prev_ub = ub
            ub_history.append({'ub': ub, 'time': time})
    return lb_history, ub_history

def read_row_samplns(row):
    lb_history, ub_history = bound_history(row['result']["iteration_info"])
    return {
        "lower_bound": row["result"]["lower_bound"],
        "upper_bound": row["result"]["upper_bound"],
        "time_used_by_yasa": row["result"]["time_used_by_yasa"],
        "timelimit_for_samplns": row["result"]["timelimit_for_samplns"],
        "runtime": row["runtime"],
        "algorithm": "SampLNS + YASA",
        "time_limit": row["parameters"]["args"]["time_limit"],
        "instance": row["parameters"]["args"]["instance_name"],
        "lb_history": lb_history,
        "ub_history": ub_history,
        "sample": row["result"]["solution"],
    }


data = algbench.read_as_pandas("./samplns_1h", read_row_samplns)

samples_for_instances = {
    instance: list() for instance in data["instance"].unique()
}

for row in data.itertuples():
    assert row.time_limit == 3600
    assert row.algorithm == "SampLNS + YASA"
    samples_for_instances[row.instance].append(row.sample)
    entry = instance_data.setdefault(row.instance, {
        "lower_bounds": [],
        "upper_bounds": [],
        "lb_histories": [],
        "ub_histories": [],
        "solve_times": [],
        "runs": 0,
        "solved": 0
    })
    entry["runs"] += 1
    if row.upper_bound and row.upper_bound > 0 and np.isfinite(row.upper_bound):
        entry["solved"] += 1
        if row.lower_bound and row.lower_bound > 0 and np.isfinite(row.lower_bound):
            entry["lower_bounds"].append(row.lower_bound)
        else:
            entry["lower_bounds"].append(0.0)
        entry["upper_bounds"].append(row.upper_bound)
        entry["solve_times"].append(row.runtime)
        entry["ub_histories"].append(row.ub_history)
        entry["lb_histories"].append(row.lb_history)

with lzma.open(Path(__file__).parent.absolute() / "02_samplns_summary_data.json.xz", "wt") as archive_out:
    json.dump(instance_data, archive_out, indent=2)

assert len(instance_data) == 51

with lzma.open(Path(__file__).parent.absolute() / "02a_samplns_samples.json.xz", "wt") as archive_out:
    json.dump(samples_for_instances, archive_out)