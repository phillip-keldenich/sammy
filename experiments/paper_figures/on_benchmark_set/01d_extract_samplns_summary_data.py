if __name__ != "__main__":
    raise ImportError("This script is not meant to be imported.")


import algbench
from pathlib import Path
import lzma
import json
import numpy as np


baseline_data_dir = Path(__file__).parent.parent.parent / "04_baseline_1h" / \
                    "samplns_1h"
summary_data_archive = Path(__file__).parent / "02_samplns_summary_data.json.xz"
EXPECTED_TIME_LIMIT = 3600.0


def bound_history(result):
    prev_lb = None
    prev_ub = None
    lb_history = []
    ub_history = []
    for iteration_info in result['iteration_info']:
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
    lb_history, ub_history = bound_history(row['result'])
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
    }


data = algbench.read_as_pandas(str(baseline_data_dir), read_row_samplns)
instance_data = {}


for row in data.itertuples():
    assert row.time_limit == EXPECTED_TIME_LIMIT
    assert row.algorithm == "SampLNS + YASA"
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

with lzma.open(summary_data_archive, "wt") as archive_out:
    json.dump(instance_data, archive_out, indent=2)
