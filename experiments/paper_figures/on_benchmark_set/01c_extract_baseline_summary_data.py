if __name__ != "__main__":
    raise ImportError("This script is not meant to be imported.")


import algbench
from pathlib import Path
import lzma
import json
import numpy as np
import math


baseline_data_dir = Path(__file__).parent.parent.parent / "03_baseline" / \
                    "baseline_1h"
summary_data_archive = Path(__file__).parent / "02_baseline_summary_data.json.xz"
EXPECTED_TIME_LIMIT = 3600.0


def read_row(row):
    return {
        "n_samples": len(row["result"]["samples"])\
                     if row["result"]["samples"] else None,
        "solve_time": row["result"]["solve_time"],
        "algorithm": row["parameters"]["args"]["alg_params"]["algorithm"],
        "time_limit": row["parameters"]["args"]["alg_params"]["time_limit"],
        "instance": row["parameters"]["args"]["instance_name"]
    }


data = algbench.read_as_pandas(str(baseline_data_dir), read_row)
algorithm_data = {
    'YASA10': dict(),
    'YASA': dict(),
    'IPOF_FT': dict(),
    'IPOF_CSP': dict(),
    'YASA5': dict(),
    'IL': dict(),
    'IC': dict(),
    'IPOG_CSP': dict(),
    'CH': dict(),
    'IPOG_FT': dict(),
    'YASA3': dict()
}

for row in data.itertuples():
    assert(row.time_limit == EXPECTED_TIME_LIMIT)
    instance = row.instance
    algorithm = row.algorithm
    entry = algorithm_data[algorithm].setdefault(
        instance, {"n_samples": [], "solve_times": [], "runs": 0, "solved": 0})
    entry["runs"] += 1
    if row.n_samples and row.n_samples > 0 and np.isfinite(row.n_samples):
        entry["n_samples"].append(row.n_samples)
        entry["solve_times"].append(row.solve_time)
        entry["solved"] += 1


for algorithm, alg_data in algorithm_data.items():
    for alg_entry in alg_data.values():
        if alg_entry["n_samples"]:
            alg_entry["avg_n_samples"] = np.mean(alg_entry["n_samples"])
            alg_entry["med_n_samples"] = np.median(alg_entry["n_samples"])
            alg_entry["min_n_samples"] = min(alg_entry["n_samples"])
            alg_entry["max_n_samples"] = max(alg_entry["n_samples"])
            alg_entry["avg_time"] = np.mean(alg_entry["solve_times"])
            alg_entry["med_time"] = np.median(alg_entry["solve_times"])
            alg_entry["min_time"] = min(alg_entry["solve_times"])
            alg_entry["max_time"] = max(alg_entry["solve_times"])
        else:
            alg_entry.update({
                "avg_n_samples": None,
                "med_n_samples": None,
                "min_n_samples": None,
                "max_n_samples": None,
                "avg_time": None,
                "med_time": None,
                "min_time": None,
                "max_time": None
            })


with lzma.open(summary_data_archive, "wt") as archive_out:
    json.dump(algorithm_data, archive_out, indent=2)
