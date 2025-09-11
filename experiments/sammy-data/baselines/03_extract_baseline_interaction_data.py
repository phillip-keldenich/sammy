from sample_analyzer.verify import count_interactions

if __name__ != "__main__":
    raise ImportError("This script is not meant to be imported.")

import algbench
from pathlib import Path
import lzma
import json
import numpy as np

baseline_data_dir = Path(__file__).parent / "baseline_1h"
yasa_1h_data_dir = Path(__file__).parent / "yasa_1h"
yasa_fixed_m_data_dir = Path(__file__).parent / "yasa_fixed_m"

summary_data_archive = Path(__file__).parent / "03a_baseline_summary_data_with_interactions.json.xz"
EXPECTED_TIME_LIMIT = 3600.0

CONCRETE_FEATURES_FILE = Path(__file__).parent / "concrete_features_for_instance.json"

with open(CONCRETE_FEATURES_FILE, "r") as f:
    concrete_features_for_instance = json.load(f)


def read_row(row):
    return {
        "n_samples": len(row["result"]["samples"]) \
            if row["result"]["samples"] else None,
        "sample": row["result"]["samples"],
        "solve_time": row["result"]["solve_time"],
        "algorithm": row["parameters"]["args"]["alg_params"]["algorithm"],
        "time_limit": row["parameters"]["args"]["alg_params"]["time_limit"],
        "instance": row["parameters"]["args"]["instance_name"],
        "n_interactions": count_interactions(
            row["result"]["samples"],
            concrete_features_for_instance[row["parameters"]["args"]["instance_name"]]
        ) if row["result"]["samples"] else None,
        "m": row["parameters"]["args"]["alg_params"]["m"] if "m" in row["parameters"]["args"]["alg_params"] else None,
    }


algorithm_data = {
    'YASA (1h)': dict(),
    'YASA': dict(),
    'YASA3': dict(),
    'YASA5': dict(),
    'YASA10': dict(),
    'IL': dict(),
    'IC': dict(),
    'IPOG_CSP': dict(),
    'CH': dict(),
    'IPOF_FT': dict(),
    'IPOF_CSP': dict(),
    'IPOG_FT': dict()
}


def add_data_from_dataframe(data, desired_interactions=None):
    for row in data.itertuples():
        assert (row.time_limit == EXPECTED_TIME_LIMIT)
        instance = row.instance
        algorithm = row.algorithm

        if (desired_interactions is not None and instance in desired_interactions
                and np.isfinite(row.n_interactions) and
                row.n_interactions != desired_interactions[instance]):
            print(f"Warning: Instance {instance} has {row.n_interactions} "
                  f"interactions for algorithm {algorithm}, expected {desired_interactions[instance]}.")
            continue

        entry = algorithm_data[algorithm].setdefault(
            instance, {"n_samples": [], "solve_times": [], "runs": 0, "solved": 0, "n_interactions": []})
        entry["runs"] += 1
        if row.n_samples and row.n_samples > 0 and np.isfinite(row.n_samples):
            entry["n_samples"].append(row.n_samples)
            entry["solve_times"].append(row.solve_time)
            entry["n_interactions"].append(row.n_interactions)
            entry["solved"] += 1


yasa_1h_data = algbench.read_as_pandas(str(yasa_1h_data_dir), read_row)
yasa_1h_data["algorithm"] = "YASA (1h)"
print("Found", len(yasa_1h_data), "YASA (1h) entries")

add_data_from_dataframe(yasa_1h_data)

desired_interactions = {
    instance: int(algorithm_data["YASA (1h)"][instance]["n_interactions"][0])
    for instance in algorithm_data["YASA (1h)"]
    if algorithm_data["YASA (1h)"][instance]["n_interactions"]
}

yasa_fixed_m_data = algbench.read_as_pandas(str(yasa_fixed_m_data_dir), read_row)
yasa_fixed_m_data["algorithm"] = yasa_fixed_m_data["m"].apply(lambda m: f"YASA{m}" if m > 1 else "YASA")
print("Found", len(yasa_fixed_m_data), "YASA fixed m entries")

add_data_from_dataframe(yasa_fixed_m_data)

baseline_data = algbench.read_as_pandas(str(baseline_data_dir), read_row)
print("Found", len(baseline_data), "baseline entries")
print("Removing YASA runs (as they are replaced by newer YASA runs)")
baseline_data = baseline_data[~baseline_data['algorithm'].isin(["YASA", "YASA3", "YASA5", "YASA10"])]
print("Pruned baseline data has", len(baseline_data), "entries")

add_data_from_dataframe(baseline_data, desired_interactions=desired_interactions)

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
