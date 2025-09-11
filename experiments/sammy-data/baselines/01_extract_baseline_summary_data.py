if __name__ != "__main__":
    raise ImportError("This script is not meant to be imported.")

import algbench
from pathlib import Path
import lzma
import json
import numpy as np

# Contains FIDE YASA runs that should be removed. Contains IPOG runs that have issues.
baseline_data_dir = Path(__file__).parent / "baseline_1h"

# we have two YASA data directories as the experiments were run separately. In the experiment script everything is
# merged into a single directory, but here we keep them separate for clarity.
yasa_1h_data_dir = Path(__file__).parent / "yasa_1h"
yasa_fixed_m_data_dir = Path(__file__).parent / "yasa_fixed_m"

summary_data_archive = Path(__file__).parent / "02_baseline_summary_data.json.xz"
EXPECTED_TIME_LIMIT = 3600.0


def read_row_with_excluded_algorithms(excluded_algorithms):
    def read_row(row):
        if row["parameters"]["args"]["alg_params"]["algorithm"] in excluded_algorithms:
            return None

        return {
            "n_samples": len(row["result"]["samples"]) \
                if row["result"]["samples"] else None,
            "sample": row["result"]["samples"],
            "solve_time": row["result"]["solve_time"],
            "algorithm": row["parameters"]["args"]["alg_params"]["algorithm"],
            "time_limit": row["parameters"]["args"]["alg_params"]["time_limit"],
            "instance": row["parameters"]["args"]["instance_name"],
            "m": row["parameters"]["args"]["alg_params"]["m"] if "m" in row["parameters"]["args"][
                "alg_params"] else None,
        }

    return read_row


algorithm_data = {
    'YASA (1h)': dict(),
    'YASA': dict(),
    'YASA3': dict(),
    'YASA5': dict(),
    'YASA10': dict(),
    'IL': dict(),
    'IC': dict(),
    'CH': dict(),
}


def add_data_from_dataframe(data):
    for row in data.itertuples():
        assert (row.time_limit == EXPECTED_TIME_LIMIT)
        instance = row.instance
        algorithm = row.algorithm

        entry = algorithm_data[algorithm].setdefault(
            instance, {"n_samples": [], "solve_times": [], "runs": 0, "solved": 0, "n_interactions": []})
        entry["runs"] += 1
        if row.n_samples and row.n_samples > 0 and np.isfinite(row.n_samples):
            entry["n_samples"].append(row.n_samples)
            entry["solve_times"].append(row.solve_time)
            entry["solved"] += 1


yasa_1h_data = algbench.read_as_pandas(str(yasa_1h_data_dir), read_row_with_excluded_algorithms([]))
yasa_1h_data["algorithm"] = "YASA (1h)"
print("Found", len(yasa_1h_data), "YASA (1h) entries")

add_data_from_dataframe(yasa_1h_data)

yasa_fixed_m_data = algbench.read_as_pandas(str(yasa_fixed_m_data_dir), read_row_with_excluded_algorithms([]))
yasa_fixed_m_data["algorithm"] = yasa_fixed_m_data["m"].apply(lambda m: f"YASA{m}" if m > 1 else "YASA")
print("Found", len(yasa_fixed_m_data), "YASA fixed m entries")

add_data_from_dataframe(yasa_fixed_m_data)

baseline_data = algbench.read_as_pandas(str(baseline_data_dir), read_row_with_excluded_algorithms([
    "YASA", "YASA3", "YASA5", "YASA10",  # these are replaced by newer YASA runs above
    "IPOG_FT", "IPOG_CSP", "IPOF_FT", "IPOF_CSP"  # these are not used in the paper as they contain bad samples.
]))
print("Found", len(baseline_data), "baseline entries")

add_data_from_dataframe(baseline_data)

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
