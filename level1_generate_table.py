import lzma
from pathlib import Path
import json
import pandas as pd
import math
import sys
import algbench

SAMPLNS_PATH = Path(__file__).parent / "experiments" / "sammy-data" / "samplns" / "samplns_1h"


def read_row_samplns(row):
    return {
        "lower_bound": row["result"]["lower_bound"],
        "upper_bound": row["result"]["upper_bound"],
        "time_used_by_yasa": row["result"]["time_used_by_yasa"],
        "runtime": row["runtime"],
        "instance": row["parameters"]["args"]["instance_name"],
    }


def generate_expected_runs(instances_dir: Path):
    expected_runs = set()
    for instance_file in instances_dir.glob("*.json.xz"):
        instance_with_scm = instance_file.name[:-8]  # Remove .json.xz
        expected_run_names = [f"{instance_with_scm}.out{i}.json.xz"
                              for i in range(1, 6)]
        expected_runs.update(expected_run_names)
    return expected_runs


def instance_name_from_expected_run(result_name: str):
    return result_name.split(".scm")[0]


def first_event(events, event_type):
    for event in events:
        if event["type"] == event_type:
            return event
    return None


def records_to_table(records, output_file: Path):
    print("Generating table to", output_file)
    records = pd.DataFrame(records)

    samplns_path = Path(SAMPLNS_PATH)
    samplns_data = None

    if samplns_path.exists():
        print("Found SampLNS results and reading from", samplns_path)
        samplns_data = algbench.read_as_pandas(str(samplns_path), read_row_samplns)

    instance_data = []
    instances = set(records["instance"])
    for instance in instances:
        instance_records = records[records["instance"] == instance]
        assert len(instance_records) == 5
        entry = {
            "instance": instance,
            "successful_runs": sum(instance_records["success"]),
            "total_runs": 5,
            "median_time": instance_records["time"].median(),
            "min_lb": instance_records["lb"].min(),
            "med_lb": instance_records["lb"].median(),
            "max_lb": instance_records["lb"].max(),
            "min_ub": instance_records["ub"].min(),
            "med_ub": instance_records["ub"].median(),
            "max_ub": instance_records["ub"].max(),
            "features": instance_records["features"].min(),
            "concrete": instance_records["concrete"].min(),
            "clauses": instance_records["clauses"].min(),
            "interactions": instance_records["interactions"].min(),
        }

        if samplns_data is not None:
            samplns_records = samplns_data[samplns_data["instance"] == instance]
            entry.update({
                "samplns_runs": len(samplns_records),
                "samplns_min_lb": samplns_records["lower_bound"].min() if len(samplns_records) > 0 else None,
                "samplns_med_lb": samplns_records["lower_bound"].median() if len(samplns_records) > 0 else None,
                "samplns_max_lb": samplns_records["lower_bound"].max() if len(samplns_records) > 0 else None,
                "samplns_min_ub": samplns_records["upper_bound"].min() if len(samplns_records) > 0 else None,
                "samplns_med_ub": samplns_records["upper_bound"].median() if len(samplns_records) > 0 else None,
                "samplns_max_ub": samplns_records["upper_bound"].max() if len(samplns_records) > 0 else None,
                "samplns_median_time": samplns_records["runtime"].median() if len(samplns_records) > 0 else None,
            })

        instance_data.append(entry)
    instance_data = pd.DataFrame(instance_data)
    instance_data = instance_data.sort_values(by="interactions")
    if output_file.suffix.lower() == ".csv":
        instance_data.to_csv(output_file, index=False,
                             float_format="%.3f", header=True)
    elif output_file.suffix.lower() in [".xlsx", ".xls"]:
        instance_data.to_excel(output_file, sheet_name="Results (Table E.1)",
                               index=False, header=True)


def generate_table(instances_dir: Path, results_dir: Path, output_file: Path):
    expected_runs = generate_expected_runs(instances_dir)
    records = []

    for expected_run in expected_runs:
        instance_name = instance_name_from_expected_run(expected_run)
        run_index = int(expected_run.split(".out")[-1][0])
        expected_file = results_dir / expected_run
        record = {
            "instance": instance_name,
            "run_index": run_index,
            "success": False,
            "time": math.inf,
            "lb": 0,
            "ub": math.inf,
            "features": None,
            "concrete": None,
            "clauses": None,
            "interactions": None
        }
        if expected_file.exists():
            with lzma.open(expected_file, "rt") as f:
                data = json.load(f)
                record["success"] = True
                record["ub"] = len(data["best_solution"])
                record["lb"] = data["lb"]
                record["time"] = data["events"][-1]["time"]
                input_event = first_event(data["events"], "INPUT_READ")
                if record["features"] is None and input_event:
                    record["features"] = input_event["num_vars"]
                    record["concrete"] = input_event["num_concrete"]
                    record["clauses"] = input_event["num_clauses"]
                elif input_event:
                    assert (record["features"] == input_event["num_vars"] and
                            record["concrete"] == input_event["num_concrete"] and
                            record["clauses"] == input_event["num_clauses"])
                initial_phase_event = first_event(data["events"],
                                                  "INITIAL_PHASE_DONE")
                if record["interactions"] is None and initial_phase_event:
                    record["interactions"] = initial_phase_event["num_interactions"]
                elif initial_phase_event:
                    assert record["interactions"] == initial_phase_event["num_interactions"]
        records.append(record)
    records_to_table(records, output_file)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python level1_generate_table.py <results_dir> <output_file>")
        sys.exit(1)
    this_dir = Path(__file__).absolute().parent
    instances_dir = this_dir / "sammy_benchmark_instances"
    results_dir = Path(sys.argv[1])
    output_file = Path(sys.argv[2])
    generate_table(instances_dir, results_dir, output_file)
