import lzma
import json

with lzma.open("04_sammy_n_interactions.json.xz", "rt") as f:
    n_interactions_for_instance = json.load(f)

with lzma.open("../baselines/03a_baseline_summary_data_with_interactions.json.xz", "rt") as f:
    baseline_data = json.load(f)

baseline_n_interactions = {
    instance: int(baseline_data["YASA"][instance]["n_interactions"][0]) if baseline_data["YASA"][instance]["n_interactions"] else None
    for instance in baseline_data["YASA"]
}

count_invalid = 0
count_checked = 0

for instance in n_interactions_for_instance:
    for count in n_interactions_for_instance[instance]:
        if not (baseline_n_interactions[instance] == count
                or baseline_n_interactions[instance] is None):
            print(f"Mismatch for instance {instance}: baseline {baseline_n_interactions[instance]}, sammy {count}")
            count_invalid += 1
        else:
            count_checked += 1

if count_invalid > 0:
    print("Found", count_invalid, "instances with mismatched interaction counts.")
else:
    print("All", count_checked, "interaction counts match between baseline and sammy data.")
