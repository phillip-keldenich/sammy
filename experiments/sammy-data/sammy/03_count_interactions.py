from collections import defaultdict

from sample_analyzer.verify import count_interactions
from pathlib import Path
import json
import lzma

CONCRETE_FEATURES_FILE = Path(__file__).parent.parent / "instance_data" / "concrete_features_for_instance.json"

with open(CONCRETE_FEATURES_FILE, "r") as f:
    concrete_features_for_instance = json.load(f)

with lzma.open("02_sammy_samples.json.xz", "rt") as f:
    samples_for_instance = json.load(f)

n_interactions_for_instance = defaultdict(list)

for instance in samples_for_instance:

    print(instance)

    for sample in samples_for_instance[instance]:
        try:
            n_interactions_for_instance[instance].append(count_interactions(
                sample,
                concrete_features_for_instance[instance]
            ))
        except ValueError:
            print(f"Error processing instance {instance}. Skipping.")
            n_interactions_for_instance[instance] = None

with lzma.open("04_sammy_n_interactions.json.xz", "wt") as f:
    json.dump(n_interactions_for_instance, f, indent=4)