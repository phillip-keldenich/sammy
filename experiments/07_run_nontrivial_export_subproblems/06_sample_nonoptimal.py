import json
import lzma
import random

random.seed(42)

with lzma.open("05_amended_outcomes.json.xz", "rt") as f:
    data = json.load(f)
nonoptimal = data['out_of_iterations'] + data['mes_optimal']
random.shuffle(nonoptimal)
with open("07_shuffled_nonoptimal.json", "w") as f:
    json.dump(nonoptimal, f, indent=2)
