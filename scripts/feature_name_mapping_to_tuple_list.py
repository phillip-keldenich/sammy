import itertools
import sys
import json


def canon_tuple(name1, val1, name2, val2):
    if name1 < name2:
        return (name1, val1, name2, val2)
    elif name1 > name2:
        return (name2, val2, name1, val1)
    else:
        raise RuntimeError("BUG: same name twice")


def generate_tuples(feature_name_mapping: list[dict[str, bool]]):
    tuples: set[tuple[str, bool, str, bool]] = set()
    for mapping in feature_name_mapping:
        for name1, name2 in itertools.combinations(mapping.keys(), 2):
            val1 = mapping[name1]
            val2 = mapping[name2]
            tuples.add(canon_tuple(name1, val1, name2, val2))
    return tuples


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: feature_name_mapping_to_tuple_list feature_name_mapping")
        sys.exit(1)

    with open(sys.argv[1], "r") as f:
        fnm = json.load(f)

    print(json.dumps(list(generate_tuples(fnm))))
