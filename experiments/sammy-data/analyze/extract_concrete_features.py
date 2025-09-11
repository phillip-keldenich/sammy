from parser.parse_xml import parse
from parser.parse_dimacs import parse_dimacs
from pathlib import Path
import json

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "directory",
        type=str,
        help="Path to the instances.",
    )

    args = parser.parse_args()

    features_for_instance = {}



    for file_path in list(Path(args.directory).glob("*.xml")) + list(Path(args.directory).glob("*.dimacs")):
        instance = parse(str(file_path)) if file_path.suffix == ".xml" else parse_dimacs(str(file_path))
        features_for_instance[file_path.stem] = instance

    print(len(features_for_instance))
    print(features_for_instance.keys())

    with open("concrete_features_for_instance.json", "w") as f:
        json.dump(
            {k: [str(f) for f in v.features] for k, v in features_for_instance.items()},
            f,
            indent=4,
        )