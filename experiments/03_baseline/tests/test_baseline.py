import tempfile
from pathlib import Path
import sys
import json

sys.path.append(str(Path(__file__).parent.parent))

from solver import get_instance_from_archive, BaselineAlgorithm


def test_baseline_solver():
    instance_name = "toybox"
    archive_path = Path(__file__).parent / "../../../instances/benchmark_models.zip"

    with tempfile.TemporaryDirectory() as tmpdir:
        instance_path = get_instance_from_archive(instance_name=instance_name, out_dir=tmpdir,
                                                  archive_path=archive_path)

        alg = BaselineAlgorithm(file_path=instance_path,
                                algorithm="YASA", seed=42)

        sample = alg.optimize(10)
        assert sample is not None, "Sample should not be None"

    with open(Path(__file__).parent / "toybox_sol.json", "r") as f:
        desired_solution = json.load(f)

    assert len(sample) == len(desired_solution), "Sample length does not match desired solution length"

    for i, s in enumerate(sample):
        for key, val in s.items():
            assert key in desired_solution[i], f"Key {key} not found in desired solution"
            assert val == desired_solution[i][key], f"Value for key {key} does not match desired solution"
