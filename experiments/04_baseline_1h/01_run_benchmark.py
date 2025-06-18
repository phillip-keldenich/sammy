import subprocess
import time
from pathlib import Path

import pandas as pd
from algbench import Benchmark

from solver.instance import get_instance_from_path
import tempfile
import slurminade

slurminade.update_default_configuration(
    partition="alg",
    constraint="alggen05",
    mail_type="FAIL",
    exclusive=True,
)
slurminade.set_dispatch_limit(300)

instances_path = "../../sammy_benchmark_instances"
output_dir = "yasa_1h"
instances = [
    # "ChatClient",
    # "DMIE",
    # "E-Shop",
    "Automotive01",
    "Automotive02_V1",
    "Automotive02_V2",
    "Automotive02_V3",
    "Automotive02_V4",
    "financial-services-2018-04-23",
    "linux_2_6_28_6",
    "linux_2_6_33_3"
]
time_limit = 3600  # seconds

benchmark = Benchmark(output_dir, hide_output=False)


def parse_sample(out_csv):
    # This matches all valid samples that were generated.

    if not out_csv.exists():
        return None

    try:
        t = pd.read_csv(out_csv, sep=";", index_col="Configuration")
    except Exception:
        print(f"Error reading {out_csv}. It may be empty or corrupted.")
        return None

    samples = []
    t.apply(
        lambda row: samples.append({k: v == "+" for k, v in row.items()}), axis=1
    )
    return samples


@slurminade.slurmify()
def solve_yasa(instance_name: str, alg_params: dict):
    def solve(instance_name, alg_params):
        assert alg_params["algorithm"] == "yasa", "Only YASA algorithm is supported in this benchmark."
        with tempfile.TemporaryDirectory() as tmpdir:
            out_csv = Path(tmpdir) / "sample.csv"

            instance_path = get_instance_from_path(instance_name=instance_name,
                                                   out_dir=tmpdir,
                                                   instances_path=instances_path)

            solve_time_start = time.time()
            runner = subprocess.run(
                [
                    "java",
                    "-jar",
                    "formula-analysis-sat4j-0.1.1-SNAPSHOT-all.jar",
                    "--command",
                    "de.featjar.formula.analysis.cli.TWiseCommand",
                    "--input",
                    str(instance_path),
                    "--output",
                    str(out_csv.absolute()),
                    "--t",
                    "2",
                    "--i",
                    "-1",
                    "--timeout",
                    str(alg_params["time_limit"]),
                    "--seed",
                    str(alg_params["seed"]),
                ],
                capture_output=False,
                check=True,
            )
            runner.check_returncode()
            solve_time = time.time() - solve_time_start
            samples = parse_sample(out_csv)
        return {
            "found_solution": samples is not None and len(samples) > 0,
            "samples": samples,
            "solve_time": solve_time,
        }

    benchmark.add(solve, instance_name=instance_name, alg_params=alg_params)


@slurminade.slurmify(
    mail_type="ALL"
)
def compress():
    """
    Compress the output directory
    """
    benchmark.compress()


if __name__ == "__main__":

    if not Path("formula-analysis-sat4j-0.1.1-SNAPSHOT-all.jar").exists():
        msg = "Please put formula-analysis-sat4j-0.1.1-SNAPSHOT-all.jar in this folder."
        raise RuntimeError(
            msg
        )

    with slurminade.JobBundling(max_size=5):  # automatically bundles up to 2 tasks
        for instance in instances:
            solve_yasa.distribute(instance_name=instance, alg_params={
                "algorithm": "yasa",
                "time_limit": time_limit,
                "seed": 42  # Fixed seed for reproducibility
            })

    slurminade.join()
    compress.distribute()
