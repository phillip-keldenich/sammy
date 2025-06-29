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
    # 46 instances from SampLNS paper
    "busybox-1_29_2",
    "busybox_2007-01-24_09-14-09",
    "busybox_2020-12-16_21-53-05",
    "fiasco_2017-09-26_11-30-56",
    "fiasco_2020-12-01_14-09-14",
    "soletta_2015-06-26_18-38-56",
    "soletta_2017-03-09_21-02-40",
    "toybox_2006-10-31_23-30-06",
    "toybox_2020-12-06_00-02-46",
    "uclibc_2008-06-05_13-46-47",
    "uclibc_2020-12-24_11-54-53",
    "APL-Model",
    "APL",
    "BankingSoftware",
    "BattleofTanks",
    "ChatClient",
    "DMIE",
    "E-Shop",
    "EMBToolkit",
    "FameDB",
    "FeatureIDE",
    "FreeBSD-8_0_0",
    "PPU",
    "SafeBali",
    "SortingLine",
    "TightVNC",
    "Violet",
    "WaterlooGenerated",
    "XSEngine",
    "aaed2000",
    "am31_sim",
    "atlas_mips32_4kc",
    "axTLS",
    "berkeleyDB1",
    "berkeleyDB2",
    "busybox-1_18_0",
    "calculate",
    "car",
    "dell",
    "eCos-3-0_i386pc",
    "ea2468",
    "email",
    "financial_services",
    "fs_2017-05-22",
    "gpl",
    "integrator_arm7",
    "lcm"
    # 8 new instances of larger size
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
n_repeats = 5

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

    for i in range(n_repeats):
        for instance in instances:
            solve_yasa.distribute(instance_name=instance, alg_params={
                "algorithm": "yasa",
                "time_limit": time_limit,
                "seed": i
            })

    slurminade.join()
    compress.distribute()
