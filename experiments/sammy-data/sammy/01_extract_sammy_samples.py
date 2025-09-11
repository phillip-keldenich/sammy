import subprocess
import json
from collections import defaultdict
from pathlib import Path

from tqdm import tqdm
import lzma

import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "sammy",
        type=str,
    )

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
        "lcm",
        # 8 new instances of larger size
        "Automotive01",
        # "Automotive02_V1",
        # "Automotive02_V2",
        # "Automotive02_V3",
        # "Automotive02_V4",
        "financial-services-2018-04-23",
        "linux_2_6_28_6",
        "linux_2_6_33_3"
    ]

    args = parser.parse_args()

    sammy_samples = defaultdict(list)
    sammy_dir = Path(args.sammy)

    for instance in tqdm(instances):

        script = sammy_dir / "scripts/solution_to_feature_name_mapping.py"
        instance_file = sammy_dir / f"sammy_benchmark_instances/{instance}.scm.json.xz"

        for i in range(1, 6):
            solution_file = sammy_dir / f"experiments/02_run_default_params/02_output/{instance}_repeat{i}.out.json.xz"

            result = subprocess.run(
                ["python", script, instance_file, solution_file],  # pass args as a list → no shell quoting issues
                capture_output=True,  # grab stdout / stderr
                text=True,  # decode bytes → str
                check=True  # raise if the script returned non-zero
            )

            json_str = result.stdout  # plain text JSON
            sample = json.loads(json_str)  # now it’s a Python dict/list, etc.
            sammy_samples[instance].append(sample)

    print("Saving to file...")
    with lzma.open("02_sammy_samples.json.xz", "wt") as f:
        json.dump(sammy_samples, f, indent=4)
