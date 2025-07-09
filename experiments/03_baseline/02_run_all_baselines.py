import time

from algbench import Benchmark

from solver import BaselineAlgorithm, get_instance_from_path
import slurminade
import tempfile
import random as rd

slurminade.update_default_configuration(
    partition="alg",
    exclusive=True,
    constraint="alggen05",
    mail_type="FAIL"
)  # global options for slurm

instances_path = "../../sammy_benchmark_instances"
output_dir = "baseline_1h"
algorithms = ["YASA", "YASA3", "YASA5", "YASA10", "IC", "CH", "IL", "IPOG_FT", "IPOG_CSP", "IPOF_FT", "IPOF_CSP"]

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
    "Automotive02_V1",
    "Automotive02_V2",
    "Automotive02_V3",
    "Automotive02_V4",
    "financial-services-2018-04-23",
    "linux_2_6_28_6",
    "linux_2_6_33_3"
]
time_limit =  3600 # seconds
n_repeats = 5  # number of repeats for each instance and algorithm

benchmark = Benchmark(output_dir, hide_output=False)
benchmark.capture_logger("baseline")


@slurminade.slurmify()
def solve_baseline(instance_name: str, alg_params: dict, seed: int):
    def solve(instance_name, alg_params, seed):
        with tempfile.TemporaryDirectory() as tmpdir:
            #instance_path = get_instance_from_archive(instance_name=instance_name, out_dir=tmpdir,
            #                                          archive_path=instances_path)

            instance_path = get_instance_from_path(instance_name=instance_name,
                                                   out_dir=tmpdir,
                                                   instances_path=instances_path)


            alg = BaselineAlgorithm(file_path=instance_path,
                                    seed=seed,
                                    algorithm=alg_params["algorithm"])

            solve_time_start = time.time()
            samples = alg.optimize(alg_params["time_limit"])
            solve_time = time.time() - solve_time_start

        return {
            "configuration": alg._configuration,
            "samples": samples,
            "solve_time": solve_time,
        }

    benchmark.add(solve, instance_name=instance_name, alg_params=alg_params, seed=seed)


@slurminade.slurmify(
    mail_type="ALL"
)
def compress():
    """
    Compress the output directory
    """
    benchmark.compress()

if __name__ == "__main__":

    all_alg_params = [
        {
            "algorithm": algo,
            "time_limit": time_limit,
        }

        for algo in algorithms
    ]

    rd.shuffle(instances)

    with slurminade.JobBundling(max_size=10):  # automatically bundles up to 2 tasks
        for i in range(n_repeats):
            for instance in instances:
                for alg_params in all_alg_params:
                    solve_baseline.distribute(instance_name=instance, alg_params=alg_params, seed=i)

    slurminade.join()
    compress.distribute()
