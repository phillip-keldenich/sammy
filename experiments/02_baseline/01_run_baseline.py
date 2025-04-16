import time
from algbench import Benchmark

from solver import BaselineAlgorithm, get_instance_from_archive
import slurminade
import tempfile

slurminade.update_default_configuration(
    partition="alg",
    exclusive=True,
    constraint="alggen05",
    mail_type="FAIL"
)  # global options for slurm

instances_path = "../../instances/benchmark_models.zip"
output_dir = "01_baseline_3h"
algorithms = ["YASA", "YASA3", "YASA5", "YASA10"]
instances = ["Automotive02_V1", "Automotive02_V2", "Automotive02_V3", "Automotive02_V4"]
time_limit = 10800 # seconds

benchmark = Benchmark(output_dir)
benchmark.capture_logger("baseline")


@slurminade.slurmify()
def solve_baseline(instance_name: str, alg_params: dict):
    def solve(instance_name, alg_params):
        with tempfile.TemporaryDirectory() as tmpdir:
            instance_path = get_instance_from_archive(instance_name=instance_name, out_dir=tmpdir,
                                                      archive_path=instances_path)

            alg = BaselineAlgorithm(file_path=instance_path,
                                    algorithm=alg_params["algorithm"])

            solve_time_start = time.time()
            samples = alg.optimize(alg_params["time_limit"])
            solve_time = time.time() - solve_time_start

        return {
            "configuration": alg._configuration,
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

    all_alg_params = [
        {
            "algorithm": algo,
            "time_limit": time_limit,
        }

        for algo in algorithms
    ]

    with slurminade.JobBundling(max_size=2):  # automatically bundles up to 2 tasks
        for instance in instances:
            for alg_params in all_alg_params:
                solve_baseline.distribute(instance_name=instance, alg_params=alg_params)

    slurminade.join()
    compress.distribute()
