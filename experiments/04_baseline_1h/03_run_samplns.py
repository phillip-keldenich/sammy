import logging
from pathlib import Path

import slurminade
from algbench import Benchmark
from samplns.baseline.algorithm import BaselineAlgorithm
from samplns.instances import parse, parse_dimacs
from samplns.lns.lns import LnsObserver
from samplns.lns.neighborhood import Neighborhood, RandomNeighborhood
from samplns.simple import SampLns
from samplns.utils import Timer
from solver.instance import get_instance_from_path
import tempfile

instances_path = "../../sammy_benchmark_instances"
output_dir = "samplns_1h"
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
samplns_config = {
    "iterations": 10_000,
    "iteration_time_limit": 60.0,
    "time_limit": time_limit
}

benchmark = Benchmark(output_dir, hide_output=False)

# ================================================
# SLURM CONFIGURATION FOR DISTRIBUTED EXECUTION
# ------------------------------------------------
slurminade.update_default_configuration(
    partition="alg",
    constraint="alggen05",
    exclusive=True,
    mail_type="FAIL",
)
slurminade.set_dispatch_limit(300)


class MyLnsLogger(LnsObserver):
    """
    A logger that will save information on the SampLNS iterations for us to evaluate
    later.
    """

    def __init__(self):
        self.timer = Timer(0)
        self.iterations = []

    def report_neighborhood_optimization(self, neighborhood: Neighborhood):
        self.iterations[-1]["nbrhd_tuples"] = len(neighborhood.missing_tuples)
        self.iterations[-1]["nbrhd_confs"] = len(neighborhood.initial_solution)

    def report_iteration_begin(self, iteration: int):
        self.iterations.append({})

    def report_iteration_end(
            self, iteration: int, runtime: float, lb: int, solution, events
    ):
        self.iterations[-1].update(
            {
                "iteration": iteration,
                "lb": lb,
                "ub": len(solution),
                "time": self.timer.time(),
                "iteration_time": runtime,
                "events": events,
            }
        )


logging.getLogger("SampLNS").addHandler(logging.StreamHandler())
logging.getLogger("SampLNS.CPSAT").setLevel(logging.WARNING)
benchmark.capture_logger("SampLNS", logging.INFO)


@slurminade.slurmify
def run_distributed(instance_name: str, rep: int):
    benchmark.add(
        run_samplns,
        instance_name,
        iterations=samplns_config["iterations"],
        iteration_time_limit=samplns_config["iteration_time_limit"],
        time_limit=samplns_config["time_limit"],
        verify=True,
        fast_verify=True,
        rep=rep,
    )


def run_samplns(
        instance_name: str,
        iterations,
        iteration_time_limit,
        time_limit,
        verify,
        fast_verify,
        rep: int,
):
    """
    Running SampLNS on an initial sample.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        instance_path = get_instance_from_path(instance_name=instance_name,
                                               out_dir=tmpdir,
                                               instances_path=instances_path)
        if instance_path.suffix == ".xml":
            instance = parse(str(instance_path))
        else:
            instance = parse_dimacs(str(instance_path))

        # measure time for running yasa
        timer = Timer()
        sample = BaselineAlgorithm(str(instance_path)).optimize(time_limit)
        _yasa_time_used = timer.time()
        assert sample

        remaining_time = time_limit - _yasa_time_used
        solver = None
        logger = MyLnsLogger()
        if remaining_time > 0:
            # setup (needs time measurement as already involves calculations)
            solver = SampLns(
                instance=instance,
                initial_solution=sample,
                neighborhood_selector=RandomNeighborhood(),
                observer=logger,
            )

            solver.optimize(
                iterations=iterations,
                iteration_timelimit=iteration_time_limit,
                timelimit=remaining_time,
            )

            solution = solver.get_best_solution(verify=verify, fast_verify=fast_verify)
        else:
            solution = sample

        # get optimized sample and verify its correctness (takes some time).
        return {
            "samplns_used": solver is not None,
            "time_used_by_yasa": _yasa_time_used,
            "timelimit_for_samplns": remaining_time,
            "solution": solution,
            "lower_bound": solver.get_lower_bound() if solver else 1,
            "upper_bound": len(solution),
            "optimal": (
                solver.get_lower_bound() == len(solver.get_best_solution())
                if solver
                else False
            ),
            "iteration_info": logger.iterations,
        }


@slurminade.slurmify(mail_type="ALL")
def pack_after_finish():
    """
    Compress the database after the experiment has finished.
    """
    benchmark.compress()


if __name__ == "__main__":

    for rep in range(5):
        for instance in instances:
            run_distributed.distribute(instance, rep)
    slurminade.join()
    pack_after_finish.distribute()
