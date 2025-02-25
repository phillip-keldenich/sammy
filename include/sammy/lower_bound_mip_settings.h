#ifndef SAMMY_LOWER_BOUND_SOLVER_STATE_H_INCLUDED_
#define SAMMY_LOWER_BOUND_SOLVER_STATE_H_INCLUDED_

#include <cstddef>
#include <limits>

namespace sammy {

enum class SolverState {
    OPTIMUM_ON_SUBGRAPH,
    NO_IMPROVEMENT_FOUND,
    IMPROVEMENT_FOUND,
    TIMEOUT_NO_IMPROVEMENT,
    TIMEOUT_IMPROVEMENT
};

struct LowerBoundMIPConfig {
    bool solve_exactly_full_mip = false;
    bool solve_exactly_cut_and_price = false;
    bool quiet_gurobi = false;
    std::size_t stop_cuts_iteration_window = 5;
    std::size_t separation_rounds_before_pricing = 5;
    double stop_cuts_below_relative_improvement = 0.03;
    double total_mip_timeout = std::numeric_limits<double>::infinity();
};

struct SATColoringConfig {
    bool initial_sat_coloring = false;
    double initial_sat_coloring_timeout = 10.0;
};

}

#endif
