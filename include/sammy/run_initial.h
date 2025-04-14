#ifndef SAMMY_RUN_INITIAL_H_INCLUDED_
#define SAMMY_RUN_INITIAL_H_INCLUDED_

#include "fast_clique.h"
#include "initial_coloring_heuristic.h"
#include "input.h"
#include "learn_infeasibilities.h"
#include "output.h"
#include "simplification.h"
#include <limits>

namespace sammy {

/**
 * @brief Configuration structure for the initial run of the algorithm.
 * 
 * This structure contains parameters to configure the behavior of the 
 * initial run, including time limits, iteration goals, and other options.
 */
struct RunInitialConfig {
    /**
     * @brief Minimum time (in seconds) for the initial run.
     * Default value is 60.0 seconds.
     */
    double min_time = 60.0,

    /**
     * @brief Maximum time (in seconds) for the initial run.
     * Default value is infinity.
     */
    double max_time = std::numeric_limits<double>::infinity();

    /**
     * @brief Number of random clique restarts per iteration.
     * Default value is 40.
     */
    std::size_t random_clique_restarts_per_iteration = 40;

    /**
     * @brief Target number of iterations for the initial run.
     * Default value is 20.
     */
    std::size_t goal_iterations = 20;

    /**
     * @brief Flag to enable or disable simplification during the initial run.
     * Default value is true (enabled).
     */
    bool simplify = true;

    /**
     * @brief Flag to enable or disable quiet mode (suppress output).
     * Default value is true (quiet mode enabled).
     */
    bool quiet = true;
};

class InitialHeuristicData {
  public:
    InitialHeuristicData(const ClauseDB& simplified, Var n_concrete,
                         EventRecorder* recorder,
                         const RunInitialConfig& config,
                         ThreadGroup<void>* thread_pool)
        : simplified(simplified), inf_map(n_concrete),
          col_solver(&this->simplified, n_concrete, thread_pool, &inf_map),
          clq_solver(SharedDBPropagator(&this->simplified), thread_pool),
          recorder(recorder), config(config) {
        col_solver.set_quiet(config.quiet);
    }

    void first_solution() {
        recorder->store_event("FIRST_SOLVE_BEGINS");
        col_solver.initialize_feasibilities();
        col_solver.color_lazy();
        recorder->store_event("FIRST_COLORING_DONE");
        col_solver.extract_feasibilities();
        recorder->store_event("FEASIBILITIES_EXTRACTED");
        learn_infeasibilities(simplified, &inf_map);
        recorder->store_event("INFEASIBILITIES_LEARNED");
        update_best_ub(0, col_solver.class_spawners());
        recorder->store_event(
            "FIRST_SOLVE_DONE",
            OutputObject{{"num_configs", col_solver.all_classes().size()}});
        all_spawners = col_solver.class_spawners();
        recorder->store_event("FIRST_LB_BEGINS");
        auto r = config.random_clique_restarts_per_iteration;
        update_best_lb(
            0, clq_solver.random_multistart_best_clique(r, all_spawners));
        recorder->store_event("FIRST_LB_DONE");
    }

    void follow_up_iteration(std::size_t iteration) {
        col_solver.reset_coloring();
        col_solver.color_lazy(best_lb);
        update_best_ub(iteration, col_solver.class_spawners());
        auto new_spawners = col_solver.class_spawners();
        add_spawners(new_spawners);
        auto r = config.random_clique_restarts_per_iteration;
        update_best_lb(iteration, clq_solver.random_multistart_best_clique(
                                      r / 2, all_spawners));
        update_best_lb(iteration, clq_solver.random_multistart_best_clique(
                                      r / 2, new_spawners));
    }

    void rerun_primal(const std::vector<Vertex>& initial) {
        auto iteration = std::numeric_limits<std::size_t>::max();
        col_solver.reset_coloring();
        col_solver.color_lazy(initial);
        update_best_ub(iteration, col_solver.class_spawners());
        add_spawners(col_solver.class_spawners());
    }

    void rerun_primal_preset(const std::vector<Vertex>& initial) {
        auto iteration = std::numeric_limits<std::size_t>::max();
        col_solver.color_lazy(initial);
        update_best_ub(iteration, col_solver.class_spawners());
        add_spawners(col_solver.class_spawners());
    }

    void add_spawners(const std::vector<Vertex>& new_spawners) {
        all_spawners.insert(all_spawners.end(), new_spawners.begin(),
                            new_spawners.end());
        std::sort(all_spawners.begin(), all_spawners.end());
        all_spawners.erase(
            std::unique(all_spawners.begin(), all_spawners.end()),
            all_spawners.end());
    }

    void update_best_ub(std::size_t iteration,
                        const std::vector<Vertex>& spawners) {
        update_best_ub(iteration, col_solver.all_classes(), spawners);
    }

    void update_best_ub(std::size_t iteration,
                        const std::vector<SharedDBPropagator>& propagators,
                        const std::vector<Vertex>& spawners) {
        if (best_ub.empty() || propagators.size() < best_ub.size()) {
            best_ub.clear();
            std::transform(propagators.begin(), propagators.end(),
                           std::back_inserter(best_ub),
                           [](const SharedDBPropagator& prop) {
                               return prop.extract_assignment();
                           });
            best_spawners = spawners;
            recorder->store_event(
                "UB_IMPROVED",
                OutputObject{{"num_configs", propagators.size()},
                             {"iteration", iteration + 1}},
                "num_configs", "iteration");
        }
    }

    void update_best_lb(std::size_t iteration, std::vector<Vertex> clique) {
        if (clique.size() > best_lb.size()) {
            recorder->store_event(
                "LB_IMPROVED",
                OutputObject{{"num_interactions", clique.size()},
                             {"iteration", iteration + 1}},
                "num_interactions", "iteration");
            best_lb = std::move(clique);
        }
    }

    std::vector<std::vector<bool>>
    reconstruct_solution(SimplifyDatastructure& simplify_ds,
                         SimplifiedInstance& simplified_inst) {
        return simplify_ds.reconstruct_sample(simplified_inst, best_ub);
    }

    std::vector<Vertex> reconstruct_lb(SimplifyDatastructure& simplify_ds,
                                       SimplifiedInstance& simplified_inst) {
        return simplify_ds.reconstruct_lb(simplified_inst, best_lb);
    }

    const std::vector<std::vector<bool>>& get_solution() const noexcept {
        return best_ub;
    }

    const std::vector<Vertex>& get_bound() const noexcept { return best_lb; }

    std::size_t best_lb_value() const noexcept {
        return (std::max)(std::size_t(1), best_lb.size());
    }

    std::size_t best_ub_value() const noexcept { return best_ub.size(); }

    ColoringHeuristicSolver& primal_solver() noexcept { return col_solver; }

    ParallelFastCliqueBuilder& dual_solver() noexcept { return clq_solver; }

    const ColoringHeuristicSolver& primal_solver() const noexcept {
        return col_solver;
    }

    const ParallelFastCliqueBuilder& dual_solver() const noexcept {
        return clq_solver;
    }

    PairInfeasibilityMap& infeasibility_map() noexcept { return inf_map; }

    const PairInfeasibilityMap& infeasibility_map() const noexcept {
        return inf_map;
    }

    const std::vector<Vertex>& get_all_spawners() const noexcept {
        return all_spawners;
    }

    ClauseDB& formula() noexcept { return simplified; }
    const ClauseDB& formula() const noexcept { return simplified; }
    Var num_concrete() const noexcept { return inf_map.get_n_concrete(); }

    const std::vector<Vertex>& get_best_spawners() const noexcept {
        return best_spawners;
    }

  private:
    ClauseDB simplified;
    PairInfeasibilityMap inf_map;
    ColoringHeuristicSolver col_solver;
    ParallelFastCliqueBuilder clq_solver;
    EventRecorder* recorder;
    std::vector<std::vector<bool>> best_ub;
    std::vector<Vertex> best_spawners;
    std::vector<Vertex> best_lb;
    std::vector<Vertex> all_spawners;
    RunInitialConfig config;
};

inline std::pair<std::size_t, std::size_t> run_initial_heuristic(
    const InputData& input, OutputObject& output, EventRecorder& recorder,
    RunInitialConfig config, ThreadGroup<void>& tgroup,
    std::optional<InitialHeuristicData>& initial_data,
    std::optional<SimplifyDatastructure>& simplifier,
    std::optional<SimplifiedInstance>& simplified,
    std::function<std::vector<std::vector<bool>>()>& get_best_solution_fn,
    std::function<std::vector<Vertex>()>& get_best_mut_ex_set) {
    const ClauseDB* formula = &input.formula;
    Var num_concrete = input.num_concrete;
    if (config.simplify) {
        SimplificationStats simplification_stats;
        simplification_stats.capture_before(input.formula, input.num_concrete);
        simplifier.emplace(input.formula, input.num_concrete);
        simplified.emplace(
            sammy::run_simplifier(*simplifier, simplification_stats));
        num_concrete = simplified->num_concrete;
        formula = &simplified->formula;
        simplification_stats.capture_after(*formula, num_concrete);
        recorder.store_event("SIMPLIFICATION_DONE");
        export_simplification(output, *simplifier, *simplified,
                              &simplification_stats);
    }
    auto before_initial = Clock::now();
    initial_data.emplace(*formula, num_concrete, &recorder, config, &tgroup);
    auto should_stop = [&](std::size_t iteration) -> bool {
        double current_time = seconds_between(before_initial, Clock::now());
        if (iteration >= config.goal_iterations &&
            current_time >= config.min_time)
            return true;
        if (initial_data->best_lb_value() >= initial_data->best_ub_value())
            return true;
        return current_time >= config.max_time;
    };
    if (config.simplify) {
        get_best_solution_fn = [&initial_data, &simplifier, &simplified]() {
            return initial_data->reconstruct_solution(*simplifier, *simplified);
        };
        get_best_mut_ex_set = [&initial_data, &simplifier, &simplified]() {
            return initial_data->reconstruct_lb(*simplifier, *simplified);
        };
    } else {
        get_best_solution_fn = [&initial_data]() {
            return initial_data->get_solution();
        };
        get_best_mut_ex_set = [&initial_data]() {
            return initial_data->get_bound();
        };
    }
    initial_data->first_solution();
    std::size_t actual_iterations;
    for (actual_iterations = 1; !should_stop(actual_iterations);
         ++actual_iterations)
    {
        initial_data->follow_up_iteration(actual_iterations);
    }
    export_solution(output, get_best_solution_fn(), "initial_heuristic",
                    OutputObject{{"iterations", actual_iterations}});
    export_bound(output, get_best_mut_ex_set(), "initial_heuristic",
                 OutputObject{{"random_clique_restarts",
                               config.random_clique_restarts_per_iteration}});
    return {initial_data->best_lb_value(), initial_data->best_ub_value()};
}

inline std::pair<std::size_t, std::size_t>
run_initial_heuristic(const InputData& input, OutputObject& output,
                      EventRecorder& recorder, RunInitialConfig config,
                      ThreadGroup<void>* thread_pool = nullptr) {
    std::optional<ThreadGroup<void>> tg;
    if (!thread_pool) {
        tg.emplace();
        thread_pool = &*tg;
    }
    std::optional<InitialHeuristicData> initial_data;
    std::optional<SimplifyDatastructure> simplifier;
    std::optional<SimplifiedInstance> simplified;
    std::function<std::vector<std::vector<bool>>()> get_best_solution;
    std::function<std::vector<Vertex>()> get_best_mut_ex_set;
    return run_initial_heuristic(input, output, recorder, config, *thread_pool,
                                 initial_data, simplifier, simplified,
                                 get_best_solution, get_best_mut_ex_set);
}

} // namespace sammy

#endif
