#ifndef SAMMY_SUBPROBLEM_SOLVER_WITH_MES_H_INCLUDED_
#define SAMMY_SUBPROBLEM_SOLVER_WITH_MES_H_INCLUDED_

#include "global_stat_info.h"
#include "gurobi_clique_solver_g2.h"
#include "incremental_sat_lns.h"
#include <boost/circular_buffer.hpp>

namespace sammy {

constexpr std::size_t MAX_CUT_ROUNDS_WITHOUT_PRICING = 5;
constexpr std::size_t SUB_MES_SOLVER_FULL_GRAPH_THRESHOLD = 200;
constexpr std::size_t PRICING_SAMPLE_SIZE_THRESHOLD = 10'000;
constexpr std::size_t DEFAULT_ITERATIONS_LIMIT = 10;

/**
 * Solver for finding good MESs in a subproblem
 * before running another solver.
 */
class SubproblemMESSolver {
  public:
    SubproblemMESSolver(PortfolioSolver* portfolio, LNSSubproblem&& subproblem,
                        SharedDBPropagator prop, EventRecorder* recorder,
                        std::size_t worker_id,
                        std::size_t iterations_without_improvement = 30,
                        std::size_t iterations_limit = DEFAULT_ITERATIONS_LIMIT)
        : m_portfolio(portfolio), m_recorder(recorder), m_worker_id(worker_id),
          m_initial_best_global_mes_size(m_portfolio->get_best_mes_size()),
          m_single_thread(0), m_clauses(&portfolio->get_clauses()),
          m_subproblem(std::move(subproblem)),
          // make & extend subgraph; on interrupt, move subproblem back
          m_subgraph(p_make_universe_subgraph(subproblem)),
          m_base_solver(&m_subgraph, m_clauses->num_vars(), m_recorder,
                        m_subproblem.removed_configurations,
                        m_subproblem.mutually_exclusive_set, p_get_config()),
          m_iterations_without_improvement_limit(
              iterations_without_improvement),
          m_iterations_limit(iterations_limit) {
        m_base_solver.make_abortable();
    }

    void abort() { m_base_solver.abort(); }

    const std::vector<Vertex>& mes_vertices() const {
        return m_base_solver.get_best_mes();
    }

    void return_subproblem(LNSSubproblem&& subproblem) {
        m_subproblem = std::move(subproblem);
    }

    LNSSubproblem move_out_subproblem() { return std::move(m_subproblem); }

    const std::vector<DynamicBitset>& get_solution() const {
        throw std::logic_error(
            "SubproblemMESSolver asked for solution; this should not happen");
    }

    std::size_t num_lp_variables() const {
        return m_base_solver.num_lp_variables();
    }

    std::size_t num_lp_constraints() const {
        return m_base_solver.num_lp_constraints();
    }

    std::size_t num_lp_nonzeros() const {
        return m_base_solver.num_lp_nonzeros();
    }

    std::size_t num_lp_solve_calls() const {
        return m_base_solver.num_lp_solve_calls();
    }

    std::size_t num_graph_vertices() const noexcept {
        return m_base_solver.num_graph_vertices();
    }

    template<typename IterationCallback/*(SubproblemMESSolver&, bool was_interrupted, bool proved_mes_optimality, bool failed)*/>
    std::optional<bool> extended_solve(IterationCallback&& callback) {
        try {
            for (;;) {
                if (get_and_clear_interrupt_flag()) {
                    return std::nullopt;
                }
                auto iteration_result = p_solve_iteration();
                if (!iteration_result) {
                    callback(*this, true, false, false);
                    return std::nullopt;
                }
                if (iteration_result->improved) {
                    if (!p_handle_improvement_should_continue()) {
                        callback(*this, false, true, false);
                        return false;
                    }
                }
                if (iteration_result->optimal_on_subgraph) {
                    if (!p_optimal_on_subgraph_can_continue()) {
                        callback(*this, false, true, false);
                        return true;
                    }
                } else {
                    if (!p_non_optimal_on_subgraph_can_continue()) {
                        callback(*this, false, false, true);
                        return true;
                    }
                }
                if (!callback(*this, false, false, false)) {
                    return true;
                }
            }
        } catch (InterruptError&) {
            return std::nullopt;
        }
    }

    std::optional<bool> solve() {
        m_iterations = 0;
        try {
            for (;;) {
                if (get_and_clear_interrupt_flag()) {
                    return std::nullopt;
                }
                auto iteration_result = p_solve_iteration();
                if (!iteration_result) {
                    // interrupted
                    return std::nullopt;
                }
                if (iteration_result->improved) {
                    if (!p_handle_improvement_should_continue()) {
                        return false;
                    }
                }
                if (iteration_result->optimal_on_subgraph) {
                    if (!p_optimal_on_subgraph_can_continue()) {
                        return true;
                    }
                } else {
                    if (!p_non_optimal_on_subgraph_can_continue()) {
                        return true;
                    }
                }
                if (!p_should_continue()) {
                    return true;
                }
            }
        } catch (InterruptError&) {
            return std::nullopt;
        }
    }

    /**
     * Whether the last change was a pricing round,
     * i.e., addition of variables.
     */
    bool last_change_was_pricing() const { return m_last_was_pricing; }

    /**
     * Whether the last change was a cutting round,
     * i.e., addition or strengthening of constraints.
     */
    bool last_change_was_cutting() const { return m_last_was_cutting; }

    /**
     * Get the last relaxation value from the clique solver.
     */
    double get_relaxation_value() const noexcept {
        return m_base_solver.get_last_value();
    }

    /**
     * Get the integral upper bound on the current subgraph.
     */
    std::size_t get_subgraph_ub() const noexcept {
        return m_base_solver.get_subgraph_ub();
    }

    /**
     * Get the number of non-trivial vertex values in the current relaxation.
     */
    std::size_t fractional_mes_support_size() const noexcept {
        return m_base_solver.get_fractional_mes_support().size();
    }

  private:
    /**
     * Check if we should continue the LNS
     * subproblem solver MES computation.
     */
    bool p_should_continue() {
        if (m_iterations >= m_iterations_limit) {
            m_recorder->store_event("LNS_MES_STOPPING_DUE_TO_MAX_ITERATIONS",
                                    {{"worker_id", m_worker_id}}, "worker_id");
            return false;
        }
        if (m_iterations_without_improvement >=
            m_iterations_without_improvement_limit)
        {
            m_recorder->store_event(
                "LNS_MES_STOPPING_DUE_TO_MISSING_IMPROVEMENT",
                {{"worker_id", m_worker_id}}, "worker_id");
            return false;
        }
        return true;
    }

    /**
     * Price a sample of vertices from the universe.
     */
    bool p_price_sample() {
        if (m_subproblem.uncovered_universe.size() <=
            PRICING_SAMPLE_SIZE_THRESHOLD)
        {
            return p_price_all();
        }
        std::size_t sample_goal_size = 1000;
        sample_goal_size =
            (std::max)(sample_goal_size,
                       std::size_t(0.1 *
                                   m_subproblem.uncovered_universe.size()));
        if (sample_goal_size >= m_subproblem.uncovered_universe.size()) {
            return p_price_all();
        }
        auto sample = sample_from_range(m_subproblem.uncovered_universe,
                                        sample_goal_size, sammy::rng());
        return m_base_solver.price_vertices(sample.begin(), sample.end()) > 0;
    }

    /**
     * Price all vertices of the universe.
     */
    bool p_price_all() {
        return m_base_solver.price_vertices(
                   m_subproblem.uncovered_universe.begin(),
                   m_subproblem.uncovered_universe.end()) > 0;
    }

    /**
     * Price a sample or, if necessary, all vertices of the universe.
     */
    bool p_price_sample_or_all() {
        if (m_subproblem.uncovered_universe.size() <=
            PRICING_SAMPLE_SIZE_THRESHOLD)
        {
            return p_price_all();
        }
        return p_price_sample() || p_price_all();
    }

    /**
     * Handle the case when we know that we have
     * the optimum MES on the current subgraph.
     * Here, we have to begin by pricing.
     */
    bool p_optimal_on_subgraph_can_continue() {
        m_cut_rounds_without_pricing = 0;
        if (p_price_sample_or_all()) {
            m_last_was_cutting = false;
            m_last_was_pricing = true;
            return true;
        }
        return false;
    }

    /**
     * Handle the case when we don't know if we have the
     * optimum MES on the current subgraph.
     * Here, we start by cuts; if separation fails or
     * every few rounds, we also perform pricing.
     */
    bool p_non_optimal_on_subgraph_can_continue() {
        if (++m_cut_rounds_without_pricing == 6) {
            m_cut_rounds_without_pricing = 0;
            if (p_price_sample()) {
                m_last_was_pricing = true;
                m_last_was_cutting = false;
                return true;
            }
        }
        if (!p_cheap_cut_round()) {
            return p_cheap_cuts_failed();
        }
        m_last_was_pricing = false;
        m_last_was_cutting = true;
        return true;
    }

    /**
     * Run a cheap cutting plane round;
     * return true if a cut was found or
     * an existing constraint strengthened
     * to cut of the current relaxed solution.
     */
    bool p_cheap_cut_round() {
        return m_base_solver.greedy_add_to_cutting_planes() ||
               m_base_solver.greedy_generate_cutting_planes();
    }

    bool p_cheap_cuts_failed() {
        // TODO expensive cuts?
        m_cut_rounds_without_pricing = 0;
        if (p_price_sample_or_all()) {
            m_last_was_pricing = true;
            m_last_was_cutting = false;
            return true;
        }
        return false;
    }

    bool p_handle_improvement_should_continue() {
        std::size_t new_mes_size = m_base_solver.get_best_mes().size();
        m_recorder->store_event(
            "LNS_MES_IMPROVED_BY_CNP",
            {{"mes_size", new_mes_size}, {"worker_id", m_worker_id}},
            "mes_size", "worker_id");
        if (new_mes_size > m_initial_best_global_mes_size) {
            m_portfolio->report_mes(m_base_solver.get_best_mes(),
                                    "Cut & Price on LNS subproblem");
            m_initial_best_global_mes_size = new_mes_size;
        }
        if (new_mes_size >= m_subproblem.removed_configurations.size()) {
            m_recorder->store_event("LNS_MES_BY_CNP_PRECLUDES_IMPROVEMENT",
                                    {{"worker_id", m_worker_id}}, "worker_id");
            return false;
        }
        return true;
    }

    struct IterationResult {
        bool improved;
        bool optimal_on_subgraph;
    };

    std::optional<IterationResult> p_solve_iteration() {
        std::size_t old_mes_size = m_base_solver.get_best_mes().size();
        IterationResult result{false, false};
        ++m_iterations_without_improvement;
        ++m_iterations;
        auto before_solve_full = std::chrono::steady_clock::now();
        auto sfr_result = m_base_solver.solve_full_relaxation();
        auto after_solve_full = std::chrono::steady_clock::now();
        get_global_stats().double_stat_add(
            "MESImprovementSolveFullRelaxationTime",
            seconds_between(before_solve_full, after_solve_full));
        switch (sfr_result) {
        default:
        case SolverState::TIMEOUT_IMPROVEMENT:
        case SolverState::TIMEOUT_NO_IMPROVEMENT:
            return std::nullopt;

        case SolverState::OPTIMUM_ON_SUBGRAPH:
            result.optimal_on_subgraph = true;
            if (old_mes_size < m_base_solver.get_best_mes().size()) {
                result.improved = true;
                m_iterations_without_improvement = 0;
            }
            return result;

        case SolverState::IMPROVEMENT_FOUND:
            result.improved = true;
            m_iterations_without_improvement = 0;
            return result;

        case SolverState::NO_IMPROVEMENT_FOUND:
            return result;
        }
    }

    std::vector<Vertex> p_get_initial_vertices() {
        if (m_subproblem.uncovered_universe.size() <=
            SUB_MES_SOLVER_FULL_GRAPH_THRESHOLD)
        {
            return m_subproblem.uncovered_universe;
        }
        return p_get_initial_vertices_sample(
            SUB_MES_SOLVER_FULL_GRAPH_THRESHOLD);
    }

    std::vector<Vertex> p_get_initial_vertices_sample(std::size_t sample_size) {
        return sample_from_range(m_subproblem.uncovered_universe, sample_size,
                                 sammy::rng());
    }

    static LowerBoundMIPConfig p_get_config() {
        LowerBoundMIPConfig config;
        config.quiet_gurobi = true;
        return config;
    }

    UniverseSubgraph
    p_make_universe_subgraph(LNSSubproblem& restore_subproblem_to) {
        try {
            std::vector<Vertex> initial_vertices = p_get_initial_vertices();
            UniverseSubgraph subgraph{m_clauses, &m_single_thread,
                                      &m_portfolio->get_infeasibility_map(),
                                      std::move(initial_vertices)};
            subgraph.extend_matrix_by_propagation(true);
            return subgraph;
        } catch (InterruptError&) {
            restore_subproblem_to = std::move(m_subproblem);
            throw;
        }
    }

    PortfolioSolver* m_portfolio;
    EventRecorder* m_recorder;
    std::size_t m_worker_id;
    std::size_t m_initial_best_global_mes_size;
    ThreadGroup<void> m_single_thread;
    ClauseDB* m_clauses;
    LNSSubproblem m_subproblem;
    UniverseSubgraph m_subgraph;
    GurobiCliqueSolverG2 m_base_solver;
    std::size_t m_cut_rounds_without_pricing = 0;
    std::size_t m_iterations_without_improvement = 0;
    std::size_t m_iterations = 0;
    std::size_t m_iterations_without_improvement_limit;
    std::size_t m_iterations_limit;
    bool m_last_was_cutting{false};
    bool m_last_was_pricing{false};
};

template <typename WrappedSubproblemSolver> class SubproblemSolverWithMES {
  public:
    static std::string name() {
        return std::string("MESImprovement<") +
               WrappedSubproblemSolver::name() + ">";
    }

    std::string strategy_name() const {
        if (!m_sub_solver) {
            return "mes_only";
        } else {
            return m_sub_solver->strategy_name();
        }
    }

    SubproblemSolverWithMES(PortfolioSolver* portfolio,
                            LNSSubproblem&& subproblem, SharedDBPropagator prop,
                            EventRecorder* recorder, std::size_t worker_id)
        : m_portfolio(portfolio), m_recorder(recorder), m_worker_id(worker_id),
          m_time_before_construction(std::chrono::steady_clock::now()),
          m_mes_solver(std::make_unique<SubproblemMESSolver>(
              portfolio, std::move(subproblem), std::move(prop), recorder,
              worker_id)),
          m_sub_solver() {
        get_global_stats().double_stat_add(
            "MESImprovementConstructionTime",
            seconds_between(m_time_before_construction,
                            std::chrono::steady_clock::now()));
    }

    void abort() {
        // we only need the lock when:
        //  - checking the pointers from control thread
        //  - changing the pointers in the worker thread
        std::unique_lock l{m_abort_mutex};
        if (m_mes_solver) {
            m_mes_solver->abort();
        } else {
            m_sub_solver->abort();
        }
    }

    std::optional<bool> solve() {
        auto time_before = std::chrono::steady_clock::now();
        std::optional<bool> mes_result = m_mes_solver->solve();
        auto time_after = std::chrono::steady_clock::now();
        get_global_stats().double_stat_add(
            "MESImprovementTime", seconds_between(time_before, time_after));
        if (!mes_result || !*mes_result) {
            get_global_stats().int_stat_add("MESImprovementDidResolveLNS", 1);
            return mes_result;
        } else {
            get_global_stats().int_stat_add("MESImprovementDidNotResolveLNS",
                                            1);
        }
        LNSSubproblem subproblem = m_mes_solver->move_out_subproblem();
        m_original_mes = std::move(subproblem.mutually_exclusive_set);
        subproblem.mutually_exclusive_set = m_mes_solver->mes_vertices();
        try {
            p_switch_phase(std::move(subproblem));
        } catch (InterruptError&) {
            subproblem.mutually_exclusive_set = std::move(m_original_mes);
            m_mes_solver->return_subproblem(std::move(subproblem));
            return std::nullopt;
        }
        auto time_after_phase_switch = std::chrono::steady_clock::now();
        get_global_stats().double_stat_add(
            "CoreSolverConstructionTime",
            seconds_between(time_after, time_after_phase_switch));
        auto solver_result = m_sub_solver->solve();
        get_global_stats().double_stat_add(
            "CoreSolverSolveTime",
            seconds_between(time_after_phase_switch,
                            std::chrono::steady_clock::now()));
        return solver_result;
    }

    LNSSubproblem move_out_subproblem() noexcept {
        if (m_mes_solver) {
            return m_mes_solver->move_out_subproblem();
        } else {
            LNSSubproblem subproblem = m_sub_solver->move_out_subproblem();
            subproblem.mutually_exclusive_set = std::move(m_original_mes);
            return subproblem;
        }
    }

    const std::vector<Vertex>& mes_vertices() const {
        if (m_mes_solver) {
            return m_mes_solver->mes_vertices();
        } else {
            return m_sub_solver->mes_vertices();
        }
    }

    const std::vector<DynamicBitset>& get_solution() const {
        if (m_mes_solver) {
            return m_mes_solver->get_solution();
        } else {
            return m_sub_solver->get_solution();
        }
    }

  private:
    void p_switch_phase(LNSSubproblem&& subproblem) {
        assert(m_mes_solver);
        assert(!m_sub_solver);
        auto new_sub_solver = std::make_unique<WrappedSubproblemSolver>(
            m_portfolio, std::move(subproblem),
            SharedDBPropagator(&m_portfolio->get_clauses()), m_recorder,
            m_worker_id);
        std::unique_ptr<SubproblemMESSolver>
            tmp; // let destructor run outside of locked region
        {
            std::unique_lock l{m_abort_mutex};
            m_sub_solver = std::move(new_sub_solver);
            tmp = std::move(m_mes_solver);
        }
        assert(m_sub_solver);
        assert(!m_mes_solver);
    }

    PortfolioSolver* m_portfolio;
    EventRecorder* m_recorder;
    std::size_t m_worker_id;
    std::chrono::steady_clock::time_point m_time_before_construction;

    std::mutex m_abort_mutex;
    std::unique_ptr<SubproblemMESSolver> m_mes_solver;
    std::vector<Vertex> m_original_mes;
    std::unique_ptr<WrappedSubproblemSolver> m_sub_solver;
};

} // namespace sammy

#endif
