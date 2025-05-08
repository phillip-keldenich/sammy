#ifndef SAMMY_BARRAGE_WORKER_LNS_H_INCLUDED_
#define SAMMY_BARRAGE_WORKER_LNS_H_INCLUDED_

#include "algorithm_ex.h"
#include "barrage.h"
#include "barrage_worker_with_core.h"
#include "clique_sat_dsatur.h"
#include "gurobi_clique_solver_g2.h"
#include "implied_vertices.h"
#include "lns_destroy.h"
#include "satdsatur_colorer.h"
#include "thread_interrupt.h"

namespace sammy {

static constexpr double RANDOM_DESTRUCTION_PROB = 0.3;
static constexpr double RANDOM_WITH_TABLE_PROB = 0.2;
static constexpr double RANDOM_AND_LEAST_UNIQUE_DESTRUCTION_PROB = 0.5;

/**
 * LNS solver core that wraps different subproblem solvers
 * with unified mechanisms for destroying and reporting improvements,
 * aborting/timeouts.
 */
template <typename SubproblemSolverType> class SubproblemLNSSolverCore {
  public:
    using SubproblemSolver = SubproblemSolverType;

    static std::unique_ptr<SubproblemLNSSolverCore>
    factory(PortfolioSolver* solver,
            PortfolioElementWithCore<SubproblemLNSSolverCore>* element) {
        return std::make_unique<SubproblemLNSSolverCore>(solver, element);
    }

    SubproblemLNSSolverCore(
        PortfolioSolver* solver,
        PortfolioElementWithCore<SubproblemLNSSolverCore>* element)
        : m_portfolio(solver), m_element_mutex(&element->get_mutex()),
          m_termination_flag(get_interrupt_flag_ptr()), m_interrupt_flag(false),
          m_local_recorder(element->get_mutable_recorder()),
          m_source("LNSCore<" + SubproblemSolverType::name() + ">"),
          m_clauses(&solver->get_clauses()), m_propagator(m_clauses),
          m_destroy(m_local_recorder, m_worker_counter++,
                    &solver->implied_cache(), solver->get_best_solution(),
                    solver->get_best_mes(), m_propagator),
          m_destroy_seen_solution(m_destroy.total_num_configs()) {
        set_interrupt_flag_ptr(&m_interrupt_flag);
        if (m_termination_flag->load()) {
            throw InterruptError();
        }
    }

    /**
     * Called by the PortfolioElementWithCore
     * if it finds the termination flag to be set.
     */
    void termination_flag_set() {
        assert(m_termination_flag->load());
        p_interrupt();
    }

    void interrupt_if_necessary(const InterruptionCheckInfo& info) {
        if (m_termination_flag->load()) {
            p_interrupt();
        } else if (info.best_upper_bound < m_destroy_seen_solution) {
            m_destroy_seen_solution = info.best_upper_bound;
            p_interrupt();
        }
    }

    void main() {
        while (!m_termination_flag->load()) {
            auto phase_result = p_destroy_phase();
            if (phase_result == DestroyPhaseResult::GLOBAL_OPT) {
                break;
            } else if (phase_result == DestroyPhaseResult::INTERRUPTED) {
                continue;
            }
            LNSSubproblem subproblem = m_destroy.move_out_subproblem();
            if (m_portfolio->subproblem_reporting_enabled()) {
                m_portfolio->report_subproblem(
                    subproblem, m_destroy.get_partial(),
                    m_destroy.best_global_mes().size(),
                    m_portfolio->get_best_lower_bound(), m_source.c_str());
            }
            auto before = std::chrono::steady_clock::now();
            try {
                p_create_subsolver(std::move(subproblem));
            } catch (InterruptError&) {
                m_destroy.return_subproblem_on_abort(std::move(subproblem));
                continue;
            }
            auto res = m_solver->solve();
            double time_taken =
                seconds_between(before, std::chrono::steady_clock::now());
            LNSSubproblem tmp = m_solver->move_out_subproblem();
            if (!res) {
                m_portfolio->lns_report_aborted(
                    tmp.removed_configurations.size(), time_taken,
                    m_solver->strategy_name());
                m_destroy.return_subproblem_on_abort(std::move(tmp));
                continue;
            }
            if (!*res) {
                auto old_lb = m_portfolio->get_best_lower_bound();
                if (old_lb < tmp.removed_configurations.size()) {
                    m_portfolio->report_lower_bound(
                        tmp.removed_configurations.size(),
                        tmp.uncovered_universe, m_source.c_str());
                }
                m_portfolio->lns_report_failure(
                    tmp.removed_configurations.size(), time_taken,
                    m_solver->strategy_name());
                m_destroy.improvement_impossible(std::move(tmp),
                                                 m_solver->mes_vertices());
            } else {
                m_portfolio->lns_report_success(
                    tmp.removed_configurations.size(), time_taken,
                    m_solver->strategy_name());
                const auto& improvement = m_solver->get_solution();
                PartialSolution improved =
                    m_destroy.improve_destroyed(improvement);
                if (m_portfolio->report_solution(improved, m_source.c_str())) {
                    m_destroy.return_subproblem_on_success(std::move(tmp),
                                                           std::move(improved));
                    std::size_t new_size = m_destroy.total_num_configs();
                    {
                        std::unique_lock l{*m_element_mutex};
                        if (new_size < m_destroy_seen_solution) {
                            m_destroy_seen_solution = new_size;
                        }
                    }
                } else {
                    m_destroy.return_subproblem_on_abort(std::move(tmp));
                }
            }
        }
        m_local_recorder->store_event("LNS_WORKER_EXITING",
                                      {{"worker_id", m_destroy.worker_index()}},
                                      "worker_id");
    }

  private:
    enum class DestroyPhaseResult { GLOBAL_OPT, INTERRUPTED, SUBPROBLEM_FOUND };

    DestroyPhaseResult p_destroy_phase() {
        {
            std::unique_lock l{*m_element_mutex};
            if (m_solver) {
                m_solver.reset();
            }
        }
        try {
            p_update_destroy_if_needed();
            if (get_and_clear_interrupt_flag() || m_termination_flag->load()) {
                return DestroyPhaseResult::INTERRUPTED;
            }
            std::size_t num_removed = m_destroy.destroy(
                m_portfolio->lns_select_removed_class_count());
            if (num_removed == 0) {
                m_portfolio->report_lower_bound(
                    m_destroy.total_num_configs(),
                    m_portfolio->get_infeasibility_map().collect_vertices(),
                    "LNS whole solution destruction");
                return DestroyPhaseResult::GLOBAL_OPT;
            }
            return DestroyPhaseResult::SUBPROBLEM_FOUND;
        } catch (InterruptError&) {
            return DestroyPhaseResult::INTERRUPTED;
        }
    }

    void p_update_destroy_if_needed() {
        if (p_mes_update_looks_necessary()) {
            m_destroy.update_global_mes_if_better(m_portfolio->get_best_mes());
        }
        if (!p_update_looks_necessary())
            return;
        PartialSolution new_full = m_portfolio->get_best_solution();
        {
            std::unique_lock l{*m_element_mutex};
            m_destroy_seen_solution = new_full.size();
        }
        m_destroy.update_full_solution_if_better(std::move(new_full));
    }

    bool p_mes_update_looks_necessary() const {
        return m_portfolio->get_best_mes_size() >
               m_destroy.best_global_mes().size();
    }

    bool p_update_looks_necessary() const {
        std::unique_lock l{*m_element_mutex};
        return m_portfolio->get_best_solution_size() <
               m_destroy.total_num_configs();
    }

    void p_interrupt() {
        m_interrupt_flag.store(true);
        if (m_solver) {
            m_solver->abort();
        }
    }

    void p_create_subsolver(LNSSubproblem&& subproblem) {
        SharedDBPropagator prop{m_clauses};
        // relies on the subsolver guarantee that,
        // on interrupt exceptions, subproblem remains usable
        auto subsolver = std::make_unique<SubproblemSolver>(
            m_portfolio, std::move(subproblem), std::move(prop),
            m_local_recorder, m_destroy.worker_index());
        std::unique_lock l{*m_element_mutex};
        m_solver = std::move(subsolver);
    }

    static std::atomic<std::size_t> m_worker_counter;
    PortfolioSolver* m_portfolio;
    std::mutex* m_element_mutex;
    std::atomic<bool>* m_termination_flag;
    std::atomic<bool> m_interrupt_flag;
    EventRecorder* m_local_recorder;
    std::string m_source;
    ClauseDB* m_clauses;
    SharedDBPropagator m_propagator;
    LNSDestroy m_destroy;
    std::size_t m_destroy_seen_solution;
    std::unique_ptr<SubproblemSolver> m_solver;
};

template <typename S>
std::atomic<std::size_t> SubproblemLNSSolverCore<S>::m_worker_counter{0};

} // namespace sammy

#endif
