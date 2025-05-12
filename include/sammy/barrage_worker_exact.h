#ifndef SAMMY_BARRAGE_WORKER_EXACT_H_INCLUDED_
#define SAMMY_BARRAGE_WORKER_EXACT_H_INCLUDED_

#include "barrage.h"
#include "clique_sat_dsatur.h"
#include "implied_vertices.h"
#include "sat_lns.h"

#include <sammy/cadical_solver.h>
#include <sammy/kissat_solver.h>

namespace sammy {

/**
 * Exact solver based on encoding the problem as a sequence of
 * SAT instances, considering all non-eliminated interactions
 * at once.
 */
class FixedSatExactSolverCore {
  public:
    static std::unique_ptr<FixedSatExactSolverCore>
    factory(PortfolioSolver* solver,
            PortfolioElementWithCore<FixedSatExactSolverCore>* element) {
        return std::make_unique<FixedSatExactSolverCore>(
            solver, &element->get_mutex(), element->get_mutable_recorder());
    }

    using SatSolver = KissatSolver;
    using ExactQuerySolver = FixedBoundSATSolver<SatSolver>;

    explicit FixedSatExactSolverCore(PortfolioSolver* solver, std::mutex* mut,
                                     EventRecorder* local_recorder)
        : m_element_mutex(mut), m_solver(solver),
          m_clauses(solver->get_clauses()), m_propagator(&m_clauses),
          m_local_recorder(*local_recorder),
          m_source(std::string("ExactFixedSat<") + SatSolver::name() + ">"),
          m_icache(p_check_cache()) {}

    /**
     * Called by the PortfolioElementWithCore
     * if it finds the termination flag to be set.
     */
    void termination_flag_set() { p_interrupt(); }

    void interrupt_if_necessary(const InterruptionCheckInfo& info) {
        if (!m_working_on_value) {
            return;
        }
        if (*m_working_on_value < info.best_lower_bound ||
            *m_working_on_value >= info.best_upper_bound)
        {
            p_interrupt();
        }
    }

    void main() {
        for (;;) {
            if (get_and_clear_interrupt_flag()) {
                break;
            }
            auto mes = m_solver->get_best_mes();
            std::size_t ub = m_solver->get_best_solution_size();
            std::size_t lb =
                (std::max)(mes.size(), m_solver->get_best_lower_bound());
            if (ub <= lb) {
                break;
            }
            std::size_t query;
            if (ub - lb > 3) {
                query = std::size_t(0.5 * (ub + lb));
            } else {
                query = lb;
            }
            auto eqs = std::make_unique<ExactQuerySolver>(
                &m_propagator, &m_icache->get_reduced_universe(), mes, query);
            {
                std::unique_lock l{*m_element_mutex};
                m_working_on_value = query;
                m_solver_instance = std::move(eqs);
            }
            auto query_result = m_solver_instance->solve();
            {
                std::unique_lock l{*m_element_mutex};
                m_working_on_value.reset();
                eqs = std::move(m_solver_instance);
            }
            if (!query_result) {
                continue;
            }
            if (*query_result) {
                m_local_recorder.store_event("EXACT_SAT_FOUND_SOLUTION",
                                             {{"query", query}}, "query");
                auto& inf = m_solver->get_infeasibility_map();
                m_solver->report_solution(eqs->get_partial(&inf),
                                          m_source.c_str());
            } else {
                m_solver->report_lower_bound(query + 1,
                                             m_icache->get_reduced_universe(),
                                             m_source.c_str());
            }
        }
        m_local_recorder.store_event("EXACT_SAT_EXITING",
                                     {{"source", m_source}}, "source");
    }

  private:
    void p_interrupt() {
        if (m_solver_instance) {
            m_solver_instance->abort();
        }
    }

    const ImpliedVertexCache* p_check_cache() {
        const ImpliedVertexCache& icache = m_solver->implied_cache();
        if (!icache.have_reduced_universe()) {
            throw std::logic_error(
                "Did not reduce universe in exact portfolio element!");
        }
        return &icache;
    }

    std::mutex* m_element_mutex;
    PortfolioSolver* m_solver;
    ClauseDB& m_clauses;
    SharedDBPropagator m_propagator;
    EventRecorder& m_local_recorder;
    std::string m_source;
    const ImpliedVertexCache* m_icache;

    // only change state/identity under element mutex
    std::optional<std::size_t> m_working_on_value;
    std::unique_ptr<ExactQuerySolver> m_solver_instance;
};

/**
 * Exact solver based on CliqueSatDSatur.
 */
class CliqueSatDSaturExactSolverCore {
  public:
    static std::unique_ptr<CliqueSatDSaturExactSolverCore> factory(
        PortfolioSolver* solver,
        PortfolioElementWithCore<CliqueSatDSaturExactSolverCore>* our_element) {
        return std::make_unique<CliqueSatDSaturExactSolverCore>(
            solver, our_element->get_mutable_recorder());
    }

    using IncrementalSolver = CadicalSolver;
    using CSDSSolver = CliqueSatDSaturSolver<IncrementalSolver>;

    explicit CliqueSatDSaturExactSolverCore(PortfolioSolver* solver,
                                            EventRecorder* local_recorder)
        : m_solver(solver), m_clauses(solver->get_clauses()),
          m_local_recorder(*local_recorder), m_icache(p_check_cache()),
          m_exact(m_icache->get_reduced_universe(),
                  &solver->get_infeasibility_map(), solver->get_clauses(),
                  solver->get_best_mes(), solver->get_best_mes().size(),
                  solver->get_best_lower_bound(),
                  solver->get_best_solution().assignments()) {
        m_exact.set_clique_candidate_callback([this]() -> std::vector<Vertex> {
            return m_solver->get_best_mes();
        });
        m_exact.set_lower_bound_callback(
            [this](std::size_t lb, const std::vector<Vertex>& subgraph) {
                if (lb == subgraph.size()) {
                    m_solver->report_mes(subgraph,
                                         "C&P/SATDSatur Exact Solver");
                } else {
                    m_solver->report_lower_bound(lb, subgraph,
                                                 "C&P/SATDSatur Exact Solver");
                }
            });
        m_exact.set_event_recorder(&m_local_recorder);
    }

    /**
     * Aside from handling the termination flag,
     * we do not need to check for other interrupt reasons.
     */
    void interrupt_if_necessary(const InterruptionCheckInfo& /*info*/) {}

    /**
     * Called by the PortfolioElementWithCore
     * if it finds the termination flag to be set.
     */
    void termination_flag_set() { m_exact.abort(); }

    void main() {
        auto status = m_exact.solve();
        if (status == CSDSSolver::SolveResult::ABORTED) {
            m_local_recorder.store_event("EXACT_ABORTED");
            return;
        }
        // the lower bound is reported via the callback;
        // the solution is reported here
        if (status == CSDSSolver::SolveResult::IMPROVED_SOLUTION) {
            auto solution = m_exact.get_partial_solution();
            m_solver->report_solution(solution, "C&P/SATDSatur Exact Solver");
        }
        m_local_recorder.store_event("EXACT_DONE");
    }

  private:
    const ImpliedVertexCache* p_check_cache() {
        const ImpliedVertexCache& icache = m_solver->implied_cache();
        if (!icache.have_reduced_universe()) {
            throw std::logic_error(
                "Did not reduce universe in exact portfolio element!");
        }
        return &icache;
    }

    PortfolioSolver* m_solver;
    ClauseDB& m_clauses;
    EventRecorder& m_local_recorder;
    const ImpliedVertexCache* m_icache;
    CSDSSolver m_exact;
};

} // namespace sammy

#endif
