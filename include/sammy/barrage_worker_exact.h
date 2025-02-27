#ifndef SAMMY_BARRAGE_WORKER_EXACT_H_INCLUDED_
#define SAMMY_BARRAGE_WORKER_EXACT_H_INCLUDED_

#include "barrage.h"
#include "cadical_solver.h"
#include "clique_sat_dsatur.h"
#include "implied_vertices.h"
#include "lingeling_solver.h"

namespace sammy {

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
