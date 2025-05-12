#ifndef SAMMY_VARIANT_SUBPROBLEM_SOLVER_H_INCLUDED_
#define SAMMY_VARIANT_SUBPROBLEM_SOLVER_H_INCLUDED_

#include "incremental_sat_lns.h"
#include "sat_dsatur.h"
#include "sat_lns.h"

#include <sammy/cadical_solver.h>
#include <sammy/kissat_solver.h>
#include <variant>

namespace sammy {

/**
 * Solver type that can dynamically choose a
 * strategy (fixed sat, incremental sat, satdsatur)
 * and subsolver to solve a given LNS subproblem.
 */
class VariantSubproblemSolver {
    using FixedSATCadical = FixedMESSATImprovementSolver<CadicalSolver>;
    using FixedSATKissat = FixedMESSATImprovementSolver<KissatSolver>;
    using IncSATCadical =
        FixedMESIncrementalSATImprovementSolver<CadicalSolver>;
    using SatDSaturCadical = FixedMESSatDSaturSolver<CadicalSolver>;
    using Variant = std::variant<FixedSATCadical, FixedSATKissat, IncSATCadical,
                                 SatDSaturCadical>;

    std::string m_strategy;
    std::size_t m_worker_id;
    std::optional<Variant> m_wrapped;
    LNSSubproblem m_subproblem;

  public:
    static std::string name() { return "Variant"; }

    std::string strategy_name() const { return m_strategy; }

    VariantSubproblemSolver(PortfolioSolver* portfolio,
                            LNSSubproblem&& subproblem, SharedDBPropagator prop,
                            EventRecorder* recorder, std::size_t worker_id)
        : m_strategy(portfolio->lns_select_strategy()), m_worker_id(worker_id),
          m_subproblem(std::move(subproblem)) {
        try {
            if (m_strategy == "satdsatur|cadical") {
                p_emplace_solver<SatDSaturCadical>(portfolio, std::move(prop),
                                                   recorder);
            } else if (m_strategy == "fixed_sat|cadical") {
                p_emplace_solver<FixedSATCadical>(portfolio, std::move(prop),
                                                  recorder);
            } else if (m_strategy == "fixed_sat|kissat") {
                p_emplace_solver<FixedSATKissat>(portfolio, std::move(prop),
                                                 recorder);
            } else if (m_strategy == "inc_sat|cadical") {
                p_emplace_solver<IncSATCadical>(portfolio, std::move(prop),
                                                recorder);
            } else {
                throw std::logic_error("Unknown LNS strategy: " + m_strategy);
            }
        } catch (...) {
            subproblem = move_out_subproblem();
            throw;
        }
    }

    LNSSubproblem move_out_subproblem() noexcept {
        if (m_wrapped) {
            return std::visit(
                [&](auto& solver) -> LNSSubproblem {
                    return solver.move_out_subproblem();
                },
                *m_wrapped);
        } else {
            return std::move(m_subproblem);
        }
    }

    const std::vector<Vertex>& mes_vertices() const {
        assert(m_wrapped);
        return std::visit(
            [&](auto& solver) -> const std::vector<Vertex>& {
                return solver.mes_vertices();
            },
            *m_wrapped);
    }

    void abort() {
        assert(m_wrapped);
        std::visit([&](auto& solver) { solver.abort(); }, *m_wrapped);
    }

    const std::vector<DynamicBitset>& get_solution() const {
        assert(m_wrapped);
        return std::visit(
            [&](auto& solver) -> const std::vector<DynamicBitset>& {
                return solver.get_solution();
            },
            *m_wrapped);
    }

    std::optional<bool> solve() {
        assert(m_wrapped);
        return std::visit(
            [&](auto& solver) -> std::optional<bool> { return solver.solve(); },
            *m_wrapped);
    }

  private:
    template <typename Solver>
    void p_emplace_solver(PortfolioSolver* portfolio, SharedDBPropagator&& prop,
                          EventRecorder* recorder) {
        m_wrapped.emplace(std::in_place_type<Solver>, portfolio,
                          std::move(m_subproblem), std::move(prop), recorder,
                          m_worker_id);
    }
};

} // namespace sammy

#endif
