#ifndef SAMMY_INCREMENTAL_SAT_LNS_H_INCLUDED_
#define SAMMY_INCREMENTAL_SAT_LNS_H_INCLUDED_

#include "clause_db.h"
#include "shared_db_propagator.h"
#include "lazy_g2_adjacency_matrix.h"
#include "output.h"
#include "barrage_lns_subproblem.h"
#include "barrage.h"
#include <variant>
#include <algorithm>

namespace sammy {

/**
 * Generalized subproblem solver interface:
 *  - creation from LNSSubproblem, Propagator and EventRecorder*,
 *  - interruptible solve method,
 *  - solution extraction method,
 *  - abort method,
 *  - MES extraction method.
 */

template<typename IncrementalSatSolver>
class FixedMESIncrementalSATImprovementSolver {
  public:
    using SatSolver = IncrementalSatSolver;
    using Lit = typename SatSolver::Lit;
    using SLit = sammy::Lit;
    using LitOrVal = std::variant<Lit, bool>;

    static std::string name() {
        return std::string("IncrementalSAT<") + SatSolver::name() + ">";
    }

  private:
    // basic incremental SAT-solver backend
    SatSolver m_solver;
    // subproblem information
    LNSSubproblem m_subproblem;
    // event recorder
    EventRecorder* m_recorder;
    // worker id for events
    std::size_t m_worker_id;
    // propagator for the SAT model
    SharedDBPropagator m_propagator;
    // flag to set when, during construction, we find infeasibility
    bool m_infeasible_by_construction;
    // contiguous array of all second literals of vertices
    std::vector<SLit> m_literal_partners_array;
    // bitmap that tracks, for each entry in m_literal_partners_array, whether the corresponding vertex is covered.
    DynamicBitset m_covered_literal_partners;
    // tracks which literal partner entries are explicitly covered by clauses/variables
    DynamicBitset m_explicitized_literal_partners;
    // how many vertices are explicitly covered?
    std::size_t m_num_explicitly_covered = 0;
    // for each first literal of a vertex, a [begin_index, end_index) of its second literals
    // as they appear in m_literal_partners_array
    std::vector<std::pair<std::size_t, std::size_t>> m_literal_partners_of;
    // (outer) index k: configuration index [0, allowed_configurations());
    // inner index: feature index [0, propagator().db().num_vars())
    std::vector<std::vector<LitOrVal>> m_config_vars;
    // current solution: for each configuration, assignment of variables
    std::vector<DynamicBitset> m_current_solution;
    // buffer for solver literals/clauses
    std::vector<Lit> m_buffer;
    // buffer for selecting uncovered vertices to be explicitized
    std::vector<std::pair<Vertex, std::size_t>> m_selected_uncovered;
    // set by abort method in addition to aborting the SAT solver itself
    std::atomic<bool> m_abort{false};

  public:
    FixedMESIncrementalSATImprovementSolver(PortfolioSolver* portfolio, LNSSubproblem&& subproblem, SharedDBPropagator prop,
                                            EventRecorder* recorder, std::size_t worker_id) :
        m_solver(),
        m_subproblem(std::move(subproblem)),
        m_recorder(recorder),
        m_worker_id(worker_id),
        m_propagator(std::move(prop)),
        m_infeasible_by_construction(allowed_configurations() < mes_size())
    {
        try {
            if(m_infeasible_by_construction) return;
            p_prepare_vertex_set();
            throw_if_interrupted();
            p_make_literal_partners();
            throw_if_interrupted();
            p_make_config_vars();
            if(m_infeasible_by_construction) return;
            throw_if_interrupted();
            p_make_vertical_clauses();
            throw_if_interrupted();
        } catch(...) {
            subproblem = std::move(m_subproblem);
            throw;
        }
    }

    std::size_t allowed_configurations() const noexcept {
        return m_subproblem.removed_configurations.size() - 1;
    }

    std::size_t mes_size() const noexcept {
        return m_subproblem.mutually_exclusive_set.size();
    }

    std::size_t num_uncovered() const noexcept {
        return m_subproblem.uncovered_universe.size();
    }

    const std::vector<Vertex>& mes_vertices() const noexcept {
        return m_subproblem.mutually_exclusive_set;
    }

    std::optional<bool> solve() {
        if(m_infeasible_by_construction) {
            p_store_event("INCREMENTAL_SAT_INFEASIBLE_BY_CONSTRUCTION");
            return false;
        }
        for(;;) {
            OutputObject event_data{
                {"universe_size", num_uncovered()},
                {"num_removed", m_subproblem.removed_configurations.size()},
                {"mes_size", mes_size()},
                {"num_explicitly_covered", m_num_explicitly_covered}
            };
            p_store_event("INCREMENTAL_SAT_LNS_BEGIN_SAT_SOLVE", event_data, 
                          "universe_size", "num_removed", "mes_size", "num_explicitly_covered");
            std::optional<bool> res = p_solve_iteration();
            if(!res) {
                p_store_event("INCREMENTAL_SAT_LNS_SOLUTION_ABORTED", event_data,
                              "universe_size", "num_removed", "mes_size", "num_explicitly_covered");
                return std::nullopt;
            }
            if(!*res) {
                p_store_event("INCREMENTAL_SAT_LNS_SOLUTION_WAS_OPTIMAL", event_data, 
                              "universe_size", "num_removed", "mes_size", "num_explicitly_covered");
                return false;
            }
            p_store_event("INCREMENTAL_SAT_LNS_SATISFIABLE", event_data,
                          "universe_size", "num_removed", "mes_size", "num_explicitly_covered");
            p_extract_solution();
            p_identify_covered();
            p_sample_uncovered();
            if(m_selected_uncovered.empty()) {
                p_store_event("INCREMENTAL_SAT_LNS_SOLUTION_IMPROVED", event_data,
                              "universe_size", "num_removed", "mes_size", "num_explicitly_covered");
                return true;
            }
            for(auto entry : m_selected_uncovered) {
                p_make_explicit_unfixed(entry.second, entry.first);
            }
            if(m_num_explicitly_covered > 0.33 * m_subproblem.uncovered_universe.size()) {
                p_make_all_explicit();
            } else if(m_selected_uncovered.size() < 0.025 * m_subproblem.uncovered_universe.size()) {
                p_sample_more_explicit();
            }
            if(m_infeasible_by_construction) {
                p_store_event("INCREMENTAL_SAT_LNS_SOLUTION_WAS_OPTIMAL", event_data, 
                              "universe_size", "num_removed", "mes_size", "num_explicitly_covered");
                return false;
            }
        }
    }

    LNSSubproblem move_out_subproblem() noexcept {
        return std::move(m_subproblem);
    }

    const std::vector<DynamicBitset>& get_solution() const {
        return m_current_solution;
    }

    void abort() {
        m_abort.store(true);
        m_solver.terminate();
    }

  private:
    void p_prepare_vertex_set() {
        if(m_subproblem.uncovered_universe.empty()) {
            throw std::logic_error("Improvement algorithm run without uncovered vertices!");
        }
        auto canonical = [] (Vertex v) -> Vertex {
            return {std::min(v.first, v.second), std::max(v.first, v.second)};
        };
        std::transform(m_subproblem.uncovered_universe.begin(), m_subproblem.uncovered_universe.end(),
                       m_subproblem.uncovered_universe.begin(), canonical);
        std::sort(m_subproblem.uncovered_universe.begin(), m_subproblem.uncovered_universe.end());
        std::transform(m_subproblem.mutually_exclusive_set.begin(), m_subproblem.mutually_exclusive_set.end(),
                       m_subproblem.mutually_exclusive_set.begin(), canonical);
        std::sort(m_subproblem.mutually_exclusive_set.begin(), m_subproblem.mutually_exclusive_set.end());
    }

    void p_make_literal_partners() {
        const auto& vertices = m_subproblem.uncovered_universe;
        const std::size_t nconcrete = m_subproblem.num_concrete;
        const std::size_t nconclit = 2 * nconcrete;
        m_literal_partners_array.reserve(vertices.size());
        m_literal_partners_of.reserve(nconclit);
        SLit last_first = 0;
        std::size_t current_begin = 0, current_size = 0;
        for(Vertex v : vertices) {
            if(v.first != last_first) {
                m_literal_partners_of.emplace_back(current_begin, current_size);
                while(++last_first < v.first) {
                    m_literal_partners_of.emplace_back(current_size, current_size);
                }
                current_begin = current_size;
            }
            m_literal_partners_array.push_back(v.second);
            ++current_size;
        }
        m_literal_partners_of.emplace_back(current_begin, current_size);
        while(++last_first < nconclit) {
            m_literal_partners_of.emplace_back(current_size, current_size);
        }
        m_explicitized_literal_partners.assign(m_literal_partners_array.size(), false);
        m_covered_literal_partners.assign(m_literal_partners_array.size(), false);
    }

    std::size_t p_find_index(Vertex v) {
        auto p_of = m_literal_partners_of[v.first];
        auto p_of_begin = p_of.first;
        auto p_of_end = p_of.second;
        auto it = m_literal_partners_array.begin() + p_of_begin;
        auto end = m_literal_partners_array.begin() + p_of_end;
        std::size_t result(std::lower_bound(it, end, v.second) - m_literal_partners_array.begin());
        assert(result < p_of_end);
        assert(m_literal_partners_array[result] == v.second);
        return result;
    }

    void p_make_config_vars() {
        const auto allowed = allowed_configurations();
        m_config_vars.reserve(allowed);
        for(Vertex clique_vertex : mes_vertices()) {
            p_create_config_with_clique(clique_vertex);
            std::size_t index = p_find_index(clique_vertex);
            m_explicitized_literal_partners[index].set();
            ++m_num_explicitly_covered;
        }
        while(m_config_vars.size() < allowed) {
            p_create_config_with_clique({NIL, NIL});
        }
    }

    void p_create_config_with_clique(Vertex clique_vertex) {
        const auto nall = m_propagator.db().num_vars();
        m_config_vars.emplace_back();
        auto& config = m_config_vars.back();
        config.reserve(nall);
        m_propagator.reset_to_zero();
        if(clique_vertex.first != clique_vertex.second && push_vertex(m_propagator, clique_vertex) < 0) {
            throw std::logic_error("Infeasible interaction in MES!");
        }
        for(SLit var = 0, n = nall; var < n; ++var) {
            SLit pos = lit::positive_lit(var);
            if(m_propagator.is_true(pos)) {
                // true
                config.emplace_back(std::in_place_type<bool>, true);
            } else if(m_propagator.is_false(pos)) {
                // false
                config.emplace_back(std::in_place_type<bool>, false);
            } else {
                // open
                config.emplace_back(std::in_place_type<Lit>, m_solver.new_var());
            }
        }
        m_propagator.reset_to_zero();
        p_config_insert_binaries();
        p_config_insert_long_clauses();
        throw_if_interrupted();
    }

    LitOrVal p_config_lit_for(const std::vector<LitOrVal>& config, SLit sl) const {
        auto var = lit::var(sl);
        bool is_pos = !lit::negative(sl);
        auto& entry = config[var];
        return std::visit(overloaded{
            [&] (bool b) -> LitOrVal {
                return LitOrVal{std::in_place_type<bool>, is_pos ? b : !b};
            },
            [&] (Lit l) -> LitOrVal {
                return LitOrVal{std::in_place_type<Lit>, is_pos ? l : -l};
            }
        }, entry);
    }

    void p_config_insert_binaries() {
        const auto& cdb = m_propagator.db();
        const auto& config = m_config_vars.back();
        for(auto bin_clause : cdb.binary_clauses()) {
            SLit l1 = bin_clause.first, l2 = bin_clause.second;
            LitOrVal v1 = p_config_lit_for(config, l1);
            LitOrVal v2 = p_config_lit_for(config, l2);
            if(std::holds_alternative<bool>(v1) || std::holds_alternative<bool>(v2)) {
                if(std::holds_alternative<bool>(v1) && std::holds_alternative<bool>(v2)) {
                    if(!*std::get_if<bool>(&v1) && !*std::get_if<bool>(&v2)) {
                        // both false, construction is infeasible;
                        // this should not actually happen.
                        m_infeasible_by_construction = true;
                        break;
                    }
                }
                // otherwise, one literal is true (either by propagation or directly);
                // need not add clause.
                continue;
            } else {
                m_solver.add_short_clause(*std::get_if<Lit>(&v1), *std::get_if<Lit>(&v2));
            }
        }
    }

    void p_config_insert_long_clauses() {
        const auto& cdb = m_propagator.db();
        const auto& config = m_config_vars.back();
        for(CRef clause_ref = 1, n = cdb.literal_db_size(); 
            clause_ref < n; clause_ref = cdb.next_clause(clause_ref)) 
        {
            m_buffer.clear();
            bool got_true = false;
            for(SLit l : cdb.lits_of(clause_ref)) {
                LitOrVal v = p_config_lit_for(config, l);
                std::visit(overloaded{
                    [&] (bool b) {
                        if(b) {
                            got_true = true;
                        }
                    },
                    [&] (Lit l) {
                        m_buffer.push_back(l);
                    }
                }, v);
                if(got_true) break;
            }
            if(got_true) continue;
            if(m_buffer.empty()) {
                // empty clause, construction is infeasible.
                m_infeasible_by_construction = true;
                break;
            }
            if(m_buffer.size() > 1) {
                m_solver.add_clause(m_buffer.begin(), m_buffer.end());
            }
        }
    }

    void p_make_vertical_clauses() {
        HashSet<SLit> occurring;
        for(Vertex v : m_subproblem.uncovered_universe) {
            occurring.insert(v.first);
            occurring.insert(v.second);
        }
        const auto allowed = allowed_configurations();
        for(SLit l : occurring) {
            m_buffer.clear();
            bool got_true = false;
            for(std::size_t config_index = 0; config_index < allowed; ++config_index) {
                LitOrVal lv = p_config_lit_for(m_config_vars[config_index], l);
                std::visit(overloaded{
                    [&] (bool b) {
                        if(b) {
                            got_true = true;
                        }
                    },
                    [&] (Lit l) {
                        m_buffer.push_back(l);
                    }
                }, lv);
                if(got_true) break;
            }
            if(got_true) continue;
            if(m_buffer.empty()) {
                m_infeasible_by_construction = true;
                break;
            }
            m_solver.add_clause(m_buffer.begin(), m_buffer.end());
        }
    }

    void p_extract_solution() {
        const Var nall(m_propagator.db().num_vars());
        const auto& model_map = m_solver.get_model();
        if(m_current_solution.empty()) {
            m_current_solution.reserve(allowed_configurations());
            for(std::size_t i = 0; i < allowed_configurations(); ++i) {
                m_current_solution.emplace_back(nall, false);
            }
        }
        for(std::size_t i = 0; i < allowed_configurations(); ++i) {
            auto& current_config = m_config_vars[i];
            auto& current_solution_config = m_current_solution[i];
            for(Var j = 0; j != nall; ++j) {
                std::visit(overloaded{
                        [&] (bool b) {
                            current_solution_config[j] = b;
                        },
                        [&] (Lit l) {
                            current_solution_config[j] = model_map[l];
                        }
                    }, current_config[j]);
            }
        }
    }

    void p_identify_covered() {
        const auto allowed = allowed_configurations();
        m_covered_literal_partners = m_explicitized_literal_partners;
        for(SLit i = 0, n = m_literal_partners_of.size(); i < n; ++i) {
            for(std::size_t config_index = 0; config_index < allowed; ++config_index) {
                const auto& solution = m_current_solution[config_index];
                if(!lit::is_true_in(i, solution)) {
                    continue;
                }
                auto indices = m_literal_partners_of[i];
                for(std::size_t p = indices.first; p < indices.second; ++p) {
                    if(m_covered_literal_partners[p]) continue;
                    SLit partner = m_literal_partners_array[p];
                    if(lit::is_true_in(partner, solution)) {
                        m_covered_literal_partners[p].set();
                    }
                }
            }
        }
    }

    void p_make_all_explicit() {
        for(SLit i = 0, n = m_literal_partners_of.size(); i < n; ++i) {
            auto indices = m_literal_partners_of[i];
            for(std::size_t p = indices.first; p < indices.second; ++p) {
                if(m_explicitized_literal_partners[p]) continue;
                SLit partner = m_literal_partners_array[p];
                p_make_explicit_possibly_fixed(p, {i, partner});
            }
        }
    }

    void p_sample_more_explicit() {
        std::geometric_distribution<std::size_t> dist(0.025);
        auto& rng = sammy::rng();
        std::size_t skip = dist(rng);
        for(SLit i = 0, n = m_literal_partners_of.size(); i < n; ++i) {
            auto indices = m_literal_partners_of[i];
            for(std::size_t p = indices.first; p < indices.second; ++p) {
                if(skip > 0) {
                    --skip;
                    continue;
                }
                skip = dist(rng);
                if(m_explicitized_literal_partners[p]) continue;
                SLit partner = m_literal_partners_array[p];
                p_make_explicit_possibly_fixed(p, {i, partner});
            }
        }
    }

    void p_sample_uncovered() {
        std::size_t goal = (std::max)(m_num_explicitly_covered + 1, std::size_t(500));
        m_selected_uncovered.clear();
        for(SLit i = 0, n = m_literal_partners_of.size(); i < n; ++i) {
            std::size_t q = m_literal_partners_of[i].second;
            for(std::size_t p = m_literal_partners_of[i].first; p < q; ++p) {
                if(m_covered_literal_partners[p]) continue;
                m_selected_uncovered.push_back({{i, m_literal_partners_array[p]}, p});
            }
        }
        if(m_selected_uncovered.empty()) {
            return;
        }
        p_store_event("DETECTED_UNCOVERED_INTERACTIONS", {{"explicitly_covered", m_num_explicitly_covered},
                                                         {"num_uncovered", m_selected_uncovered.size()}}, 
                      "explicitly_covered", "num_uncovered");
        if(goal >= 0.75 * m_selected_uncovered.size()) return;
        std::shuffle(m_selected_uncovered.begin(), m_selected_uncovered.end(), sammy::rng());
        m_selected_uncovered.resize(goal);
    }

    std::optional<bool> p_solve_iteration() {
        if(m_abort.load()) return std::nullopt;
        auto res = m_solver.solve();
        if(m_abort.load()) return std::nullopt;
        return res;
    }

    /**
     * Called to explicitly require a currently
     * uncovered interaction to be covered.
     */
    void p_make_explicit_unfixed(std::size_t index, Vertex vertex) {
        const auto visit_combination = overloaded{
            [&] (bool b1, bool b2) {
                if(!b1 || !b2) return;
                throw std::logic_error("Asked to explicitize already fixed interaction!");
            },
            [&] (bool b1, Lit l2) {
                if(!b1) return;
                m_buffer.push_back(l2);
            },
            [&] (Lit l1, bool b2) {
                if(!b2) return;
                m_buffer.push_back(l1);
            },
            [&] (Lit l1, Lit l2) {
                Lit coverage_var = m_solver.new_var();
                m_solver.add_short_clause(coverage_var, -l1, -l2);
                m_solver.add_short_clause(-coverage_var, l1);
                m_solver.add_short_clause(-coverage_var, l2);
                m_buffer.push_back(coverage_var);
            }
        };

        ++m_num_explicitly_covered;
        m_explicitized_literal_partners[index].set();
        m_buffer.clear();
        for(const auto& config : m_config_vars) {
            LitOrVal lv1 = p_config_lit_for(config, vertex.first);
            LitOrVal lv2 = p_config_lit_for(config, vertex.second);
            std::visit(visit_combination, lv1, lv2);
        }
        if(m_buffer.empty()) {
            p_store_event("EMPTY_COVERAGE_CLAUSE");
            m_infeasible_by_construction = true;
        } else {
            m_solver.add_clause(m_buffer.begin(), m_buffer.end());
            m_buffer.clear();
        }
    }

    void p_make_explicit_possibly_fixed(std::size_t index, Vertex vertex) {
        for(const auto& config : m_config_vars) {
            LitOrVal lv1 = p_config_lit_for(config, vertex.first);
            LitOrVal lv2 = p_config_lit_for(config, vertex.second);
            if(std::holds_alternative<bool>(lv1) && std::holds_alternative<bool>(lv2)) {
                bool b1 = *std::get_if<bool>(&lv1);
                bool b2 = *std::get_if<bool>(&lv2);
                if(b1 && b2) {
                    ++m_num_explicitly_covered;
                    m_explicitized_literal_partners[index].set();
                    return;
                }
            }
        }
        p_make_explicit_unfixed(index, vertex);
    }

    template<typename... Args>
    void p_store_event(const char* name, OutputObject data, Args&&... keys) {
        if(!m_recorder) return;
        data["worker_id"] = m_worker_id;
        m_recorder->store_event(name, std::move(data), std::forward<Args>(keys)...);
    }

    void p_store_event(const char* name) {
        p_store_event(name, {{"worker_id", m_worker_id}});
    }
};

}

#endif
