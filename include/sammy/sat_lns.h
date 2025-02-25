#ifndef SAMMY_SAT_LNS_H_INCLUDED_
#define SAMMY_SAT_LNS_H_INCLUDED_

#include "clause_db.h"
#include "shared_db_propagator.h"
#include "barrage_lns_subproblem.h"
#include "barrage.h"
#include <variant>
#include <typeinfo>

namespace sammy {

/**
 * Create a non-incremental SAT solver
 * that attempts to improve the current
 * solution by at least one configuration.
 */
template<typename BasicSatSolver>
class FixedMESSATImprovementSolver {
  public:
    using SatSolver = BasicSatSolver;
    using Lit = typename SatSolver::Lit;
    using SLit = sammy::Lit;
    using LitOrVal = std::variant<Lit, bool>;

    static std::string name() {
        return std::string("SAT<") + BasicSatSolver::name() + ">";
    }

    /**
     * Create a new solver for the given subproblem.
     */
    FixedMESSATImprovementSolver(PortfolioSolver* portfolio, LNSSubproblem&& subproblem,
                                 SharedDBPropagator prop, EventRecorder* recorder, 
                                 std::size_t worker_id) :
        m_solver(),
        m_subproblem(std::move(subproblem)),
        m_recorder(recorder),
        m_worker_id(worker_id),
        m_propagator(std::move(prop)),
        m_num_configs_allowed(m_subproblem.removed_configurations.size() - 1),
        m_config_vars()
    {
        try {
            OutputObject event_data{
                {"num_uncovered", m_subproblem.uncovered_universe.size()},
                {"num_removed", m_subproblem.removed_configurations.size()},
                {"mes_size", m_subproblem.mutually_exclusive_set.size()},
                {"solver_name", m_solver.name()},
                {"num_configs_allowed", m_num_configs_allowed}
            };
            p_store_event("SAT_LNS_CONSTRUCTION_BEGIN", event_data,
                          "num_uncovered", "num_removed", "mes_size", "solver_name", "num_configs_allowed");
            if(m_num_configs_allowed < m_subproblem.mutually_exclusive_set.size()) {
                m_infeasible_by_construction = true;
                return;
            }
            // make variables for configurations and ensure their consistency
            p_make_config_vars();
            assert(m_config_vars.size() == m_num_configs_allowed);
            assert(std::all_of(m_config_vars.begin(), m_config_vars.end(), [&] (const auto& v) {
                return v.size() == m_propagator.db().num_vars();
            }));
            // make clauses for all individual literals in uncovered interactions
            p_make_vertical_clauses();
            // ensure that every uncovered universe element is covered
            p_make_coverage_clauses();
            p_store_event("SAT_LNS_CONSTRUCTION_DONE", std::move(event_data),
                          "num_uncovered", "num_removed", "mes_size", "solver_name", "num_configs_allowed");
        } catch(InterruptError&) {
            subproblem = std::move(m_subproblem);
            throw;
        }
    }

    /**
     * Get the number of allowed configurations.
     */
    std::size_t allowed_configurations() const noexcept {
        return m_subproblem.removed_configurations.size() - 1;
    }

    /**
     * Get the number of vertices in the mutually exclusive set.
     */
    std::size_t mes_size() const noexcept {
        return m_subproblem.mutually_exclusive_set.size();
    }

    /**
     * Get the number of interactions in the uncovered universe.
     */
    std::size_t num_uncovered() const noexcept {
        return m_subproblem.uncovered_universe.size();
    }

    /**
     * Get the vertices in the mutually exclusive set.
     */
    const std::vector<Vertex>& mes_vertices() const noexcept {
        return m_subproblem.mutually_exclusive_set;
    }

    /**
     * Move out the subproblem, including the
     * (potentially huge) list of interactions.
     */
    LNSSubproblem move_out_subproblem() noexcept {
        return std::move(m_subproblem);
    }

    /**
     * Get the solution of the solver, after solve()
     * returned true.
     */
    const std::vector<DynamicBitset>& get_solution() const {
        return m_current_solution;
    }

    /**
     * Abort the solve.
     */
    void abort() {
        m_solver.terminate();
    }

    /**
     * Solve the problem; return
     * nullopt if interrupted, false if
     * an improvement is impossible,
     * and true if an improvement was found.
     */
    std::optional<bool> solve() {
        if(m_infeasible_by_construction) {
            p_store_event("SAT_LNS_INFEASIBLE_BY_CONSTRUCTION");
            return false;
        }
        if(m_subproblem.uncovered_universe.empty()) {
            p_store_event("SAT_LNS_EMPTY_UNIVERSE");
            return true;
        }
        OutputObject event_data{
            {"num_uncovered", num_uncovered()},
            {"num_removed", m_subproblem.removed_configurations.size()},
            {"mes_size", mes_size()},
            {"solver_name", m_solver.name()}
        };
        if(get_and_clear_interrupt_flag()) {
            return std::nullopt;
        }
        p_store_event("SAT_LNS_BEGIN_SAT_SOLVE", event_data,
                      "num_uncovered", "num_removed", "mes_size", "solver_name");
        std::optional<bool> res = m_solver.solve();
        const char* event_name = !res ? "SAT_LNS_ABORT_SAT_SOLVE" : 
                                (*res ? "SAT_LNS_IMPROVEMENT_FOUND" :
                                 "SAT_LNS_SOLUTION_WAS_OPTIMAL");
        p_store_event(event_name, event_data,
                      "num_uncovered", "num_removed", "mes_size", "solver_name");
        if(!res) {
            return std::nullopt;
        }
        if(!*res) {
            return false;
        }
        p_extract_solution();
        return true;
    }

  private:
    // basic SAT-solver backend
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
    bool m_infeasible_by_construction{false};
    // how many configurations are needed?
    std::size_t m_num_configs_allowed;
    // (outer) index k: configuration index [0, allowed_configurations());
    // inner index: feature index [0, propagator().db().num_vars())
    std::vector<std::vector<LitOrVal>> m_config_vars;
    // buffer for solver literals/clauses
    std::vector<Lit> m_buffer;
    // buffer for the current solution
    std::vector<DynamicBitset> m_current_solution;

    // --------------------------------- IMPLEMENTATION ---------------------------------
    /**
     * Store an event in the recorder.
     */
    template<typename... Args>
    void p_store_event(const char* name, OutputObject data, Args&&... keys) {
        if(!m_recorder) return;
        data["worker_id"] = m_worker_id;
        m_recorder->store_event(name, std::move(data), std::forward<Args>(keys)...);
    }

    /**
     * Store an event in the recorder.
     */
    void p_store_event(const char* name) {
        p_store_event(name, {{"worker_id", m_worker_id}});
    }

    void p_make_config_vars() {
        for(std::size_t config_index = 0; config_index < m_num_configs_allowed; ++config_index) {
            m_propagator.reset_or_throw();
            if(config_index < m_subproblem.mutually_exclusive_set.size()) {
                if(push_vertex(m_propagator, m_subproblem.mutually_exclusive_set[config_index]) < 0) {
                    throw std::logic_error("Invalid interaction in MES!");
                }
            }
            p_make_config_from_propagator();
        }
        m_propagator.reset_or_throw();
        assert(m_num_configs_allowed == m_subproblem.removed_configurations.size() - 1);
        assert(m_num_configs_allowed == m_config_vars.size());
    }

    void p_make_config_from_propagator() {
        const auto nall = m_propagator.db().num_vars();
        std::vector<LitOrVal> config;
        config.reserve(nall);
        for(Var var = 0; var < nall; ++var) {
            SLit pos = lit::positive_lit(var);
            if(m_propagator.is_true(pos)) {
                config.emplace_back(std::in_place_type<bool>, true);
            } else if(m_propagator.is_false(pos)) {
                config.emplace_back(std::in_place_type<bool>, false);
            } else {
                config.emplace_back(std::in_place_type<Lit>, m_solver.new_var());
            }
        }
        assert(config.size() == nall);
        m_config_vars.emplace_back(std::move(config));
        assert(m_config_vars.back().size() == nall);
        p_config_insert_binaries();
        p_config_insert_long_clauses();
    }

    void p_config_insert_binaries() {
        const auto& db = m_propagator.db();
        const auto& config = m_config_vars.back();
        assert(db.num_vars() == config.size());
        for(auto clause : db.binary_clauses()) {
            LitOrVal lv1 = p_config_lit_for(config, clause.first);
            LitOrVal lv2 = p_config_lit_for(config, clause.second);
            std::visit(overloaded{
                [&] (bool b1, bool b2) {
                    // clause is either violated (indicating a logic error)
                    // or already satisfied by propagation
                    if(!b1 && !b2) {
                        throw std::logic_error("Infeasible binary clause!");
                    }
                },
                [&] (bool b1, Lit l2) {
                    // if b1 is true, clause is satisfied;
                    // if b1 is false, propagation should have made l2 true.
                    if(!b1) {
                        throw std::logic_error("Propagation should have made l2 true!");
                    }
                },
                [&] (Lit l1, bool b2) {
                    // if b2 is true, clause is satisfied;
                    // if b2 is false, propagation should have made l1 true.
                    if(!b2) {
                        throw std::logic_error("Propagation should have made l1 true!");
                    }
                },
                [&] (Lit l1, Lit l2) {
                    m_solver.add_short_clause(l1, l2);
                }
            }, lv1, lv2);
        }
    }

    void p_config_insert_long_clauses() {
        const auto& db = m_propagator.db();
        const auto& config = m_config_vars.back();
        for(CRef clause = 1, dbs = db.literal_db_size(); clause < dbs; clause = db.next_clause(clause)) {
            m_buffer.clear();
            bool got_true = false;
            for(SLit l : db.lits_of(clause)) {
                LitOrVal v = p_config_lit_for(config, l);
                if(std::holds_alternative<bool>(v)) {
                    if(*std::get_if<bool>(&v)) {
                        got_true = true;
                        break;
                    }
                    continue;
                }
                m_buffer.push_back(*std::get_if<Lit>(&v));
            }
            if(got_true) continue;
            if(m_buffer.empty()) {
                throw std::logic_error("Infeasible interaction in MES!");
            }
            if(m_buffer.size() == 1) {
                throw std::logic_error("Unit clause in MES - propagation should have assigned a value!");
            }
            m_solver.add_clause(m_buffer.begin(), m_buffer.end());
        }
    }

    void p_make_vertical_clauses() {
        if(m_infeasible_by_construction) return;
        HashSet<SLit> occurring;
        for(Vertex v : m_subproblem.uncovered_universe) {
            occurring.insert(v.first);
            occurring.insert(v.second);
        }
        for(SLit l : occurring) {
            m_buffer.clear();
            bool got_true = false;
            for(const auto& config : m_config_vars) {
                assert(config.size() == m_propagator.db().num_vars());
                LitOrVal lv = p_config_lit_for(config, l);
                if(std::holds_alternative<bool>(lv)) {
                    if(*std::get_if<bool>(&lv)) {
                        got_true = true;
                        break;
                    }
                    continue;
                }
                m_buffer.push_back(*std::get_if<Lit>(&lv));
            }
            if(got_true) continue;
            if(m_buffer.empty()) {
                m_infeasible_by_construction = true;
                break;
            }
            m_solver.add_clause(m_buffer.begin(), m_buffer.end());
        }
    }

    LitOrVal p_config_lit_for(const std::vector<LitOrVal>& config, SLit literal) const {
        auto var = lit::var(literal);
        bool is_pos = !lit::negative(literal);
        assert(var < config.size());
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

    /**
     * Extract the solution from the solver
     * and store it in m_current_solution.
     */
    void p_extract_solution() {
        const auto nclasses = m_config_vars.size();
        const auto nvars = m_propagator.db().num_vars();
        auto model_map = m_solver.get_model();
        std::vector<DynamicBitset> result(nclasses);
        for(std::size_t class_index = 0; class_index < nclasses; ++class_index) {
            DynamicBitset config(nvars, false);
            for(Var v = 0; v < nvars; ++v) {
                config[v] = p_get_class_value(model_map, class_index, v);
            }
            result[class_index] = std::move(config);
        }
        m_current_solution = std::move(result);
    }

    /**
     * Get the value of a literal or boolean
     * constant, using a model map.
     */
    template<typename ModelMap>
    bool p_get_value(const ModelMap& model, LitOrVal v) const noexcept {
        return std::visit(overloaded{
            [&] (bool b) { return b; },
            [&] (Lit l) { return model[l]; }
        }, v);
    }

    /**
     * Get the value of a variable in a class.
     */
    template<typename ModelMap>
    bool p_get_class_value(const ModelMap& model, std::size_t class_index, Var var) const noexcept {
        assert(class_index < m_config_vars.size());
        assert(var < m_config_vars[class_index].size());
        return p_get_value(model, m_config_vars[class_index][var]);
    }

    /**
     * Add variables and clauses to ensure that
     * every uncovered universe element is covered.
     */
    void p_make_coverage_clauses() {
        std::size_t count = 0;
        for(Vertex v : m_subproblem.uncovered_universe) {
            if(++count == 32768) {
                count = 0;
                throw_if_interrupted();
            }
            bool got_true = false;
            m_buffer.clear();
            // first scan partially-assigned configurations;
            // this avoids variable creation for vertices that
            // are already covered in a configuration
            for(const auto& config : m_config_vars) {
                LitOrVal l1 = p_config_lit_for(config, v.first);
                LitOrVal l2 = p_config_lit_for(config, v.second);
                if(std::holds_alternative<bool>(l1) || std::holds_alternative<bool>(l2)) {
                    bool l1true = false;
                    bool l2true = false;
                    if(std::holds_alternative<bool>(l1)) {
                        if(*std::get_if<bool>(&l1)) {
                            l1true = true;
                        } else {
                            // fixed false in this configuration
                            continue;
                        }
                    }
                    if(std::holds_alternative<bool>(l2)) {
                        if(*std::get_if<bool>(&l2)) {
                            l2true = true;
                        } else {
                            // fixed false in this configuration
                            continue;
                        }
                    }
                    if(l1true && l2true) {
                        got_true = true;
                        break;
                    }
                    if(l1true) {
                        m_buffer.push_back(*std::get_if<Lit>(&l2));
                    } else {
                        m_buffer.push_back(*std::get_if<Lit>(&l1));
                    }
                }
            }
            if(got_true) continue;
            // m_buffer contains the beginnings of a clause ensuring
            // that some configuration covers v; we now complete it
            // using extra variables for the other classes where
            // v is covered but needs two literals to be set properly
            for(const auto& config : m_config_vars) {
                LitOrVal l1v = p_config_lit_for(config, v.first);
                LitOrVal l2v = p_config_lit_for(config, v.second);
                if(std::holds_alternative<bool>(l1v) || std::holds_alternative<bool>(l2v)) {
                    continue;
                }
                Lit l1 = *std::get_if<Lit>(&l1v);
                Lit l2 = *std::get_if<Lit>(&l2v);
                Lit coverage_var = m_solver.new_var();
                m_solver.add_short_clause(coverage_var, -l1, -l2);
                m_solver.add_short_clause(-coverage_var, l1);
                m_solver.add_short_clause(-coverage_var, l2);
                m_buffer.push_back(coverage_var);
            }
            if(m_buffer.empty()) {
                m_infeasible_by_construction = true;
                break;
            }
            m_solver.add_clause(m_buffer.begin(), m_buffer.end());
        }
    }
};


/**
 * Non-incremental use of a SAT solver for improvement-by-one search.
 */
template<typename BaseSolver> class ImprovementSolver {
  public:
    using Lit = typename BaseSolver::Lit;
    using SLit = sammy::Lit;

    ImprovementSolver(SharedDBPropagator* propagator,
                      std::vector<Vertex> vertices,
                      std::vector<std::size_t> clique_indices,
                      std::size_t k) :
        m_base_solver(),
        m_propagator(propagator),
        m_all_vertices(std::move(vertices)),
        m_clique_indices(std::move(clique_indices)),
        m_class_vars(k),
        m_vertex_vars(m_all_vertices.size())
    {
        for(std::size_t i = 0; i < k; ++i) {
            if(i <= m_clique_indices.size()) m_propagator->reset_or_throw();
            if(i < m_clique_indices.size()) {
                p_make_class_with_clique_vertex(i, m_all_vertices[m_clique_indices[i]]);
            } else {
                p_make_unconstrained_class(i);
            }
            if(!p_clauses_for_class(i)) {
                m_construction_infeasible = true;
                break;
            }
        }
        if(!m_construction_infeasible) {
            p_clauses_for_lits();
        }
        if(!m_construction_infeasible) {
            for(std::size_t i = 0, n = m_all_vertices.size(); i < n; ++i) {
                if(!p_vertex_clause(i)) {
                    m_construction_infeasible = true;
                    break;
                }
            }
        }
    }

    std::optional<bool> solve(double time_limit = std::numeric_limits<double>::infinity()) {
        return m_base_solver.solve(time_limit);
    }

    std::vector<std::vector<bool>> get_solution() const {
        const auto nclasses = m_class_vars.size();
        const auto nvars = m_propagator->db().num_vars();
        auto model_map = m_base_solver.get_model();
        std::vector<std::vector<bool>> result(nclasses);
        for(std::size_t class_index = 0; class_index < nclasses; ++class_index) {
            std::vector<bool> config(nvars, false);
            for(Var v = 0; v < nvars; ++v) {
                config[v] = p_get_class_value(model_map, class_index, v);
            }
            result[class_index] = std::move(config);
        }
        return result;
    }

  private:
    template<typename ModelMap>
    bool p_get_value(const ModelMap& model, std::variant<Lit, bool> v) const {
        return std::visit(overloaded{
            [&] (bool& b) { return b; },
            [&] (Lit& l) { return model[l]; }
        }, v);
    }

    template<typename ModelMap>
    bool p_get_class_value(const ModelMap& model, std::size_t class_index, Var var) const {
        return p_get_value(model, m_class_vars[class_index][var]);
    }

    bool p_vertex_clause(std::size_t vertex_index) {
        auto& vvars = m_vertex_vars[vertex_index];
        m_clause_buffer.clear();
        for(auto& var : vvars) {
            if(!std::visit(
                overloaded{
                    [&] (bool& b) -> bool { return !b; },
                    [&] (Lit& l) -> bool {
                        m_clause_buffer.push_back(l);
                        return true;
                    }
                }, var)) 
            {
                return true;
            }
        }
        if(m_clause_buffer.empty()) return false;
        m_base_solver.add_clause(m_clause_buffer.begin(), m_clause_buffer.end());
        return true;
    }

    std::variant<Lit, bool> p_class_var_for(std::size_t class_index, SLit literal) {
        Var v = lit::var(literal);
        bool neg = lit::negative(literal);
        std::variant<Lit, bool> cv = m_class_vars[class_index][v];
        std::visit(overloaded{
            [&] (bool& b) { if(neg) b = !b; },
            [&] (Lit& l) { if(neg) l = -l; }
        }, cv);
        return cv;
    }

    void p_class_vars_from_prop(std::size_t i) {
        const auto& clauses = m_propagator->db();
        auto& cvars = m_class_vars[i];
        for(Var v = 0, nv = clauses.num_vars(); v < nv; ++v) {
            if(m_propagator->is_open(lit::positive_lit(v))) {
                cvars.emplace_back(std::in_place_index<0>, m_base_solver.new_var());
            } else {
                cvars.emplace_back(std::in_place_index<1>, m_propagator->is_true(lit::positive_lit(v)));
            }
        }
    }

    void p_vertex_vars_from_prop(std::size_t class_index) {
        for(std::size_t i = 0, n = m_all_vertices.size(); i < n; ++i) {
            Vertex v = m_all_vertices[i];
            auto v1 = p_class_var_for(class_index, v.first);
            auto v2 = p_class_var_for(class_index, v.second);
            auto& vvars = m_vertex_vars[i];
            auto one_fixed = [&] (Lit& l1, bool& b2) {
                if(!b2) {
                    vvars.emplace_back(std::in_place_index<1>, false);
                } else {
                    vvars.emplace_back(std::in_place_index<0>, l1);
                }
            };
            std::visit(overloaded{
                [&] (bool& b1, bool& b2) {
                    vvars.emplace_back(std::in_place_index<1>, b1 && b2);
                },
                [&] (Lit& l1, Lit& l2) {
                    Lit vvar = m_base_solver.new_var();
                    vvars.emplace_back(std::in_place_index<0>, vvar);
                    m_base_solver.add_short_clause(-vvar, l1);
                    m_base_solver.add_short_clause(-vvar, l2);
                    m_base_solver.add_short_clause(-l1, -l2, vvar);
                },
                [&] (Lit& l1, bool& b2) { one_fixed(l1, b2); },
                [&] (bool& b1, Lit& l2) { one_fixed(l2, b1); }
            }, v1, v2);
        }
    }

    void p_make_unconstrained_class(std::size_t i) {
        p_class_vars_from_prop(i);
        p_vertex_vars_from_prop(i);
    }

    void p_make_class_with_clique_vertex(std::size_t class_index, Vertex clique_vertex) {
        if(m_propagator->is_open(clique_vertex.first))
            m_propagator->push_level(clique_vertex.first);
        if(m_propagator->is_open(clique_vertex.second))
            m_propagator->push_level(clique_vertex.second);
        p_class_vars_from_prop(class_index);
        p_vertex_vars_from_prop(class_index);
    }

    bool p_clauses_for_class(std::size_t i) {
        const auto& db = m_propagator->db();
        for(auto [l1,l2] : db.binary_clauses()) {
            SLit ls[2] = {l1, l2};
            if(!p_to_clause(i, +ls, ls+2)) return false;
        }
        for(CRef c = 1, n = db.literal_db_size(); c < n; c = db.next_clause(c)) {
            auto lits = db.lits_of(c);
            if(!p_to_clause(i, lits.begin(), lits.end())) return false;
        }
        return true;
    }

    template<typename Iterator>
    bool p_to_clause(std::size_t class_index, Iterator slit_begin, Iterator slit_end)
    {
        bool found_true = false;
        m_lit_buffer.clear();
        for(SLit sl : IteratorRange{slit_begin, slit_end}) {
            if(m_propagator->is_open(sl)) {
                m_lit_buffer.push_back(sl);
            } else if(m_propagator->is_true(sl)) {
                found_true = true; break;
            }
        }
        if(found_true) return true;
        if(m_lit_buffer.empty()) return false;
        m_clause_buffer.clear();
        std::transform(m_lit_buffer.begin(), m_lit_buffer.end(),
                       std::back_inserter(m_clause_buffer),
                       [&] (SLit sl) { return std::get<Lit>(p_class_var_for(class_index, sl)); });
        m_base_solver.add_clause(m_clause_buffer.begin(), m_clause_buffer.end());
        return true;
    }

    /**
     * For each literal that occurs in any
     * vertex that we have to cover: add a clause
     * that ensures that at least one configuration
     * has that literal set to true, unless that
     * is already trivially satisfied.
     */
    void p_clauses_for_lits() {
        {
            HashSet<SLit> occurring;
            for(Vertex v : m_all_vertices) {
                occurring.insert(v.first);
                occurring.insert(v.second);
            }
            m_lit_buffer.reserve(occurring.size());
            m_lit_buffer.assign(occurring.begin(), occurring.end());
        }
        std::sort(m_lit_buffer.begin(), m_lit_buffer.end());
        for(SLit l : m_lit_buffer) {
            m_clause_buffer.clear();
            bool got_true = false;
            for(std::size_t config_index = 0, num_configs = m_class_vars.size(); 
                config_index < num_configs; ++config_index) 
            {
                auto lv = p_class_var_for(config_index, l);
                std::visit(overloaded{
                    [&] (bool b) {
                        if(b) {
                            got_true = true;
                        }
                    },
                    [&] (Lit l) {
                        m_clause_buffer.push_back(l);
                    }
                }, lv);
                if(got_true) break;
            }
            if(got_true) continue;
            if(m_clause_buffer.empty()) {
                m_construction_infeasible = true;
                break;
            }
            m_base_solver.add_clause(m_clause_buffer.begin(), m_clause_buffer.end());
        }
        m_lit_buffer.clear();
    }
    
    BaseSolver m_base_solver;
    SharedDBPropagator* m_propagator;
    std::vector<Vertex> m_all_vertices;
    std::vector<std::size_t> m_clique_indices;
    // m_class_vars[c][v] -> variable (or fixed value) for class c and basic variable v
    std::vector<std::vector<std::variant<Lit,bool>>> m_class_vars;
    // m_vertex_vars[vi][c] -> variable (or fixed value) for class c and vertex with index vi
    std::vector<std::vector<std::variant<Lit,bool>>> m_vertex_vars;
    bool m_construction_infeasible = false;
    std::vector<SLit> m_lit_buffer;
    std::vector<Lit> m_clause_buffer;
};

}

#endif
