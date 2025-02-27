#ifndef SAMMY_ORTOOLS_SOLVER_H_INCLUDED_
#define SAMMY_ORTOOLS_SOLVER_H_INCLUDED_
#ifdef SAMMY_ORTOOLS_SUPPORTED
#if SAMMY_ORTOOLS_SUPPORTED

#include <ortools/sat/cp_model.h>
#include <ortools/sat/cp_model.pb.h>
#include <ortools/sat/cp_model_solver.h>
#include <ortools/sat/model.h>
#include <ortools/sat/sat_parameters.pb.h>
#include <ortools/util/time_limit.h>

#include "clause_db.h"
#include "dynamic_bitset.h"
#include "output.h"
#include "pair_infeasibility_map.h"
#include "partial_solution.h"
#include "vertex_operations.h"

#include <atomic>

namespace sammy {

class ORToolsSolver {
  public:
    /**
     * Create an ORTools solver.
     */
    ORToolsSolver(std::vector<Vertex> considered_vertices,
                  PairInfeasibilityMap* infeasibility_map, ClauseDB& clause_db,
                  std::vector<Vertex> best_local_mes,
                  std::size_t /*best_global_mes*/, std::size_t best_global_lb,
                  std::vector<DynamicBitset> covering_assignments,
                  std::optional<DynamicBitset> remaining_configuration)
        : m_all_vertices(std::move(considered_vertices)),
          m_infeasibility_map(infeasibility_map), m_clauses(&clause_db),
          m_propagator(m_clauses), m_best_local_mes(std::move(best_local_mes)),
          m_covering_assignments(std::move(covering_assignments)),
          m_remaining_configuration(std::move(remaining_configuration)),
          m_lower_bound(m_best_local_mes.size()),
          m_initial_global_lower_bound(best_global_lb),
          m_partial_solution(m_clauses->num_vars(), m_infeasibility_map) {}

    /**
     * Turn on 'single improvement' mode, where we do not
     * minimize a objective function but simply ask for
     * an improvement by at least one configuration.
     */
    void set_find_single_improvement(bool find_single_improvement) {
        if (m_initialized_objective)
            throw std::logic_error(
                "Called set_find_single_improvement too late!");
        m_want_single_improvement = find_single_improvement;
    }

    /**
     * Abort the solve (on timeout, optimality, LNS with better solution, ...).
     */
    void abort() noexcept { m_abort_flag.store(true); }

    /**
     * The size of the graph (all vertices) we are working on.
     */
    std::size_t get_graph_size() const noexcept {
        return m_all_vertices.size();
    }

    /**
     * Enum describing the outcome of the solve process.
     */
    enum class SolveResult {
        IMPROVED_SOLUTION, //< We found a new, optimal solution.
        ABORTED,           //< We were aborted (timeout, ...) before completion.
        SOLUTION_WAS_OPTIMAL //< We completed but the initial solution was
                             //optimal.
    };

    /**
     * Run the solution process.
     * @return SolveResult describing the outcome.
     */
    SolveResult solve() {
        p_construct_model();
        if (m_abort_flag.load())
            return SolveResult::ABORTED;
        Model model = m_model.Build();
        if (m_abort_flag.load())
            return SolveResult::ABORTED;
        p_report_event("CPSAT_MODEL_BUILT", {});
        operations_research::sat::Model model_config;
        model_config.GetOrCreate<operations_research::TimeLimit>()
            ->RegisterExternalBooleanAsLimit(&m_abort_flag);
        p_report_event("CPSAT_SOLVE_BEGIN", {});
        operations_research::sat::CpSolverResponse response;
        response = operations_research::sat::SolveCpModel(model, &model_config);
        switch (response.status()) {
        case operations_research::sat::CpSolverStatus::OPTIMAL:
            p_report_event("CPSAT_SOLVE_OPTIMAL", {});
            break;

        case operations_research::sat::CpSolverStatus::FEASIBLE:
            p_report_event("CPSAT_SOLVE_FEASIBLE", {});
            break;

        default:
        case operations_research::sat::CpSolverStatus::UNKNOWN:
            p_report_event("CPSAT_SOLVER_ABORTED", {});
            return SolveResult::ABORTED;

        case operations_research::sat::CpSolverStatus::MODEL_INVALID:
            throw std::logic_error("Model invalid error!");

        case operations_research::sat::CpSolverStatus::INFEASIBLE:
            p_report_event("CPSAT_SOLVE_INFEASIBLE", {});
            m_lower_bound = m_covering_assignments.size();
            if (m_lower_bound > m_initial_global_lower_bound &&
                m_lower_bound_callback)
            {
                m_lower_bound_callback(m_lower_bound, m_all_vertices);
            }
            return SolveResult::SOLUTION_WAS_OPTIMAL;
        }
        if (!m_want_single_improvement) {
            std::size_t new_bound =
                std::size_t(response.best_objective_bound() + 0.01);
            if (new_bound > m_lower_bound) {
                m_lower_bound = new_bound;
                if (m_lower_bound > m_initial_global_lower_bound &&
                    m_lower_bound_callback)
                {
                    m_lower_bound_callback(m_lower_bound, m_all_vertices);
                }
            }
        }
        p_extract_partial_solution(response);
        return SolveResult::IMPROVED_SOLUTION;
    }

    /**
     * Get the best solution found. Only valid if
     * the SolveResult was IMPROVED_SOLUTION.
     *
     * @return The best solution found.
     */
    PartialSolution get_partial_solution() const { return m_partial_solution; }

    /**
     * @return A reference to the covering assignments given initially.
     */
    const std::vector<DynamicBitset>&
    get_covering_assignments() const noexcept {
        return m_covering_assignments;
    }

    /**
     * Check whether we may have found an improvement
     * on the global lower bound.
     */
    bool improved_global_bound() const noexcept {
        return m_lower_bound > m_initial_global_lower_bound;
    }

    /**
     * Check whether we may have found an improvement
     * on the global MES size. This method exists for
     * compatibility with the other solvers and always returns false.
     */
    bool improved_mes() const noexcept { return false; }

    /**
     * Get a copy of the best clique found.
     */
    std::vector<Vertex> get_best_mes() const { return m_best_local_mes; }

    /**
     * Get the best lower bound on the number of configurations
     * needed for the subgraph.
     */
    std::size_t get_best_bound() const noexcept { return m_lower_bound; }

    /**
     * Set the event recorder to report events to;
     * this is purely optional.
     */
    void set_event_recorder(EventRecorder* recorder) noexcept {
        m_local_recorder = recorder;
    }

    /**
     * Get the vertices inducing the subgraph on which
     * we established the lower bound.
     */
    std::vector<Vertex> get_best_bound_subgraph() const {
        return m_all_vertices;
    }

    /**
     * Set a callback to be invoked each time we try to search for a new clique.
     * It can return a clique (containing any vertices); the solver will filter
     * out all irrelevant vertices and only consider the ones in m_all_vertices.
     * This method exists for compatibility with other solvers and does nothing.
     */
    void set_clique_candidate_callback(std::function<std::vector<Vertex>()>) {}

    /**
     * Set a callback to be invoked each time this solver finds a new lower
     * bound. The callback is involved with the lower bound and the vertices of
     * the subgraph inducing it.
     */
    void set_lower_bound_callback(
        std::function<void(std::size_t, const std::vector<Vertex>&)> callback) {
        m_lower_bound_callback = std::move(callback);
    }

  private:
    /**
     * All vertices in the subproblem we are considering.
     */
    std::vector<Vertex> m_all_vertices;

    /**
     * The pair infeasibility map.
     */
    PairInfeasibilityMap* m_infeasibility_map;

    /**
     * The clause DB.
     */
    ClauseDB* m_clauses;

    /**
     * A propagator connected to m_clauses.
     */
    SharedDBPropagator m_propagator;

    /**
     * Symmetry-breaking MES.
     */
    std::vector<Vertex> m_best_local_mes;

    /**
     * The removed configurations constituting this subproblem.
     */
    std::vector<DynamicBitset> m_covering_assignments;

    /**
     * A single remaining configuration.
     * If it is available, can help breaking symmetries.
     */
    std::optional<DynamicBitset> m_remaining_configuration;

    /**
     * The best lower bound we have found.
     */
    std::size_t m_lower_bound;

    /**
     * The best lower bound we have found for the complete problem.
     */
    std::size_t m_initial_global_lower_bound;

    /**
     * Pointer to the event recorder to store our events in.
     */
    EventRecorder* m_local_recorder{nullptr};

    /**
     * A flag to control whether we do minimization
     * or just find any improvement.
     */
    bool m_want_single_improvement{false};

    /**
     * A flag to store whether we have initialized
     * the objective function.
     */
    bool m_initialized_objective{false};

    /**
     * Callback invoked with new lower bounds
     * and the subgraphs on which they were obtained.
     */
    std::function<void(std::size_t, const std::vector<Vertex>&)>
        m_lower_bound_callback;

    /**
     * A flag to abort the CP-SAT solve.
     */
    std::atomic<bool> m_abort_flag{false};

    /**
     * Key types from CP-SAT.
     */
    using ModelBuilder = operations_research::sat::CpModelBuilder;
    using Model = operations_research::sat::CpModelProto;
    using SolverResponse = operations_research::sat::CpSolverResponse;
    using BoolVar = operations_research::sat::BoolVar;
    using VarOrVal = std::variant<BoolVar, bool>;

    /**
     * The CP-SAT model (or rather, its builder).
     */
    ModelBuilder m_model;

    /**
     * The CP-SAT variables (or fixed values) for the configurations.
     * Index [c][x] has the value of variable x in configuration c.
     */
    std::vector<std::vector<VarOrVal>> m_configuration_variables;

    /**
     * The interaction coverage variables (or fixed values).
     * Index [c][i] checks for coverage of the ith interaction in configuration
     * c.
     */
    std::vector<std::vector<VarOrVal>> m_interaction_coverage_variables;

    /**
     * Variables for the usage/activity of configurations.
     */
    std::vector<VarOrVal> m_configuration_usage_variables;

    /**
     * Buffer for clause building.
     */
    std::vector<BoolVar> m_clause_buffer;

    /**
     * Partial solution (initially empty).
     */
    PartialSolution m_partial_solution;

    /**
     * Construct the CP-SAT model.
     */
    void p_construct_model() {
        std::size_t num_copies_needed = m_covering_assignments.size() - 1;
        if (!m_want_single_improvement) {
            // turn on single improvement if we are in a -1 situation
            if (num_copies_needed == m_lower_bound) {
                m_want_single_improvement = true;
            }
        }
        m_initialized_objective = true;
        p_report_event(
            "BEGIN_CONSTRUCT_CPSAT_MODEL",
            {{"max_num_configurations", num_copies_needed},
             {"best_lb", m_lower_bound},
             {"best_global_lb", m_initial_global_lower_bound},
             {"best_local_mes_size", m_best_local_mes.size()},
             {"num_interactions", m_all_vertices.size()},
             {"single_improvement_search", m_want_single_improvement}},
            "max_num_configurations", "num_interactions",
            "single_improvement_search", "best_lb", "best_global_lb",
            "best_local_mes_size");
        for (std::size_t copy_index : range(num_copies_needed)) {
            p_add_copy(copy_index);
        }
        p_add_coverage_constraints();
        p_add_objective();
        p_report_event("CPSAT_MODEL_BUILDER_POPULATED", {});
    }

    /**
     * Extract a partial solution from the CP-SAT solution.
     */
    void p_extract_partial_solution(const SolverResponse& response) {
        std::vector<DynamicBitset> new_configurations;
        std::size_t num_vars = m_clauses->num_vars();
        std::size_t num_copies = m_covering_assignments.size() - 1;
        for (std::size_t copy_index : range(num_copies)) {
            VarOrVal usage_var = m_configuration_usage_variables[copy_index];
            if (!p_get_value(usage_var, response))
                break;
            DynamicBitset next_configuration(num_vars, false);
            for (std::size_t v : range(num_vars)) {
                VarOrVal entry = m_configuration_variables[copy_index][v];
                if (p_get_value(entry, response)) {
                    next_configuration[v].set();
                }
            }
            new_configurations.emplace_back(std::move(next_configuration));
        }
        m_partial_solution = PartialSolution(Var(num_vars), m_infeasibility_map,
                                             new_configurations.begin(),
                                             new_configurations.end());
    }

    /**
     * Construct another configuration copy.
     */
    void p_add_copy(std::size_t copy_index) {
        p_add_copy_vars(copy_index);
        p_copy_add_clauses(copy_index);
        p_copy_add_usage(copy_index);
        p_copy_add_coverage(copy_index);
    }

    /**
     * Add configuration copy usage variables.
     */
    void p_copy_add_usage(std::size_t copy_index) {
        if (copy_index < m_best_local_mes.size() || m_want_single_improvement) {
            // definitely used
            m_configuration_usage_variables.emplace_back(
                std::in_place_type_t<bool>(), true);
            return;
        }
        m_configuration_usage_variables.push_back(m_model.NewBoolVar());
    }

    /**
     * Add variables for a configuration copy.
     */
    void p_add_copy_vars(std::size_t copy_index) {
        std::vector<VarOrVal> copy_vars;
        bool all_open = true;
        if (copy_index < m_best_local_mes.size()) {
            all_open = false;
            m_propagator.reset_to_zero();
            Vertex vmes = m_best_local_mes[copy_index];
            if (push_vertex(m_propagator, vmes) < 0) {
                throw std::logic_error("Invalid interaction in MES!");
            }
        }
        for (Var v = 0, nv = m_clauses->num_vars(); v < nv; ++v) {
            Lit p = lit::positive_lit(v);
            if (all_open || m_propagator.is_open(p)) {
                copy_vars.emplace_back(std::in_place_type_t<BoolVar>{},
                                       m_model.NewBoolVar());
            } else {
                bool val = m_propagator.is_true(p);
                copy_vars.emplace_back(std::in_place_type_t<bool>{}, val);
            }
        }
        m_configuration_variables.emplace_back(std::move(copy_vars));
    }

    VarOrVal p_value_of(std::size_t copy_index, Lit l) {
        Var v = lit::var(l);
        bool negated = lit::negative(l);
        return std::visit(
            overloaded{
                [&](BoolVar& b) -> VarOrVal { return negated ? b.Not() : b; },
                [&](bool& b) -> VarOrVal { return negated ? !b : b; }},
            m_configuration_variables[copy_index][v]);
    }

    /**
     * Add interaction coverage variables for the current copy.
     */
    void p_copy_add_coverage(std::size_t copy_index) {
        bool all_open = (copy_index >= m_best_local_mes.size());
        auto usage_var = m_configuration_usage_variables[copy_index];
        bool definitely_used = std::holds_alternative<bool>(usage_var) &&
                               std::get<bool>(usage_var);
        m_interaction_coverage_variables.emplace_back();
        auto& coverage_vars = m_interaction_coverage_variables[copy_index];
        for (std::size_t vi = 0, n = m_all_vertices.size(); vi < n; ++vi) {
            Vertex v = m_all_vertices[vi];
            if (all_open || (m_propagator.is_open(v.first) &&
                             m_propagator.is_open(v.second)))
            {
                BoolVar vcov = m_model.NewBoolVar();
                coverage_vars.push_back(vcov);
                BoolVar v1 = std::get<BoolVar>(p_value_of(copy_index, v.first));
                BoolVar v2 =
                    std::get<BoolVar>(p_value_of(copy_index, v.second));
                m_model.AddBoolOr({vcov.Not(), v1});
                m_model.AddBoolOr({vcov.Not(), v2});
                if (definitely_used || m_remaining_configuration) {
                    // if the class is definitely used, or if we have
                    // a remaining configuration covering none of
                    // m_all_vertices, we can just add the backwards clause
                    // ensuring vcov == v1 & v2
                    m_model.AddBoolOr({vcov, v1.Not(), v2.Not()});
                } else {
                    // otherwise, we ensure that vcov == u & v1 & v2
                    // where u is the variable indicating the class is used
                    BoolVar u = std::get<BoolVar>(usage_var);
                    m_model.AddBoolOr({vcov.Not(), u});
                    m_model.AddBoolOr({vcov, u.Not(), v1.Not(), v2.Not()});
                }
                continue;
            }
            // here, we are in a definitely used class!
            // otherwise, all_open would be true.
            if (m_propagator.is_false(v.first) ||
                m_propagator.is_false(v.second))
            {
                coverage_vars.push_back(false); // v is not covered here!
                continue;
            }
            if (m_propagator.is_true(v.first) && m_propagator.is_true(v.second))
            {
                coverage_vars.push_back(true); // v is covered here!
                continue;
            }
            // one true, one open: if the open variable is made true, v is
            // covered here!
            Lit vopen = m_propagator.is_open(v.first) ? v.first : v.second;
            coverage_vars.push_back(p_value_of(copy_index, vopen));
        }
    }

    /**
     * Add the relevant/remaining clauses for the given copy to the model.
     */
    void p_copy_add_clauses(std::size_t copy_index) {
        for (Lit u : m_clauses->unary_literals()) {
            const Lit ubuf[1] = {u};
            p_copy_add_clause(copy_index, ubuf, ubuf + 1);
        }
        for (auto b : m_clauses->binary_clauses()) {
            const Lit bbuf[2] = {b.first, b.second};
            p_copy_add_clause(copy_index, bbuf, bbuf + 2);
        }
        for (CRef c = 1, n = m_clauses->literal_db_size(); c < n;
             c = m_clauses->next_clause(c))
        {
            auto lits = m_clauses->lits_of(c);
            p_copy_add_clause(copy_index, lits.begin(), lits.end());
        }
    }

    /**
     * Add a clause attached to the given copy.
     * Only add literals that aren't set to false.
     * Only add clauses that aren't already satisfied.
     */
    template <typename LitIter>
    void p_copy_add_clause(std::size_t copy_index, LitIter lbegin,
                           LitIter lend) {
        bool is_already_satisfied = false;
        m_clause_buffer.clear();
        const auto& copy_vars = m_configuration_variables[copy_index];
        for (Lit l : IteratorRange{lbegin, lend}) {
            VarOrVal entry = copy_vars[lit::var(l)];
            if (std::holds_alternative<bool>(entry)) {
                if (std::get<bool>(entry) != lit::negative(l)) {
                    is_already_satisfied = true;
                    break;
                }
            } else {
                BoolVar bv = std::get<BoolVar>(entry);
                m_clause_buffer.push_back(lit::negative(l) ? bv.Not() : bv);
            }
        }
        if (!is_already_satisfied) {
            if (m_clause_buffer.empty()) {
                throw std::logic_error(
                    "Invalid MES vertex or unsatisfiable configuration model!");
            }
            m_model.AddBoolOr(m_clause_buffer);
        }
    }

    /**
     * Add constraints to make sure all interactions are covered.
     */
    void p_add_coverage_constraints() {
        std::size_t num_copies = m_covering_assignments.size() - 1;
        for (std::size_t vi : range(m_all_vertices.size())) {
            m_clause_buffer.clear();
            bool already_satisfied = false;
            for (std::size_t ci : range(num_copies)) {
                VarOrVal vcov = m_interaction_coverage_variables[ci][vi];
                if (std::holds_alternative<bool>(vcov)) {
                    if (std::get<bool>(vcov)) {
                        already_satisfied = true;
                        break;
                    }
                } else {
                    m_clause_buffer.push_back(std::get<BoolVar>(vcov));
                }
            }
            if (!already_satisfied) {
                if (m_clause_buffer.empty()) {
                    throw std::logic_error("Uncoverable interaction in model!");
                }
                m_model.AddBoolOr(m_clause_buffer);
            }
        }
    }

    /**
     * Add the objective function (if we are not doing single improvements).
     */
    void p_add_objective() {
        if (m_want_single_improvement)
            return;
        std::size_t num_copies = m_covering_assignments.size() - 1;
        std::optional<BoolVar> previous = std::nullopt;
        operations_research::sat::LinearExpr sum;
        for (std::size_t i = 0; i < num_copies; ++i) {
            VarOrVal usage = m_configuration_usage_variables[i];
            if (std::holds_alternative<bool>(usage)) {
                if (!std::get<bool>(usage)) {
                    throw std::logic_error("Configuration fixed to unused!");
                }
                sum += 1;
                continue;
            } else {
                BoolVar vu = std::get<BoolVar>(usage);
                if (previous) {
                    m_model.AddBoolOr({vu.Not(), *previous});
                }
                sum += vu;
                previous = vu;
                if (m_remaining_configuration) {
                    // if we have a remaining configuration, we can
                    // force unused classes to the remaining configuration
                    // (which does not cover anything in m_all_vertices)
                    p_compute_fixed_from_remaining(i);
                    m_model.AddBoolAnd(m_clause_buffer).OnlyEnforceIf(vu.Not());
                }
            }
        }
        m_model.Minimize(sum);
    }

    /**
     * Compute fixed setting for unused configurations
     * from the remaining configuration (whose covered interactions are
     * all not in m_all_vertices).
     */
    void p_compute_fixed_from_remaining(std::size_t copy_index) {
        m_clause_buffer.clear();
        const auto& copy_vars = m_configuration_variables[copy_index];
        const auto& fixed_conf = *m_remaining_configuration;
        for (Var v : range(Var(m_clauses->num_vars()))) {
            BoolVar xv = std::get<BoolVar>(copy_vars[v]);
            if (fixed_conf[v]) {
                m_clause_buffer.push_back(xv);
            } else {
                m_clause_buffer.push_back(xv.Not());
            }
        }
    }

    /**
     * Report events if we have a recorder.
     */
    template <typename... PrintArgs>
    void p_report_event(std::string event_name, OutputObject data,
                        PrintArgs&&... args) {
        if (m_local_recorder) {
            m_local_recorder->store_event(std::move(event_name),
                                          std::move(data),
                                          std::forward<PrintArgs>(args)...);
        }
    }

    /**
     * Get a boolean value in the given solution.
     */
    bool p_get_value(VarOrVal v, const SolverResponse& response) {
        return std::visit(
            overloaded{
                [&](BoolVar& b) -> bool {
                    return operations_research::sat::SolutionBooleanValue(
                        response, b);
                },
                [&](bool& b) -> bool { return b; }},
            v);
    }
};

} // namespace sammy

#endif
#endif // SAMMY_ORTOOLS_SUPPORTED
#endif
