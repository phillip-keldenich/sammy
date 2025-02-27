#ifndef SAMMY_GUROBI_CLIQUE_SOLVER_G2_H_INCLUDED_
#define SAMMY_GUROBI_CLIQUE_SOLVER_G2_H_INCLUDED_

#include "fast_clique.h"
#include "gurobi.h"
#include "initial_coloring_heuristic.h"
#include "literals.h"
#include "lower_bound_mip_settings.h"
#include "output.h"
#include "rng.h"
#include "universe_subgraph.h"
#include "vertex_operations.h"

namespace sammy {

/**
 * Gurobi cut-and-price clique solver working on an extensible
 * subgraph of the conflicts in the interaction universe.
 */
class GurobiCliqueSolverG2 {
  public:
    template <typename BitsetType>
    GurobiCliqueSolverG2(UniverseSubgraph* subgraph, Var n_all,
                         EventRecorder* recorder,
                         const std::vector<BitsetType>& initial_sample,
                         const std::vector<Vertex>& initial_mut_ex_set,
                         const LowerBoundMIPConfig& config)
        : m_subgraph(subgraph), m_recorder(recorder), m_config(config),
          m_model(sammy::gurobi_environment(config.quiet_gurobi)),
          m_nall(n_all), m_indset_builder(subgraph->independent_set_builder()),
          m_clique_builder(subgraph->clique_builder()),
          m_subgraph_upper_bound(initial_sample.size()),
          m_best_sample(p_to_sample(initial_sample)),
          m_best_mes(initial_mut_ex_set) {
        m_model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
        m_model.set(GRB_IntParam_Threads,
                    int(m_subgraph->thread_group().num_threads() + 1));
        m_model.set(GRB_IntParam_Method, GRB_METHOD_DUAL);
        p_create_initial_vars();
        add_configurations(initial_sample.begin(), initial_sample.end());
    }

    /**
     * Make this solver abortable.
     * This means that the 'abort' method can
     * be used to attempt to abort long-running methods
     * before they are completed.
     * Such methods will behave as if a timeout occurred.
     */
    void make_abortable() {
        class GurobiAbortCallback : public GRBCallback {
          public:
            GurobiAbortCallback(GurobiCliqueSolverG2* that) : that(that) {}

            void callback() override {
                if (that->m_aborted.load()) {
                    this->abort();
                }
            }

          private:
            GurobiCliqueSolverG2* that;
        };
        if (!m_abortable) {
            m_abortable = true;
            m_callback = std::make_unique<GurobiAbortCallback>(this);
            m_model.setCallback(m_callback.get());
        }
    }

    /**
     * Attempt to stop long-running methods at the next opportunity.
     * Such methods will return a timeout result.
     */
    void abort() {
        if (!m_abortable)
            throw std::logic_error("Trying to abort non-abortable solver!");
        m_aborted.store(true);
    }

    /**
     * Get the thread group for this solver.
     */
    ThreadGroup<void>& thread_group() const noexcept {
        return m_subgraph->thread_group();
    }

    std::size_t num_lp_variables() const {
        return std::size_t(m_model.get(GRB_IntAttr_NumVars));
    }

    std::size_t num_lp_constraints() const {
        return std::size_t(m_model.get(GRB_IntAttr_NumConstrs));
    }

    std::size_t num_lp_nonzeros() const {
        return std::size_t(m_model.get(GRB_IntAttr_NumNZs));
    }

    std::size_t num_lp_solve_calls() const { return m_num_solve_calls; }

    std::size_t num_graph_vertices() const noexcept { return m_subgraph->n(); }

    /**
     * Add configurations as constraints.
     */
    template <typename ConfigurationIterator>
    inline void add_configurations(ConfigurationIterator cbeg,
                                   ConfigurationIterator cend);

    /**
     * Add a configuration as constraint.
     */
    template <typename ConfigurationType>
    inline void add_configuration(const ConfigurationType& cfg);

    /**
     * Add a set of vertices to the model and the subgraph.
     * Does not add vertices to independent sets, but adds
     * vertex variables to all independent set constraints
     * for independent sets the vertices are already added to.
     */
    inline void add_vertices(const std::vector<Vertex>& vertices);

    /**
     * Add the vertices from vertices using add_vertices that are not already in
     * the graph.
     */
    inline void add_new_vertices(const std::vector<Vertex>& vertices) {
        std::vector<Vertex> actual_vertices;
        std::copy_if(vertices.begin(), vertices.end(),
                     std::back_inserter(actual_vertices),
                     [&](Vertex v) { return !m_subgraph->has_vertex(v); });
        if (!actual_vertices.empty())
            add_vertices(actual_vertices);
    }

    /**
     * Solve the linear relaxation of the current model
     * and try to extract (essentially by greedy rounding) a
     * new clique from the solution. The solution may
     * still have violated non-edges, but any extracted clique is valid.
     */
    inline SolverState solve_relaxation(
        double time_allowed = std::numeric_limits<double>::infinity());

    /**
     * Solve the full relaxation, i.e., the linear relaxation including
     * all non-edge constraints on the current subgraph.
     */
    inline SolverState solve_full_relaxation(
        double time_allowed = std::numeric_limits<double>::infinity());

    /**
     * Solve the problem as MIP.
     * Temporarily transforms the problem into a MIP.
     */
    inline SolverState
    solve_as_mip(double time_allowed = std::numeric_limits<double>::infinity());

    /**
     * Find and prohibit violated nonedge constraints.
     * Returns true if a violated nonedge was found.
     */
    inline bool prohibit_violated_nonedges();

    /**
     * Add an independent set encoded by a propagator.
     * If the propagator encodes a complete configuration,
     * this adds a complete independent set; otherwise, it
     * adds an incomplete set that may be extended later.
     */
    inline void add_potentially_incomplete_independent_set(
        const SharedDBPropagator& propagator);

    /**
     * Add an independent set given as list of vertex indices.
     */
    inline void
    add_direct_independent_set(const std::vector<std::size_t>& indices);

    /**
     * Price a given range of vertices, and add those that could potentially
     * improve the optimal solution to the current relaxation.
     * Return true if it found any potentially interesting vertices, and false
     * if the current relaxation is optimal even if all priced vertices were
     * to be included in the current subgraph.
     */
    template <typename VertexIterator>
    inline std::size_t price_vertices(VertexIterator vbeg, VertexIterator vend);

    /**
     * Get the best mutually exclusive set found so far.
     */
    const std::vector<Vertex>& get_best_mes() const noexcept {
        return m_best_mes;
    }

    /**
     * Get the best solution, i.e., the best sample.
     */
    const std::vector<std::vector<bool>>& get_best_solution() const noexcept {
        return m_best_sample;
    }

    /**
     * Update the best (primal) solution.
     */
    void
    update_best_solution(const std::vector<std::vector<bool>>& new_solution) {
        if (new_solution.size() >= m_best_sample.size())
            return;
        m_best_sample = new_solution;
        if (m_subgraph_upper_bound > m_best_sample.size()) {
            m_subgraph_upper_bound = m_best_sample.size();
        }
    }

    double get_last_value() const noexcept { return m_last_objective_value; }

    std::size_t get_subgraph_ub() const noexcept {
        return m_subgraph_upper_bound;
    }

    inline void count_compressible_constraints(
        std::size_t min_iterations_without_relevance) {
        const auto ncmp = m_complete_independent_sets.size();
        const auto ninc = m_inc_independent_sets.size();
        std::size_t ccmp = 0, cinc = 0;
        for (std::size_t i = 0; i < ncmp; ++i) {
            if (m_complete_independent_set_zero_since[i] >=
                min_iterations_without_relevance)
            {
                ++ccmp;
            }
        }
        for (std::size_t i = 0; i < ninc; ++i) {
            if (m_inc_independent_set_zero_since[i] >=
                min_iterations_without_relevance)
            {
                ++cinc;
            }
        }
        double percent_reduction = 100.0 * (ccmp + cinc) / double(ncmp + ninc);
        m_recorder->store_event("COUNTED_COMPRESSIBLE",
                                {{"incomplete_reducible", cinc},
                                 {"complete_reducible", ccmp},
                                 {"percent_reducible", percent_reduction}},
                                "incomplete_reducible", "complete_reducible",
                                "percent_reducible");
    }

    /**
     * Try to greedily add new vertices to existing independent sets,
     * strengthening the corresponding constraints. Usually should be
     * tried before generating new independent sets.
     * Returns true if one of the strengthened constraints is now violated.
     */
    inline bool greedy_add_to_cutting_planes();

    /**
     * Try to greedily generate new, violated independent set constraints.
     * Returns true if a violated new constraint could be found.
     */
    inline bool greedy_generate_cutting_planes();

    /**
     * Incorporate a mutually exclusive set found externally.
     */
    inline void external_improved_mes(std::vector<Vertex> external_mes);

    /**
     * Add the (up to) num_constraints incomplete or complete
     * independent sets with the highest dual weights to the given
     * primal solver as initial partial configurations.
     */
    inline void
    export_highest_weights_to_primal(ColoringHeuristicSolver& primal,
                                     std::size_t num_constraints);

    /**
     * Get the support of the current fractional solution,
     * i.e., the vertices with non-zero value.
     */
    inline std::vector<Vertex> get_fractional_mes_support() const;

    /**
     * Use a vector of complete configurations to try to find violated cutting
     * planes.
     */
    inline std::size_t
    cuts_from_primal(const std::vector<SharedDBPropagator>& primal_configs,
                     std::size_t max_num_cuts);

    /**
     * Use a sample of complete configurations to try to find violated cutting
     * planes.
     */
    inline std::size_t
    cuts_from_sample(const std::vector<std::vector<bool>>& sample);

  private:
    // subgraph we are working on
    UniverseSubgraph* m_subgraph;
    // event recorder
    EventRecorder* m_recorder;
    // configuration values
    LowerBoundMIPConfig m_config;

    // gurobi model
    GRBModel m_model;

    // collection of as-of-yet incomplete independent sets
    // with buffers for their dual values and constraint handles
    std::unique_ptr<double[]> m_inc_independent_dual{nullptr};
    std::vector<GRBConstr> m_inc_independent_set_constraints;
    std::vector<SharedDBPropagator> m_inc_independent_sets;
    std::vector<std::uint32_t> m_inc_independent_set_zero_since;

    // collection of complete independent sets (i.e., configurations)
    // with buffers for their dual values and constraint handles
    std::unique_ptr<double[]> m_complete_independent_dual{nullptr};
    std::vector<GRBConstr> m_complete_independent_set_constraints;
    std::vector<DynamicBitset> m_complete_independent_sets;
    std::vector<std::uint32_t> m_complete_independent_set_zero_since;

    /**
     * Encoding of a directly-represented independent set.
     */
    struct DirectIndependentSet {
        // vertices in the independent set
        HashSet<std::size_t> vertex_indices;

        // literals are incompatible if any vertex
        // implies their negation
        DynamicBitset compatible_literals;
    };

    // collection of directly-represented independent sets
    // with buffers for their dual values and constraint handles
    std::unique_ptr<double[]> m_direct_independent_dual{nullptr};
    std::vector<GRBConstr> m_direct_independent_set_constraints;
    std::vector<DirectIndependentSet> m_direct_independent_sets;

    // number of features
    Var m_nall;

    // vector of all variables and buffer of primal values in last solution
    std::vector<GRBVar> m_vars;
    std::unique_ptr<double[]> m_last_solution_values{nullptr};

    using IndependentSetBuilder = UniverseSubgraph::IndependentSetBuilder;
    using CliqueBuilder = UniverseSubgraph::CliqueBuilder;
    IndependentSetBuilder m_indset_builder;
    CliqueBuilder m_clique_builder;

    // the currently extracted solution (vertex_index, value pairs), sorted by
    // decreasing value
    using OrderedSolution = std::vector<std::pair<std::size_t, double>>;
    using OSolIter = OrderedSolution::const_iterator;
    OrderedSolution m_ordered_solution;

    // the last objective value
    double m_last_objective_value;
    // the best integral value we can hope for on the current subgraph
    std::size_t m_subgraph_upper_bound;
    // if we have added new constraints by cutting, we want to solve
    // the relaxation to have dual values
    bool m_have_new_constraints = false;

    // the best sample we currently know
    std::vector<std::vector<bool>> m_best_sample;
    // the best mutually exclusive set we currently know
    std::vector<Vertex> m_best_mes;

    // buffer for a gurobi linear expression
    GRBLinExpr m_expr_buffer;

    // buffer for violated nonedges
    std::vector<std::pair<OSolIter, OSolIter>> m_violated_nonedges;

    // possibly, a callback if we're doing a MIP solve
    // or an abortable solve
    std::unique_ptr<GRBCallback> m_callback;
    std::unique_ptr<GRBCallback> m_old_callback;

    // mutex for when we need it during parallel operations
    std::mutex ds_lock;

    // if this solver is made abortable,
    // m_abortable is true. if the search
    // shall be aborted, m_aborted is set to true
    bool m_abortable{false};
    std::atomic<bool> m_aborted{false};

    // number of calls to optimize()
    std::size_t m_num_solve_calls{0};

    /**
     * Initialize the vector of variables.
     */
    inline void p_create_initial_vars();

    /**
     * Execute callable once for each literal.
     */
    template <typename Callable>
    inline void p_foreach_literal(Callable&& callable);

    const std::vector<std::vector<bool>>&
    p_to_sample(const std::vector<std::vector<bool>>& s) {
        return s;
    }

    std::vector<std::vector<bool>>
    p_to_sample(const std::vector<DynamicBitset>& s) {
        std::vector<std::vector<bool>> result;
        std::transform(s.begin(), s.end(), std::back_inserter(result),
                       [](const DynamicBitset& b) {
                           return static_cast<std::vector<bool>>(b);
                       });
        return result;
    }

    /**
     * Extract an ordered solution from the model, i.e.,
     * refresh/replace m_ordered_solution.
     */
    inline void
    p_extract_ordered_solution(bool possible_suboptimal_mip = false);

    /**
     * Turn a fractional bound on the objective value into an integer one,
     * leaving a small gap for numerical issues.
     */
    inline std::size_t
    p_make_objective_value_integral(double fractional) const noexcept;

    /**
     * Try to use rounding to turn the current fractional solution into a
     * clique.
     */
    inline void p_ordered_solution_to_clique();

    inline void p_begin_mip();
    inline void p_end_mip();

    /**
     * Identify violated nonedges.
     */
    inline bool p_identify_violated_nonedges();

    /**
     * Attempt to prohibit a violated non-edge by extending an
     * existing constraint to cover that non-edge.
     */
    inline bool p_prohibit_extend_existing(std::size_t v_index,
                                           DynamicBitset& covered_vne);

    /**
     * Update the lists of dual values.
     */
    inline void p_update_dual_values();

    /**
     * Add the given vertices to independent sets.
     * If the score of vertices increases to a safe level,
     * may delete some vertices from the list.
     */
    inline void
    p_add_vertices_to_independent_sets(std::vector<Vertex>& vertices,
                                       VertexMapTo<double>& vertex_scores);

    /**
     * Compute the total dual value of all existing
     * independent sets already containing v.
     */
    inline double p_price_in_existing(Vertex v) const noexcept;

    /**
     * Refresh an incomplete independent set constraint,
     * i.e., regenerate it based on its independent set
     * (to handle the case of indirectly added vertices).
     * Called in greedy constraint strengthening
     * when we detect an independent set to be violated.
     */
    inline void p_refresh_inc_constraint(std::size_t inc_index);

    /**
     * Collect vertices (either 'good' vertices, which correspond
     * to a dual constraint violated by at least 0.01, or 'possible'
     * vertices which are not explicitly in the current subgraph
     * but are only within 0.01 of corresponding to a violated dual constraint).
     */
    template <typename VertexIterator>
    inline bool
    p_collect_good_or_possible(VertexIterator vbeg, VertexIterator vend,
                               std::vector<Vertex>& good_vertices,
                               std::vector<Vertex>& possible_vertices,
                               VertexMapTo<double>& interesting_scores);

    /**
     * Categorize the given good vertices depending on their reduced costs,
     * and use that information to reduce the number of good vertices.
     */
    inline void p_categorize_good_vertices(
        std::vector<Vertex>& good, VertexMapTo<double>& interesting_scores,
        std::size_t num_categories, std::size_t goal_vertices);

    /**
     * Extract a complete independent set (constraint) and add it as initially
     * present configuration to the given primal heuristic solver.
     */
    inline void p_complete_to_primal(ColoringHeuristicSolver& primal,
                                     const DynamicBitset& complete_indset);

    /**
     * Extract an incomplete independent set (constraint) and add it as
     * initially present configuration to the given primal heuristic solver.
     */
    inline void p_incomplete_to_primal(ColoringHeuristicSolver& primal,
                                       const SharedDBPropagator& propagator);

    /**
     * Generate directly-represented independent sets.
     * Used if cheap cut procedures fail.
     */
    inline bool p_greedy_generate_direct_independent_sets();

    /**
     * Check if a given index can be added to a direct independent set
     * that it is not already part of.
     */
    inline bool p_direct_is_compatible(const DirectIndependentSet& direct,
                                       std::size_t index);

    /**
     * Strengthen direct independent set constraints
     * by adding vertices to them, based on the current
     * relaxed value of the vertices.
     */
    inline std::size_t p_greedy_add_to_direct();
};

namespace detail {

inline bool is_in_complete(const DynamicBitset& c, Vertex v) {
    using sammy::lit::negative;
    using sammy::lit::var;
    return c[var(v.first)] != negative(v.first) &&
           c[var(v.second)] != negative(v.second);
}

inline bool is_in_incomplete(const SharedDBPropagator& ic, Vertex v) {
    return ic.is_true(v.first) && ic.is_true(v.second);
}

} // namespace detail

void GurobiCliqueSolverG2::p_create_initial_vars() {
    const auto n = m_subgraph->n();
    m_vars.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        m_vars.emplace_back(m_model.addVar(0.0, 1.0, 1.0, GRB_CONTINUOUS));
    }
}

template <typename ConfigurationIterator>
inline void
GurobiCliqueSolverG2::add_configurations(ConfigurationIterator cbeg,
                                         ConfigurationIterator cend) {
    using ConfigType =
        typename std::iterator_traits<ConfigurationIterator>::value_type;
    std::for_each(cbeg, cend,
                  [this](const ConfigType& cfg) { add_configuration(cfg); });
}

template <typename ConfigurationType>
inline void
GurobiCliqueSolverG2::add_configuration(const ConfigurationType& cfg) {
    auto literal_value = [&cfg](Lit l) {
        return cfg[lit::var(l)] != lit::negative(l);
    };
    m_have_new_constraints = true;
    m_complete_independent_dual.reset();
    m_complete_independent_sets.emplace_back(cfg);
    m_complete_independent_set_zero_since.push_back(0);
    const auto n = m_subgraph->n();
    m_expr_buffer.clear();
    for (std::size_t vi = 0; vi < n; ++vi) {
        Vertex v = m_subgraph->vertex(vi);
        if (literal_value(v.first) && literal_value(v.second)) {
            m_expr_buffer += m_vars[vi];
        }
    }
    m_complete_independent_set_constraints.push_back(
        m_model.addConstr(m_expr_buffer <= 1));
}

template <typename Callable>
inline void GurobiCliqueSolverG2::p_foreach_literal(Callable&& callable) {
    for (Lit l = 0, nall_lit = 2 * m_nall; l < nall_lit; ++l) {
        std::invoke(std::forward<Callable>(callable), l);
    }
}

SolverState GurobiCliqueSolverG2::solve_relaxation(double time_allowed) {
    if (m_subgraph_upper_bound <= m_best_mes.size() && !m_have_new_constraints)
    {
        return SolverState::OPTIMUM_ON_SUBGRAPH;
    }
    if (time_allowed <= 0.0) {
        return SolverState::TIMEOUT_NO_IMPROVEMENT;
    }
    if (std::isfinite(time_allowed)) {
        m_model.set(GRB_DoubleParam_TimeLimit, time_allowed);
    }
    // m_recorder->store_event("BEGIN_SOLVE_RELAXATION");
    m_model.optimize();
    ++m_num_solve_calls;
    m_have_new_constraints = false;
    if (m_model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
        m_recorder->store_event("TIMEOUT_SOLVE_RELAXATION");
        return SolverState::TIMEOUT_NO_IMPROVEMENT;
    }
    std::size_t prev_clique_size = m_best_mes.size();
    p_extract_ordered_solution();
    // m_recorder->store_event("DONE_SOLVE_RELAXATION", {{"bound",
    // m_last_objective_value}}, "bound");
    p_ordered_solution_to_clique();
    std::size_t new_clique_size = m_best_mes.size();
    bool is_now_opt = (new_clique_size >= m_subgraph_upper_bound);
    return is_now_opt ? SolverState::OPTIMUM_ON_SUBGRAPH
                      : (new_clique_size > prev_clique_size
                             ? SolverState::IMPROVEMENT_FOUND
                             : SolverState::NO_IMPROVEMENT_FOUND);
}

void GurobiCliqueSolverG2::p_extract_ordered_solution(
    bool possible_suboptimal_mip) {
    m_ordered_solution.clear();
    m_last_solution_values.reset(
        m_model.get(GRB_DoubleAttr_X, m_vars.data(), int(m_vars.size())));
    const auto n = m_subgraph->n();
    for (std::size_t i = 0; i < n; ++i) {
        if (m_last_solution_values[i] > 1.0e-6) {
            m_ordered_solution.emplace_back(i, m_last_solution_values[i]);
        }
    }
    std::sort(
        m_ordered_solution.begin(), m_ordered_solution.end(),
        [](const auto& p1, const auto& p2) { return p1.second > p2.second; });
    if (possible_suboptimal_mip) {
        m_last_objective_value = m_model.get(GRB_DoubleAttr_ObjBound);
        m_subgraph_upper_bound =
            p_make_objective_value_integral(m_last_objective_value);
    } else {
        m_last_objective_value = m_model.get(GRB_DoubleAttr_ObjVal);
        m_subgraph_upper_bound =
            p_make_objective_value_integral(m_last_objective_value);
        p_update_dual_values();
    }
}

std::size_t
GurobiCliqueSolverG2::p_make_objective_value_integral(double v) const noexcept {
    auto i = std::make_signed_t<std::size_t>(std::floor(v + 0.01));
    if (i < 0)
        return 0;
    return std::size_t(i);
}

void GurobiCliqueSolverG2::p_ordered_solution_to_clique() {
    if (m_subgraph_upper_bound > m_best_mes.size()) {
        m_clique_builder.reset_vertices();
        m_clique_builder.greedily_extend_weighted(m_ordered_solution.begin(),
                                                  m_ordered_solution.end());
        m_clique_builder.randomly_make_maximal(sammy::rng());
        if (m_clique_builder.size() > m_best_mes.size()) {
            m_recorder->store_event(
                "IMPROVED_LB",
                {{"lb", m_clique_builder.size()}, {"method", "LP rounding"}},
                "lb", "method");
            m_best_mes = m_clique_builder.get_vertices();
        }
    }
}

SolverState GurobiCliqueSolverG2::solve_full_relaxation(double time_allowed) {
    if (m_subgraph_upper_bound == m_best_mes.size() && !m_have_new_constraints)
    {
        return SolverState::OPTIMUM_ON_SUBGRAPH;
    }
    if (time_allowed <= 0.0) {
        return SolverState::TIMEOUT_NO_IMPROVEMENT;
    }
    m_recorder->store_event(
        "BEGIN_SOLVE_FULL_RELAXATION",
        {{"incomplete_sets", m_inc_independent_sets.size()},
         {"complete_sets", m_complete_independent_sets.size()},
         {"n", m_subgraph->n()}},
        "incomplete_sets", "complete_sets", "n");
    auto before = Clock::now();
    bool any_improvement_found = false;
    for (;;) {
        double time_remaining =
            time_allowed - seconds_between(before, Clock::now());
        switch (solve_relaxation(time_remaining)) {
        case SolverState::IMPROVEMENT_FOUND:
            any_improvement_found = true;
            break;
        case SolverState::NO_IMPROVEMENT_FOUND:
            break;
        case SolverState::OPTIMUM_ON_SUBGRAPH:
            m_recorder->store_event("DONE_SOLVE_FULL_RELAXATION",
                                    {{"bound", m_subgraph_upper_bound},
                                     {"status", "optimum on subgraph"},
                                     {"best_mes", m_best_mes.size()}},
                                    "bound", "status", "best_mes");
            return SolverState::OPTIMUM_ON_SUBGRAPH;
        case SolverState::TIMEOUT_NO_IMPROVEMENT:
            m_recorder->store_event("TIMEOUT_SOLVE_FULL_RELAXATION");
            return any_improvement_found ? SolverState::TIMEOUT_IMPROVEMENT
                                         : SolverState::TIMEOUT_NO_IMPROVEMENT;
        case SolverState::TIMEOUT_IMPROVEMENT:
            m_recorder->store_event("TIMEOUT_SOLVE_FULL_RELAXATION");
            return SolverState::TIMEOUT_IMPROVEMENT;
        }
        if (!prohibit_violated_nonedges()) {
            m_recorder->store_event("DONE_SOLVE_FULL_RELAXATION",
                                    {{"bound", m_subgraph_upper_bound},
                                     {"best_mes", m_best_mes.size()}},
                                    "bound", "best_mes");
            return any_improvement_found ? SolverState::IMPROVEMENT_FOUND
                                         : SolverState::NO_IMPROVEMENT_FOUND;
        }
    }
}

/**
 * Identify violated non-edges.
 */
bool GurobiCliqueSolverG2::p_identify_violated_nonedges() {
    m_violated_nonedges.clear();
    for (OSolIter c = m_ordered_solution.begin(), e = m_ordered_solution.end();
         c != e; ++c)
    {
        double weight = c->second;
        std::size_t index = c->first;
        double thresh = 1.01 - weight;
        if (thresh > weight)
            break;
        const auto& row_c = m_subgraph->matrix_row(index);
        for (OSolIter d = std::next(c); d != e; ++d) {
            double weight_d = d->second;
            if (weight_d < thresh)
                break;
            std::size_t index_d = d->first;
            if (!row_c[index_d])
                m_violated_nonedges.emplace_back(c, d);
        }
    }
    return !m_violated_nonedges.empty();
}

bool GurobiCliqueSolverG2::p_prohibit_extend_existing(
    std::size_t v_index, DynamicBitset& covered_vne) {
    auto handle_other = [&](SharedDBPropagator& target_prop, std::size_t vio) {
        if (covered_vne[vio])
            return;
        auto [vit2, wit2] = m_violated_nonedges[vio];
        Vertex v = m_subgraph->vertex(vit2->first),
               w = m_subgraph->vertex(wit2->first);
        if (push_vertex_pair(target_prop, v, w)) {
            covered_vne[vio] = true;
        }
    };
    auto [vit, wit] = m_violated_nonedges[v_index];
    for (std::size_t ind_index = 0, n_inc = m_inc_independent_sets.size();
         ind_index < n_inc; ++ind_index)
    {
        std::size_t vi = vit->first, wi = wit->first;
        Vertex v = m_subgraph->vertex(vi), w = m_subgraph->vertex(wi);
        SharedDBPropagator& target_prop = m_inc_independent_sets[ind_index];
        if (target_prop.is_false(v.first) || target_prop.is_false(v.second) ||
            target_prop.is_false(w.first) || target_prop.is_false(w.second))
        {
            continue;
        }
        if (push_vertex_pair(target_prop, v, w)) {
            covered_vne[v_index] = true;
            auto nvio = m_violated_nonedges.size();
            for (std::size_t neindex = v_index + 1; neindex < nvio; ++neindex) {
                handle_other(target_prop, neindex);
            }
            p_refresh_inc_constraint(ind_index);
            return true;
        }
    }
    return false;
}

bool GurobiCliqueSolverG2::prohibit_violated_nonedges() {
    if (!p_identify_violated_nonedges()) {
        return false;
    }
    DynamicBitset covered_vne(m_violated_nonedges.size(), false);
    OSolIter last_first = m_ordered_solution.end();
    const auto nve = m_violated_nonedges.size();
    SharedDBPropagator& prop = m_subgraph->get_propagator();
    auto handle_other_ve = [&](std::size_t j) {
        if (!covered_vne[j]) {
            auto [xit, yit] = m_violated_nonedges[j];
            Vertex x = m_subgraph->vertex(xit->first);
            Vertex y = m_subgraph->vertex(yit->first);
            if (push_vertex_pair(prop, x, y)) {
                covered_vne[j].set();
            }
        }
    };
    for (std::size_t i = 0; i < nve; ++i) {
        auto [vit, wit] = m_violated_nonedges[i];
        if (vit != last_first && !covered_vne[i]) {
            if (p_prohibit_extend_existing(i, covered_vne))
                continue;
            last_first = vit;
            Vertex v = m_subgraph->vertex(vit->first);
            Vertex w = m_subgraph->vertex(wit->first);
            reset_and_push_noresolve(prop, v);
            if (push_vertex(prop, w) < 0) {
                // can happen if we've learned new clauses;
                // essentially, this is a 'hidden' edge in our graph
                covered_vne[i].set();
                m_subgraph->nonedge_to_edge(vit->first, wit->first);
                continue;
            }
            covered_vne[i].set();
            for (std::size_t j = i + 1; j < nve; ++j) {
                handle_other_ve(j);
            }
            for (std::size_t j = 0; j < i; ++j) {
                handle_other_ve(j);
            }
            add_potentially_incomplete_independent_set(prop);
        }
    }
    m_have_new_constraints = true;
    return true;
}

void GurobiCliqueSolverG2::add_direct_independent_set(
    const std::vector<std::size_t>& indices) {
    SharedDBPropagator& propagator = m_subgraph->get_propagator();
    DynamicBitset compatible_literals(2 * m_nall, true);
    for (std::size_t index : indices) {
        propagator.reset_or_throw();
        Vertex v = m_subgraph->vertex(index);
        if (push_vertex(propagator, v) < 0) {
            throw std::logic_error("Broken vertex in independent set!");
        }
        for (Lit l : propagator.get_trail()) {
            compatible_literals[lit::negate(l)].reset();
        }
    }
    m_direct_independent_dual.reset();
    m_direct_independent_sets.push_back(
        {HashSet<std::size_t>(indices.begin(), indices.end()),
         std::move(compatible_literals)});
    m_expr_buffer.clear();
    for (std::size_t index : indices) {
        m_expr_buffer += m_vars[index];
    }
    m_direct_independent_set_constraints.push_back(
        m_model.addConstr(m_expr_buffer <= 1));
}

void GurobiCliqueSolverG2::add_potentially_incomplete_independent_set(
    const SharedDBPropagator& propagator) {
    if (propagator.get_trail().size() == m_nall) {
        std::vector<bool> new_cfg(m_nall, false);
        for (Lit l : propagator.get_trail()) {
            if (!lit::negative(l)) {
                new_cfg[lit::var(l)] = true;
            }
        }
        add_configuration(new_cfg);
    } else {
        m_inc_independent_dual.reset();
        m_inc_independent_sets.push_back(propagator);
        m_inc_independent_set_zero_since.push_back(0);
        m_expr_buffer.clear();
        for (std::size_t vi = 0, n = m_subgraph->n(); vi != n; ++vi) {
            Vertex v = m_subgraph->vertex(vi);
            if (propagator.is_true(v.first) && propagator.is_true(v.second)) {
                m_expr_buffer += m_vars[vi];
            }
        }
        m_inc_independent_set_constraints.push_back(
            m_model.addConstr(m_expr_buffer <= 1));
    }
}

void GurobiCliqueSolverG2::add_vertices(const std::vector<Vertex>& vertices) {
    if (vertices.empty())
        return;
    m_subgraph_upper_bound = m_best_sample.size();
    const auto old_size = m_subgraph->n();
    m_subgraph->add_vertices(vertices);
    const auto new_size = m_subgraph->n();
    if (old_size + vertices.size() != new_size) {
        throw std::logic_error(
            "Duplicate or already existing vertices passed to add_vertices!");
    }
    for (Vertex v : vertices) {
        GRBColumn column;
        for (std::size_t i = 0, ni = m_complete_independent_sets.size();
             i != ni; ++i)
        {
            if (detail::is_in_complete(m_complete_independent_sets[i], v)) {
                column.addTerm(1.0, m_complete_independent_set_constraints[i]);
            }
        }
        for (std::size_t i = 0, ni = m_inc_independent_sets.size(); i != ni;
             ++i)
        {
            if (detail::is_in_incomplete(m_inc_independent_sets[i], v)) {
                column.addTerm(1.0, m_inc_independent_set_constraints[i]);
            }
        }
        m_vars.push_back(m_model.addVar(0.0, 1.0, 1.0, GRB_CONTINUOUS, column));
    }
}

void GurobiCliqueSolverG2::p_categorize_good_vertices(
    std::vector<Vertex>& good, VertexMapTo<double>& interesting_scores,
    std::size_t num_categories, std::size_t goal_vertices) {
    if (goal_vertices >= good.size())
        return;
    std::vector<std::size_t> category_count(num_categories, 0);
    const double delta_cat = 1.0 / num_categories;
    auto category_of = [&](double value) -> std::size_t {
        if (value <= 0)
            return 0;
        if (value >= 1)
            return num_categories - 1;
        return std::size_t(value / delta_cat);
    };
    for (Vertex v : good) {
        double val = interesting_scores.at(v);
        ++category_count[category_of(val)];
    }
    if (category_count.front() > goal_vertices) {
        goal_vertices = (std::max)(std::size_t(0.2 * category_count.front()),
                                   goal_vertices);
    }
    std::size_t take_all_from_below = 0;
    std::size_t total_count = 0;
    double last_ratio = 1.0;
    for (std::size_t c = 0; c < num_categories; ++c) {
        std::size_t this_cat = category_count[c];
        std::size_t count_after = total_count + this_cat;
        if (count_after <= goal_vertices) {
            ++take_all_from_below;
            total_count = count_after;
        } else {
            last_ratio = double(goal_vertices - total_count) / this_cat;
            break;
        }
    }
    auto& rng = sammy::rng();
    good.erase(std::remove_if(
                   good.begin(), good.end(),
                   [&](Vertex v) {
                       std::size_t cat = category_of(interesting_scores.at(v));
                       if (cat < take_all_from_below)
                           return false;
                       if (cat > take_all_from_below)
                           return true;
                       return bool(std::uniform_real_distribution<double>{
                                       0.0, 1.0}(rng) > last_ratio);
                   }),
               good.end());
}

template <typename VertexIterator>
std::size_t GurobiCliqueSolverG2::price_vertices(VertexIterator vbeg,
                                                 VertexIterator vend) {
    std::vector<Vertex> good_vertices, possible_vertices;
    VertexMapTo<double> interesting_scores;
    if (!p_collect_good_or_possible(vbeg, vend, good_vertices,
                                    possible_vertices, interesting_scores))
    {
        return 0;
    }
    std::size_t ret_result = good_vertices.size() + possible_vertices.size();
    if (!good_vertices.empty()) {
        std::size_t total_good = good_vertices.size();
        p_categorize_good_vertices(good_vertices, interesting_scores, 20,
                                   std::size_t(0.2 * m_subgraph->n()));
        std::size_t total_taken = good_vertices.size();
        p_add_vertices_to_independent_sets(good_vertices, interesting_scores);
        std::size_t added_to_independent_sets =
            total_taken - good_vertices.size();
        m_recorder->store_event(
            "ADDED_VERTICES",
            {{"added_to_graph", good_vertices.size()},
             {"added_to_indsets", added_to_independent_sets},
             {"considered", total_taken},
             {"total_good", total_good},
             {"type", "good"}},
            "added_to_graph", "added_to_indsets", "considered", "total_good");
        add_vertices(good_vertices);
    } else {
        std::size_t analyzed = possible_vertices.size();
        p_add_vertices_to_independent_sets(possible_vertices,
                                           interesting_scores);
        m_recorder->store_event("ADDED_VERTICES",
                                {{"new_count", possible_vertices.size()},
                                 {"total_count", analyzed},
                                 {"type", "possible"}},
                                "new_count", "total_count", "type");
        add_vertices(possible_vertices);
    }
    return ret_result;
}

void GurobiCliqueSolverG2::p_update_dual_values() {
    m_complete_independent_dual.reset(m_model.get(
        GRB_DoubleAttr_Pi, m_complete_independent_set_constraints.data(),
        int(m_complete_independent_set_constraints.size())));
    m_inc_independent_dual.reset(
        m_model.get(GRB_DoubleAttr_Pi, m_inc_independent_set_constraints.data(),
                    int(m_inc_independent_set_constraints.size())));
    m_direct_independent_dual.reset(m_model.get(
        GRB_DoubleAttr_Pi, m_direct_independent_set_constraints.data(),
        int(m_direct_independent_set_constraints.size())));
    for (std::size_t i = 0; i < m_complete_independent_sets.size(); ++i) {
        if (m_complete_independent_dual[i] > 1.0e-5) {
            m_complete_independent_set_zero_since[i] = 0;
        } else {
            ++m_complete_independent_set_zero_since[i];
        }
    }
    for (std::size_t i = 0; i < m_inc_independent_sets.size(); ++i) {
        if (m_inc_independent_dual[i] > 1.0e-5) {
            m_inc_independent_set_zero_since[i] = 0;
        } else {
            ++m_inc_independent_set_zero_since[i];
        }
    }
}

double GurobiCliqueSolverG2::p_price_in_existing(Vertex v) const noexcept {
    double score = 0.0;
    for (std::size_t i = 0, ni = m_complete_independent_sets.size(); i < ni;
         ++i)
    {
        if (detail::is_in_complete(m_complete_independent_sets[i], v)) {
            score += m_complete_independent_dual[i];
            if (score >= 1.02)
                return 1.02;
        }
    }
    for (std::size_t i = 0, ni = m_inc_independent_sets.size(); i < ni; ++i) {
        if (detail::is_in_incomplete(m_inc_independent_sets[i], v)) {
            score += m_inc_independent_dual[i];
            if (score >= 1.02)
                return 1.02;
        }
    }
    return score;
}

void GurobiCliqueSolverG2::p_add_vertices_to_independent_sets(
    std::vector<Vertex>& vertices, VertexMapTo<double>& vertex_scores) {
    for (Vertex v : vertices) {
        double& weight = vertex_scores.at(v);
        for (std::size_t i = 0, ni = m_inc_independent_sets.size(); i != ni;
             ++i)
        {
            SharedDBPropagator& prop = m_inc_independent_sets[i];
            if (!detail::is_in_incomplete(prop, v) && push_vertex(prop, v) >= 0)
            {
                weight += m_inc_independent_dual[i];
                if (weight >= 1.01)
                    break;
            }
        }
    }
    auto has_high_weight = [&vertex_scores](Vertex v) {
        return vertex_scores.at(v) >= 1.01;
    };
    vertices.erase(
        std::remove_if(vertices.begin(), vertices.end(), has_high_weight),
        vertices.end());
}

bool GurobiCliqueSolverG2::greedy_add_to_cutting_planes() {
    std::size_t num_refreshed = 0;
    std::size_t num_strengthened = 0;
    for (std::size_t i = 0, ni = m_inc_independent_sets.size(); i != ni; ++i) {
        double could_add = 0.0, current_weight = 0.0;
        SharedDBPropagator& prop = m_inc_independent_sets[i];
        for (const auto& e : m_ordered_solution) {
            std::size_t vi = e.first;
            double wi = e.second;
            Vertex v = m_subgraph->vertex(vi);
            if (detail::is_in_incomplete(prop, v)) {
                current_weight += wi;
            } else if (can_push(prop, v)) {
                could_add += wi;
            }
        }
        if (current_weight > 1.01) {
            // we have, possibly indirectly, added
            // present vertices to this independent set
            p_refresh_inc_constraint(i);
            ++num_strengthened;
            ++num_refreshed;
        } else if (current_weight + could_add > 1.01) {
            for (const auto& e : m_ordered_solution) {
                Vertex v = m_subgraph->vertex(e.first);
                if (detail::is_in_incomplete(prop, v)) {
                    continue;
                }
                if (push_vertex(prop, v) >= 0) {
                    current_weight += e.second;
                    m_model.chgCoeff(m_inc_independent_set_constraints[i],
                                     m_vars[e.first], 1.0);
                }
            }
            if (current_weight > 1.01) {
                ++num_strengthened;
            }
        }
    }
    if (num_strengthened == 0) {
        num_strengthened = p_greedy_add_to_direct();
    }
    if (num_strengthened > 0) {
        m_have_new_constraints = true;
        m_recorder->store_event("STRENGTHENED_CONSTRAINTS",
                                {{"num_strengthened", num_strengthened},
                                 {"num_refreshed", num_refreshed}},
                                "num_strengthened", "num_refreshed");
    }
    return num_strengthened > 0;
}

void GurobiCliqueSolverG2::p_refresh_inc_constraint(std::size_t inc_index) {
    m_expr_buffer.clear();
    const SharedDBPropagator& prop = m_inc_independent_sets[inc_index];
    for (std::size_t vi = 0, n = m_subgraph->n(); vi != n; ++vi) {
        Vertex v = m_subgraph->vertex(vi);
        if (detail::is_in_incomplete(prop, v)) {
            m_expr_buffer += m_vars[vi];
        }
    }
    GRBConstr constraint = m_model.addConstr(m_expr_buffer <= 1);
    assert(inc_index < m_inc_independent_set_constraints.size());
    m_model.remove(m_inc_independent_set_constraints[inc_index]);
    m_inc_independent_set_constraints[inc_index] = constraint;
}

std::size_t GurobiCliqueSolverG2::p_greedy_add_to_direct() {
    std::size_t num_strengthened = 0;
    auto& propagator = m_subgraph->get_propagator();
    for (std::size_t direct_index = 0, dn = m_direct_independent_sets.size();
         direct_index < dn; ++direct_index)
    {
        auto& direct = m_direct_independent_sets[direct_index];
        double current_weight = std::accumulate(
            m_ordered_solution.begin(), m_ordered_solution.end(), 0.0,
            [&](double acc, const auto& e) {
                if (direct.vertex_indices.count(e.first)) {
                    return acc + e.second;
                }
                return acc;
            });
        for (auto [index, value] : m_ordered_solution) {
            if (direct.vertex_indices.count(index))
                continue;
            Vertex v = m_subgraph->vertex(index);
            if (!direct.compatible_literals[v.first] ||
                !direct.compatible_literals[v.second])
            {
                continue;
            }
            if (!p_direct_is_compatible(direct, index)) {
                continue;
            }
            current_weight += value;
            direct.vertex_indices.insert(index);
            reset_and_push_noresolve(propagator, v);
            for (Lit l : propagator.get_trail()) {
                direct.compatible_literals[lit::negate(l)].reset();
            }
            m_model.chgCoeff(m_direct_independent_set_constraints[direct_index],
                             m_vars[index], 1.0);
        }
        if (current_weight > 1.01) {
            ++num_strengthened;
        }
    }
    return num_strengthened;
}

bool GurobiCliqueSolverG2::p_greedy_generate_direct_independent_sets() {
    if (m_ordered_solution.size() < 6) {
        return false;
    }
    auto ordered_solution_copy = m_ordered_solution;
    std::size_t initial_part_size = m_ordered_solution.size() / 3;
    auto initial_part_begin = ordered_solution_copy.begin();
    auto second_part_begin = initial_part_begin + initial_part_size;
    auto second_part_end = ordered_solution_copy.end();
    auto& rng = sammy::rng();
    double largest_weight = 0.0;
    std::vector<std::size_t> largest_set;
    m_indset_builder.reset_vertices();
    largest_weight = m_indset_builder.greedily_extend_weighted(
        initial_part_begin, second_part_end);
    largest_set = m_indset_builder.get_indices();
    for (std::size_t it = 0; it < 100; ++it) {
        std::shuffle(initial_part_begin, second_part_begin, rng);
        std::shuffle(second_part_begin, second_part_end, rng);
        m_indset_builder.reset_vertices();
        double weight = m_indset_builder.greedily_extend_weighted(
            initial_part_begin, second_part_end);
        weight += m_indset_builder.greedily_extend_weighted(second_part_begin,
                                                            second_part_end);
        if (weight > 1.01) {
            if (weight > largest_weight) {
                largest_weight = weight;
                largest_set = m_indset_builder.get_indices();
            }
        }
    }
    if (largest_weight > 1.01) {
        add_direct_independent_set(largest_set);
        return true;
    }
    // if nothing is found, return false
    return false;
}

bool GurobiCliqueSolverG2::greedy_generate_cutting_planes() {
    using IEntry = std::pair<std::vector<Vertex>, double>;
    SharedDBPropagator& propagator = m_subgraph->get_propagator();
    std::vector<IEntry> violated;
    std::vector<Vertex> current_pushed;
    double added_weights = 0.0;
    auto handle_entry = [&](std::pair<std::size_t, double> e) {
        Vertex v = m_subgraph->vertex(e.first);
        if (push_vertex(propagator, v) >= 0) {
            added_weights += e.second;
            current_pushed.push_back(v);
        }
    };
    for (OSolIter cbegin = m_ordered_solution.begin(),
                  cend = m_ordered_solution.end();
         cbegin != cend; ++cbegin)
    {
        added_weights = 0.0;
        current_pushed.clear();
        propagator.reset_or_throw();
        for (const auto& e : IteratorRange{cbegin, cend})
            handle_entry(e);
        for (const auto& e : IteratorRange{m_ordered_solution.cbegin(), cbegin})
            handle_entry(e);
        if (added_weights >= 1.01) {
            violated.emplace_back(current_pushed, added_weights);
        }
    }
    if (violated.empty()) {
        return p_greedy_generate_direct_independent_sets();
    }
    auto pos = violated.end();
    if (violated.size() > 10) {
        pos = violated.begin() + 10;
        std::nth_element(violated.begin(), pos, violated.end(),
                         [](const IEntry& e1, const IEntry& e2) {
                             return e1.second > e2.second;
                         });
    }
    double max_vio = 1.0;
    for (auto i = violated.begin(); i != pos; ++i) {
        if (i->second > max_vio)
            max_vio = i->second;
        propagator.reset_or_throw();
        for (Vertex v : i->first) {
            push_vertex(propagator, v);
        }
        add_potentially_incomplete_independent_set(propagator);
    }
    m_recorder->store_event(
        "GREEDY_CUTTING_PLANES",
        {{"count", pos - violated.begin()}, {"max_violation", max_vio - 1.0}},
        "count", "max_violation");
    m_have_new_constraints = true;
    return true;
}

template <typename VertexIterator>
bool GurobiCliqueSolverG2::p_collect_good_or_possible(
    VertexIterator vbeg, VertexIterator vend,
    std::vector<Vertex>& good_vertices, std::vector<Vertex>& possible_vertices,
    VertexMapTo<double>& interesting_scores) {
    assert(m_inc_independent_sets.empty() ||
           m_inc_independent_dual.get() != nullptr);
    assert(m_complete_independent_sets.empty() ||
           m_complete_independent_dual.get() != nullptr);
    struct Context {
        bool have_good = false;
        bool many_good = false;
    };
    const auto n = m_subgraph->n();
    thread_group().parallel_foreach(
        vbeg, vend, Context{}, [&](Context& ctx, Vertex v) {
            // if(ctx.many_good) return;
            if (m_subgraph->has_vertex(v))
                return;
            double ex_score = p_price_in_existing(v);
            if (ex_score >= 1.01)
                return;
            if (ex_score >= 0.99 && ctx.have_good)
                return;
            std::unique_lock l{ds_lock};
            ctx.have_good = !good_vertices.empty();
            ctx.many_good = (good_vertices.size() > 0.2 * n);
            // if(ctx.many_good) return;
            if (ex_score >= 0.99 && ctx.have_good)
                return;
            if (ex_score < 0.99) {
                ctx.have_good = true;
                if (!possible_vertices.empty()) {
                    possible_vertices.clear();
                    interesting_scores.clear();
                }
                good_vertices.push_back(v);
                interesting_scores[v] = ex_score;
            } else {
                possible_vertices.push_back(v);
                interesting_scores[v] = ex_score;
            }
        });
    return !good_vertices.empty() || !possible_vertices.empty();
}

void GurobiCliqueSolverG2::export_highest_weights_to_primal(
    ColoringHeuristicSolver& primal, std::size_t num_constraints) {
    struct ConstraintInfo {
        ConstraintInfo(std::size_t i, bool c, double d) noexcept
            : index(i), is_complete(c), dual_weight(d) {}

        std::size_t index;
        bool is_complete;
        double dual_weight;

        bool operator<(const ConstraintInfo& cinfo) const noexcept {
            return dual_weight > cinfo.dual_weight;
        }
    };

    std::vector<ConstraintInfo> cinfo;
    const auto nc = m_complete_independent_sets.size(),
               ni = m_inc_independent_sets.size();
    cinfo.reserve(nc + ni);
    for (std::size_t ic = 0; ic < nc; ++ic) {
        if (m_complete_independent_dual[ic] > 0.01) {
            cinfo.emplace_back(ic, true, m_complete_independent_dual[ic]);
        }
    }
    for (std::size_t ii = 0; ii < ni; ++ii) {
        if (m_inc_independent_dual[ii] > 0.01) {
            cinfo.emplace_back(ii, false, m_inc_independent_dual[ii]);
        }
    }
    auto end = cinfo.end();
    if (num_constraints < cinfo.size()) {
        end = cinfo.begin() + num_constraints;
        std::nth_element(cinfo.begin(), end, cinfo.end());
    }
    for (ConstraintInfo& ci : IteratorRange{cinfo.begin(), end}) {
        if (ci.is_complete) {
            p_complete_to_primal(primal, m_complete_independent_sets[ci.index]);
        } else {
            p_incomplete_to_primal(primal, m_inc_independent_sets[ci.index]);
        }
    }
}

void GurobiCliqueSolverG2::p_complete_to_primal(
    ColoringHeuristicSolver& primal, const DynamicBitset& complete_indset) {
    SharedDBPropagator propagator = m_subgraph->get_propagator();
    propagator.reset_or_throw();
    for (Var v = 0; v < complete_indset.size(); ++v) {
        Lit l =
            complete_indset[v] ? lit::positive_lit(v) : lit::negative_lit(v);
        if (propagator.is_false(l) ||
            (!propagator.is_true(l) && !propagator.push_level(l)))
        {
            throw std::logic_error(
                "Failed to push literal from complete configuration!");
        }
    }
    primal.add_color_class(std::move(propagator));
}

void GurobiCliqueSolverG2::p_incomplete_to_primal(
    ColoringHeuristicSolver& primal, const SharedDBPropagator& propagator) {
    primal.add_color_class(propagator);
}

std::vector<Vertex> GurobiCliqueSolverG2::get_fractional_mes_support() const {
    std::vector<Vertex> result;
    for (const auto& e : m_ordered_solution) {
        if (e.second > 0.01) {
            result.push_back(m_subgraph->vertex(e.first));
        } else
            break;
    }
    return result;
}

std::size_t GurobiCliqueSolverG2::cuts_from_primal(
    const std::vector<SharedDBPropagator>& primal_configs,
    std::size_t max_num_cuts) {
    using CVElement = std::pair<std::size_t, double>;
    std::vector<CVElement> cut_values;
    for (std::size_t i = 0, ncfg = primal_configs.size(); i < ncfg; ++i) {
        const SharedDBPropagator& prop = primal_configs[i];
        double total = 0.0;
        for (const auto& e : m_ordered_solution) {
            Vertex v = m_subgraph->vertex(e.first);
            if (prop.is_true(v.first) && prop.is_true(v.second)) {
                total += e.second;
            }
        }
        if (total > 1.01) {
            cut_values.emplace_back(i, total);
        }
    }
    if (cut_values.empty())
        return 0;
    auto end = cut_values.end();
    if (cut_values.size() > max_num_cuts) {
        end = cut_values.begin() + max_num_cuts;
        std::nth_element(cut_values.begin(), end, cut_values.end(),
                         [&](const CVElement& p1, const CVElement& p2) {
                             return p1.second > p2.second;
                         });
    }
    for (const CVElement& e : IteratorRange{cut_values.begin(), end}) {
        add_potentially_incomplete_independent_set(primal_configs[e.first]);
    }
    m_have_new_constraints = true;
    return std::size_t(end - cut_values.begin());
}

SolverState GurobiCliqueSolverG2::solve_as_mip(double time_allowed) {
    if (time_allowed <= 0)
        return SolverState::TIMEOUT_NO_IMPROVEMENT;
    if (std::isfinite(time_allowed)) {
        m_model.set(GRB_DoubleParam_TimeLimit, time_allowed);
    }
    p_begin_mip();
    m_model.optimize();
    m_num_solve_calls++;
    auto status = m_model.get(GRB_IntAttr_Status);
    auto solcount = m_model.get(GRB_IntAttr_SolCount);
    if (solcount > 0) {
        std::size_t old_mes_size = get_best_mes().size();
        p_extract_ordered_solution(true);
        p_ordered_solution_to_clique();
        p_end_mip();
        if (status == GRB_OPTIMAL)
            return SolverState::OPTIMUM_ON_SUBGRAPH;
        return get_best_mes().size() > old_mes_size
                   ? SolverState::TIMEOUT_IMPROVEMENT
                   : SolverState::TIMEOUT_NO_IMPROVEMENT;
    }
    p_end_mip();
    return SolverState::TIMEOUT_NO_IMPROVEMENT;
}

void GurobiCliqueSolverG2::p_begin_mip() {
    class Callback : public GRBCallback {
      public:
        explicit Callback(GurobiCliqueSolverG2* that) noexcept : that(that) {}

        void callback() override {
            if (that->m_abortable && that->m_aborted.load()) {
                GRBCallback::abort();
                return;
            }
            if (where != GRB_CB_MIPSOL) {
                return;
            }
            std::unique_ptr<double[]> sol_values{this->getSolution(
                that->m_vars.data(), int(that->m_vars.size()))};
            sol_vertices.clear();
            for (std::size_t i = 0, n = that->m_vars.size(); i < n; ++i) {
                if (sol_values[i] > 0.5)
                    sol_vertices.push_back(i);
            }
            for (std::size_t vi = 0, vn = sol_vertices.size(); vi < vn; ++vi) {
                std::size_t v = sol_vertices[vi];
                for (std::size_t vj = 0; vj < vi; ++vj) {
                    std::size_t w = sol_vertices[vj];
                    if (!that->m_subgraph->is_edge(v, w)) {
                        this->addLazy(that->m_vars[v] + that->m_vars[w] <= 1);
                    }
                }
            }
        }

      private:
        GurobiCliqueSolverG2* that;
        std::vector<std::size_t> sol_vertices;
    };
    for (GRBVar v : m_vars) {
        v.set(GRB_CharAttr_VType, GRB_BINARY);
    }
    m_old_callback.reset(m_callback.release());
    m_callback = std::make_unique<Callback>(this);
    m_model.set(GRB_IntParam_LazyConstraints, 1);
    m_model.setCallback(m_callback.get());
}

void GurobiCliqueSolverG2::p_end_mip() {
    for (GRBVar v : m_vars) {
        v.set(GRB_CharAttr_VType, GRB_CONTINUOUS);
    }
    m_model.set(GRB_IntParam_LazyConstraints, 0);
    m_callback.reset(m_old_callback.release());
    m_model.setCallback(m_callback.get());
}

void GurobiCliqueSolverG2::external_improved_mes(
    std::vector<Vertex> external_mes) {
    if (external_mes.size() <= m_best_mes.size())
        return;
    add_new_vertices(external_mes);
    m_best_mes = std::move(external_mes);
}

std::size_t GurobiCliqueSolverG2::cuts_from_sample(
    const std::vector<std::vector<bool>>& sample) {
    auto literal_value = [&](const std::vector<bool>& config, Lit l) -> bool {
        bool v = config[lit::var(l)];
        return lit::negative(l) ^ v;
    };
    auto vertex_value = [&](const std::vector<bool>& config, Vertex v) {
        return literal_value(config, v.first) &&
               literal_value(config, v.second);
    };
    std::size_t result = 0;
    for (std::size_t i = 0, ncfg = sample.size(); i < ncfg; ++i) {
        const auto& config = sample[i];
        double total = 0.0;
        for (const auto& e : m_ordered_solution) {
            Vertex v = m_subgraph->vertex(e.first);
            if (vertex_value(config, v))
                total += e.second;
        }
        if (total > 1.01) {
            add_configuration(config);
            ++result;
        }
    }
    return result;
}

bool GurobiCliqueSolverG2::p_direct_is_compatible(
    const DirectIndependentSet& direct, std::size_t index) {
    return std::all_of(direct.vertex_indices.begin(),
                       direct.vertex_indices.end(),
                       [&](std::size_t other_index) {
                           return !m_subgraph->is_edge(index, other_index);
                       });
}

} // namespace sammy

#endif
