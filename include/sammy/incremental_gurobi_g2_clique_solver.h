#ifndef SAMMY_INCREMENTAL_GUROBI_G2_CLIQUE_SOLVER_H_INCLUDED_
#define SAMMY_INCREMENTAL_GUROBI_G2_CLIQUE_SOLVER_H_INCLUDED_

#include "algorithm_ex.h"
#include "best_k.h"
#include "gurobi.h"
#include "lazy_g2_adjacency_matrix.h"
#include "literals.h"
#include "thread_interrupt.h"

namespace sammy {

/**
 * Store the relaxed solution of a Gurobi clique model,
 * both ordered by value and accessible by index.
 */
class IncrementalGurobiCliqueRelaxedSolution {
  public:
    using CliqueModelIndex = std::size_t;

    struct OrderedValue {
        CliqueModelIndex clique_model_index;
        double value;

        bool operator<(const OrderedValue& other) const noexcept {
            return value > other.value;
        }
    };

    void update(GRBModel& model, const std::vector<GRBVar>& variables,
                const std::vector<GRBConstr>& constraints, bool primal_only) {
        m_nonzeros.clear();
        m_values.reset(model.get(GRB_DoubleAttr_X, variables.data(),
                                 static_cast<int>(variables.size())));
        if (!primal_only) {
            m_dual.reset(model.get(GRB_DoubleAttr_Pi, constraints.data(),
                                   static_cast<int>(constraints.size())));
        }
        for (CliqueModelIndex i = 0, s = variables.size(); i < s; ++i) {
            double value = m_values[i];
            if (value > 1e-4) {
                m_nonzeros.push_back({i, value});
            }
        }
        std::sort(m_nonzeros.begin(), m_nonzeros.end());
    }

    const std::vector<OrderedValue>& nonzeros() const noexcept {
        return m_nonzeros;
    }

    double operator[](CliqueModelIndex index) const noexcept {
        return m_values[index];
    }

    const double* dual_values() const noexcept { return m_dual.get(); }

  private:
    std::unique_ptr<double[]> m_values;
    std::unique_ptr<double[]> m_dual;
    std::vector<OrderedValue> m_nonzeros;
};

/**
 * Gurobi callback that only checks for aborts.
 */
class GurobiCallbackCheckAbort : public GRBCallback {
  public:
    GurobiCallbackCheckAbort() = default;

    void callback() override {
        if (get_and_clear_interrupt_flag()) {
            this->abort();
        }
    }
};

class ExplicitIndependentSet {
  public:
    ExplicitIndependentSet(
        const LazyG2AdjacencyMatrix* matrix,
        const std::vector<std::size_t>& full_indices,
        const std::vector<std::size_t>& clique_model_vertices)
        : matrix(matrix), vertices(full_indices.begin(), full_indices.end()),
          clique_model_vertices(clique_model_vertices),
          possible_extensions(matrix->n(), true) {
        for (std::size_t vi : vertices) {
            possible_extensions -= matrix->row(vi);
            possible_extensions[vi].reset();
        }
    }

    bool contains(Vertex v) const {
        return vertices.count(matrix->index_of(v));
    }

    bool contains(std::size_t full_vertex_index) const {
        return vertices.count(full_vertex_index);
    }

    const HashSet<std::size_t>& full_vertex_indices() const noexcept {
        return vertices;
    }

    const std::vector<std::size_t>& clique_model_indices() const noexcept {
        return clique_model_vertices;
    }

    bool extendable(std::size_t full_vertex_index) const noexcept {
        return possible_extensions[full_vertex_index];
    }

    bool extend_if_possible(std::size_t full_vertex_index,
                            std::size_t clique_model_index) noexcept {
        if (!possible_extensions[full_vertex_index]) {
            return false;
        }
        vertices.insert(full_vertex_index);
        clique_model_vertices.push_back(clique_model_index);
        possible_extensions -= matrix->row(full_vertex_index);
        possible_extensions[full_vertex_index].reset();
        return true;
    }

    bool extend_if_possible(std::size_t fv1, std::size_t ci1, std::size_t fv2,
                            std::size_t ci2) noexcept {
        if (!possible_extensions[fv1] || !possible_extensions[fv2]) {
            return false;
        }
        vertices.insert(fv1);
        clique_model_vertices.push_back(ci1);
        vertices.insert(fv2);
        clique_model_vertices.push_back(ci2);
        possible_extensions -= matrix->row(fv1);
        possible_extensions -= matrix->row(fv2);
        possible_extensions[fv1].reset();
        possible_extensions[fv2].reset();
        return true;
    }

  private:
    const LazyG2AdjacencyMatrix* matrix;
    HashSet<std::size_t> vertices;
    std::vector<std::size_t> clique_model_vertices;
    DynamicBitset possible_extensions;
};

/**
 * A constraint in the incremental Gurobi clique model;
 * can be either represented by a 'complete' assignment (bitset),
 * or by a propagator that can be used to track and change the assignment.
 */
class IncrementalGurobiCliqueConstraint {
  public:
    explicit IncrementalGurobiCliqueConstraint(const DynamicBitset& bitset)
        : m_representation(bitset) {}

    explicit IncrementalGurobiCliqueConstraint(DynamicBitset&& bitset) noexcept
        : m_representation(std::move(bitset)) {}

    explicit IncrementalGurobiCliqueConstraint(
        const SharedDBPropagator& propagator)
        : m_representation(propagator) {}

    explicit IncrementalGurobiCliqueConstraint(
        SharedDBPropagator&& propagator) noexcept
        : m_representation(std::move(propagator)) {}

    explicit IncrementalGurobiCliqueConstraint(
        const LazyG2AdjacencyMatrix* matrix,
        const std::vector<std::size_t>& vertices,
        const std::vector<std::size_t>& clique_model_vertices)
        : m_representation(std::in_place_type<ExplicitIndependentSet>, matrix,
                           vertices, clique_model_vertices) {}

    bool is_complete() const noexcept {
        return std::holds_alternative<DynamicBitset>(m_representation);
    }

    bool is_partial() const noexcept {
        return std::holds_alternative<SharedDBPropagator>(m_representation);
    }

    bool is_explicit() const noexcept {
        return std::holds_alternative<ExplicitIndependentSet>(m_representation);
    }

    bool contains(Vertex v) const {
        return std::visit(
            overloaded{[&](const SharedDBPropagator& prop) -> bool {
                           return prop.is_true(v.first) &&
                                  prop.is_true(v.second);
                       },
                       [&](const DynamicBitset& bitset) -> bool {
                           return lit::is_true_in(v.first, bitset) &&
                                  lit::is_true_in(v.second, bitset);
                       },
                       [&](const ExplicitIndependentSet& set) -> bool {
                           return set.contains(v);
                       }},
            m_representation);
    }

    SharedDBPropagator& partial() {
        return *std::get_if<SharedDBPropagator>(&m_representation);
    }

    const SharedDBPropagator& partial() const {
        return *std::get_if<SharedDBPropagator>(&m_representation);
    }

    ExplicitIndependentSet& explicit_set() {
        return *std::get_if<ExplicitIndependentSet>(&m_representation);
    }

    const ExplicitIndependentSet& explicit_set() const {
        return *std::get_if<ExplicitIndependentSet>(&m_representation);
    }

    GRBConstr& constraint() noexcept { return m_constraint; }

    const GRBConstr& constraint() const noexcept { return m_constraint; }

  private:
    GRBConstr m_constraint;
    std::variant<SharedDBPropagator, DynamicBitset, ExplicitIndependentSet>
        m_representation;
};

class IncrementalCliqueGreedyRounder {
  public:
    using FullVertexIndex = std::size_t;

    explicit IncrementalCliqueGreedyRounder(LazyG2AdjacencyMatrix* adj)
        : m_adj(adj), m_vertex_bitset_buffer(adj->n()) {}

    void round_single_iteration(
        const IncrementalGurobiCliqueRelaxedSolution& solution,
        const std::vector<FullVertexIndex>& clique_model_vertices) {
        p_reset();
        p_partial_clique_from_ordered_solution(solution, clique_model_vertices);
        p_randomly_extend_clique();
    }

    const std::vector<FullVertexIndex>& clique() const noexcept {
        return m_vertex_set_buffer;
    }

    void
    randomly_extend_given_clique(const std::vector<FullVertexIndex>& clique) {
        p_reset();
        for (FullVertexIndex full_index : clique) {
            m_vertex_bitset_buffer &= m_adj->row(full_index);
            m_vertex_set_buffer.push_back(full_index);
        }
        p_randomly_extend_clique();
    }

  private:
    /**
     * Lazy G2.
     */
    LazyG2AdjacencyMatrix* m_adj;

    /**
     * Buffer for building vertex sets.
     */
    std::vector<FullVertexIndex> m_vertex_set_buffer;

    /**
     * Buffer for building vertex bitsets.
     */
    DynamicBitset m_vertex_bitset_buffer;

    /**
     * Buffer for random extension.
     */
    std::vector<FullVertexIndex> m_random_ext_buffer;

    void p_reset() {
        m_vertex_set_buffer.clear();
        m_vertex_bitset_buffer.set();
    }

    void p_partial_clique_from_ordered_solution(
        const IncrementalGurobiCliqueRelaxedSolution& solution,
        const std::vector<FullVertexIndex>& clique_model_vertices) {
        for (const auto& ordered_value : solution.nonzeros()) {
            FullVertexIndex full_index =
                clique_model_vertices[ordered_value.clique_model_index];
            if (!m_vertex_bitset_buffer[full_index]) {
                // not a possible extension of the clique
                continue;
            }
            m_vertex_bitset_buffer &= m_adj->row(full_index);
            m_vertex_set_buffer.push_back(full_index);
        }
    }

    void p_randomly_extend_clique() {
        p_randomly_extend_clique_by_trials(100);
        p_randomly_extend_clique_with_list();
    }

    void p_randomly_extend_clique_by_trials(std::size_t num_unsuccessful) {
        auto& rng = sammy::rng();
        const std::size_t n = m_adj->n();
        std::uniform_int_distribution<std::size_t> indices(0, n - 1);
        for (std::size_t trial_counter = 1; trial_counter <= num_unsuccessful;
             ++trial_counter)
        {
            std::size_t vi = indices(rng);
            if (!m_vertex_bitset_buffer[vi])
                continue;
            m_vertex_bitset_buffer &= m_adj->row(vi);
            m_vertex_set_buffer.push_back(vi);
            trial_counter = 0;
            assert(!m_vertex_bitset_buffer[vi]);
        }
    }

    void p_randomly_extend_clique_with_list() {
        for (FullVertexIndex vi : m_vertex_bitset_buffer.ones()) {
            m_random_ext_buffer.push_back(vi);
        }
        if (m_random_ext_buffer.empty())
            return;
        std::shuffle(m_random_ext_buffer.begin(), m_random_ext_buffer.end(),
                     sammy::rng());
        for (FullVertexIndex vi : m_random_ext_buffer) {
            if (!m_vertex_bitset_buffer[vi])
                continue;
            m_vertex_set_buffer.push_back(vi);
            m_vertex_bitset_buffer &= m_adj->row(vi);
            assert(!m_vertex_bitset_buffer[vi]);
        }
        m_random_ext_buffer.clear();
    }
};

class IncrementalGurobiG2CliqueSolver {
  public:
    using CliqueModelIndex = std::size_t;
    using FullVertexIndex = std::size_t;
    using OrderedValue = IncrementalGurobiCliqueRelaxedSolution::OrderedValue;
    using OrderedSolutionIter = std::vector<OrderedValue>::const_iterator;

    IncrementalGurobiG2CliqueSolver(
        LazyG2AdjacencyMatrix* adjacency,
        std::vector<FullVertexIndex> best_local_clique,
        const std::vector<DynamicBitset>& best_covering_assignment,
        bool gurobi_quiet = true,
        const std::unordered_map<std::string, std::string>& params = {})
        : m_model(sammy::gurobi_environment(gurobi_quiet)), m_adj(adjacency),
          m_all_clique_vars(adjacency->n()),
          m_best_clique(std::move(best_local_clique)),
          m_is_present(adjacency->n(), false),
          m_present_vertices_with_concrete_literal(2 *
                                                   adjacency->num_concrete()),
          m_lazy_nonedge_callback(this), m_greedy_rounder(adjacency),
          m_temp_propagator(adjacency->temp_propagator()),
          m_greedy_cut_candidates(p_get_param<std::size_t>(
              params, "num_greedy_cut_candidates", 7)) {
        m_pricing_ignore_possible = p_get_param<bool>(
            params, "pricing_ignore_possible", m_pricing_ignore_possible);
        m_greedy_cut_extension_max_num_cuts = p_get_param<std::size_t>(
            params, "greedy_cut_extension_max_num_cuts",
            m_greedy_cut_extension_max_num_cuts);
        m_greedy_cut_extension_undo_failed_extensions = p_get_param<bool>(
            params, "greedy_cut_extension_undo_failed_extensions",
            m_greedy_cut_extension_undo_failed_extensions);
        m_greedy_cut_generation_num_candidates =
            p_get_param<std::size_t>(params, "num_greedy_cut_candidates", 7);
        m_prohibit_vne_max_new_constraints =
            p_get_param<std::size_t>(params, "prohibit_vne_max_new_constraints",
                                     m_prohibit_vne_max_new_constraints);
        m_max_cheap_cut_rounds = p_get_param<std::size_t>(
            params, "max_cheap_cut_rounds", m_max_cheap_cut_rounds);
        m_pricing_goal_fraction_existing =
            p_get_param<double>(params, "pricing_goal_fraction_existing",
                                m_pricing_goal_fraction_existing);
        m_pricing_goal_min = p_get_param<std::size_t>(
            params, "pricing_goal_min", m_pricing_goal_min);
        m_pricing_goal_max = p_get_param<std::size_t>(
            params, "pricing_goal_max", m_pricing_goal_max);
        m_max_cheap_cut_required_improvement =
            p_get_param<double>(params, "max_cheap_cut_required_improvement",
                                m_max_cheap_cut_required_improvement);
        m_grb_method = p_get_param<int>(params, "grb_method", m_grb_method);
        m_grb_presolve =
            p_get_param<int>(params, "grb_presolve", m_grb_presolve);
        m_initial_random_local_cliques =
            p_get_param<std::size_t>(params, "initial_random_local_cliques",
                                     m_initial_random_local_cliques);
        m_cheap_cut_rounds_per_gap_check =
            p_get_param<std::size_t>(params, "cheap_cut_rounds_per_gap_check",
                                     m_cheap_cut_rounds_per_gap_check);
        m_temp_propagator.reset_to_zero();
        p_init_model_params();
        p_init_extend_local_clique();
        p_initial_random_local_cliques();
        for (const DynamicBitset& assignment : best_covering_assignment) {
            add_complete_assignment(assignment);
        }
    }

    /**
     * Add a new vertex from the full graph to
     * the clique model, if it is not already present.
     */
    void add_vertex(FullVertexIndex index) {
        assert(index < m_is_present.size());
        if (m_is_present[index])
            return;
        p_add_new_vertex(index);
        m_is_present[index].set();
    }

    /**
     * Add a range of vertices from the full graph to
     * the clique model, if they are not already present.
     */
    template <typename R> void add_vertices(const R& indices) {
        for (FullVertexIndex i : indices) {
            add_vertex(i);
        }
    }

    /**
     * Add a complete assignment, represented by a bitset, as constraint.
     */
    void add_complete_assignment(const DynamicBitset& bitset);

    /**
     * Add a partial assignment, represented by a propagator, as constraint.
     */
    void add_partial_assignment(SharedDBPropagator propagator);

    /**
     * Solve the relaxation, lazily including all non-edge constraints.
     * @throws InterruptError if the computation was interrupted.
     */
    void solve_full_relaxation();

    /**
     * Run pricing on the given range of vertices,
     * targeting the addition of goal_num_vertices
     * to the clique model.
     */
    template <typename Range>
    bool price_vertex_range(Range&& range, std::size_t goal_num_vertices);

    /**
     * Run pricing on the given range of vertices,
     * targeting the addition of goal_num_vertices
     * to the clique model.
     */
    template <typename Iterator>
    bool price_vertex_range(Iterator begin, Iterator end,
                            std::size_t goal_num_vertices) {
        return price_vertex_range(IteratorRange{begin, end}, goal_num_vertices);
    }

    /**
     * Return whether the current clique is optimal
     * for the current subset of vertices.
     */
    bool optimal_on_current_vertices() const noexcept {
        return m_best_clique.size() >= m_last_rounded_objective;
    }

    /**
     * Run the greedy rounding, cutting and pricing loop for
     * a given number of iterations. Price over the given
     * vertices.
     * Returns true if the best clique was improved and false otherwise.
     */
    bool run_iterations(std::size_t num_iterations,
                        const std::vector<FullVertexIndex>& vertices);

    /**
     * Get the best clique.
     */
    const std::vector<FullVertexIndex>& get_best_clique() const noexcept {
        return m_best_clique;
    }

    /**
     * Turn the problem into a MIP with lazy non-edge constraints and solve it.
     * Returns true if the best clique was improved and false otherwise.
     * @throws InterruptError if the computation was interrupted, but not before
     * updating the best clique, if a better one was found.
     */
    bool run_as_mip();

    LazyG2AdjacencyMatrix& adjacency_matrix() noexcept { return *m_adj; }

    const LazyG2AdjacencyMatrix& adjacency_matrix() const noexcept {
        return *m_adj;
    }

  private:
    /**
     * The actual gurobi clique model.
     */
    GRBModel m_model;

    /**
     * The adjacency matrix of the subgraph.
     * Contains a list indicating which vertex index
     * corresponds to which interaction.
     */
    LazyG2AdjacencyMatrix* m_adj;

    /**
     * Array of all (potential) clique vertex variables.
     * Index i corresponds to index i in m_all_vertices.
     * The variable is present if and only if the
     * corresponding bit in m_present_in_clique_model is set.
     */
    std::vector<GRBVar> m_all_clique_vars;

    /**
     * Array of all existing clique vertex variables.
     * Necessary for querying all values at once.
     */
    std::vector<GRBVar> m_existing_clique_vars;

    /**
     * Array (parallel to m_existing_clique_vars) of
     * the indices in the full graph of the vertices
     * present in the clique model.
     */
    std::vector<FullVertexIndex> m_full_vertex_indices;

    /**
     * List of vertices in the best clique found so far.
     */
    std::vector<FullVertexIndex> m_best_clique;

    /**
     * Bitset indicating which vertices are present in the clique model.
     */
    DynamicBitset m_is_present;

    /**
     * For each concrete literal, the list of vertex indices
     * into m_existing_clique_vars containing that literal
     * (that are present in the clique).
     */
    std::vector<std::vector<CliqueModelIndex>>
        m_present_vertices_with_concrete_literal;

    /**
     * The list of constraints.
     */
    std::vector<IncrementalGurobiCliqueConstraint> m_constraints;

    /**
     * The list of raw constraints.
     */
    std::vector<GRBConstr> m_raw_constraints;

    /**
     * Entry for pricing vertices.
     */
    struct PricingEntry {
        double
            dual_weight; // the LHS of the dual >= 1 constraint.
                         // the lower, the more violated is the dual constraint.
        FullVertexIndex full_index;

        bool operator<(const PricingEntry& other) const noexcept {
            return dual_weight < other.dual_weight;
        }
    };

    /**
     * Vertices with 'good' price, i.e., dual_weight <= 0.99.
     */
    std::vector<PricingEntry> m_pricing_good;

    /**
     * Vertices with 'possible' price, i.e., dual_weight <= 1.01.
     */
    std::vector<PricingEntry> m_pricing_possible;

    /**
     * Buffer for vertices selected by pricing.
     */
    std::vector<FullVertexIndex> m_pricing_selected;

    /**
     * Buffer for vertices in greedy cut under construction.
     */
    std::vector<FullVertexIndex> m_greedy_cut_buffer;

    /**
     * PARAMETER:
     * Ignore 'possible' vertices.
     * This usually does not change anything, but it
     * might cause very rare numerical issues where
     * vertices that may be in better cliques are ignored.
     */
    bool m_pricing_ignore_possible = false;

    /**
     * PARAMETER:
     * The maximum number of constraints to turn into cuts
     * during greedy constraint expansion to generate cuts.
     */
    std::size_t m_greedy_cut_extension_max_num_cuts = 10;

    /**
     * PARAMETER:
     * On performing greedy constraint expansion to generate cuts,
     * undo the partial extension if it does not lead to a cut.
     */
    bool m_greedy_cut_extension_undo_failed_extensions = true;

    /**
     * Buffer for creating a new column.
     */
    GRBColumn m_column_buffer;

    /**
     * Buffer for creating a new row.
     */
    GRBLinExpr m_expr_buffer;

    /**
     * Buffer for the latest relaxation solution.
     */
    IncrementalGurobiCliqueRelaxedSolution m_relaxed_solution;

    /**
     * The objective of the last (relaxed) solution.
     */
    double m_last_objective = std::numeric_limits<double>::infinity();

    /**
     * The rounded-down (with some numerical buffer)
     * objective of the last solution.
     * Used as UB on clique size.
     */
    std::size_t m_last_rounded_objective =
        std::numeric_limits<std::size_t>::max();

    /**
     * Abortion callback for the Gurobi model.
     */
    GurobiCallbackCheckAbort m_abort_callback;

    class GurobiCallbackMIPLazyNonEdge : public GRBCallback {
      public:
        explicit GurobiCallbackMIPLazyNonEdge(
            IncrementalGurobiG2CliqueSolver* that) noexcept
            : that(that) {}

        void callback() override {
            if (get_and_clear_interrupt_flag()) {
                this->abort();
                return;
            }
            if (where != GRB_CB_MIPSOL) {
                return;
            }
            std::unique_ptr<double[]> sol_values{
                this->getSolution(that->m_existing_clique_vars.data(),
                                  int(that->m_existing_clique_vars.size()))};
            m_sol_vertices.clear();
            for (std::size_t i = 0, n = that->m_existing_clique_vars.size();
                 i < n; ++i)
            {
                if (sol_values[i] > 0.5)
                    m_sol_vertices.push_back(i);
            }
            p_add_constraints();
        }

      private:
        void p_add_constraints() {
            for (auto i = m_sol_vertices.begin(), e = m_sol_vertices.end();
                 i != e; ++i)
            {
                FullVertexIndex fi = that->m_full_vertex_indices[*i];
                for (auto j = std::next(i); j != e; ++j) {
                    FullVertexIndex fj = that->m_full_vertex_indices[*j];
                    if (that->m_adj->is_definitive_nonedge(fi, fj)) {
                        this->addLazy(that->m_existing_clique_vars[*i] +
                                          that->m_existing_clique_vars[*j] <=
                                      1);
                    }
                }
            }
        }

        IncrementalGurobiG2CliqueSolver* that;
        std::vector<CliqueModelIndex> m_sol_vertices;
    };

    /**
     * Callback to add lazy non-edge constraints.
     */
    GurobiCallbackMIPLazyNonEdge m_lazy_nonedge_callback;

    /**
     * Rounding implementation.
     */
    IncrementalCliqueGreedyRounder m_greedy_rounder;

    /**
     * Collection buffer for violated non-edges.
     */
    std::vector<std::pair<OrderedSolutionIter, OrderedSolutionIter>>
        m_violated_nonedges;

    /**
     * Buffer for tracking the set of violated non-edges that we already
     * covered.
     */
    DynamicBitset m_covered_violated_nonedges;

    /**
     * Temporary propagator; emptied after each use.
     */
    SharedDBPropagator m_temp_propagator;

    /**
     * Represent a greedily generated candidate cut.
     */
    struct GreedyCutCandidate {
        std::vector<FullVertexIndex> inducing_vertices;
        double value;

        GreedyCutCandidate(double score_only) : value(score_only) {}

        GreedyCutCandidate(double score, std::vector<FullVertexIndex> vertices)
            : inducing_vertices(std::move(vertices)), value(score) {}

        bool operator<(const GreedyCutCandidate& other) const noexcept {
            return value > other.value;
        }
    };

    /**
     * PARAMETER:
     * Initially generated random local cliques.
     */
    std::size_t m_initial_random_local_cliques = 7;

    /**
     * PARAMETER:
     * Maximum number of greedy cuts to generate in each round.
     */
    std::size_t m_greedy_cut_generation_num_candidates = 7;

    /**
     * Buffer for the best greedy cut candidates during greedy cut generation.
     */
    BestK<GreedyCutCandidate> m_greedy_cut_candidates;

    /**
     * PARAMETER:
     * Maximum number of new constraints to generate
     * when prohibiting violated non-edges.
     */
    std::size_t m_prohibit_vne_max_new_constraints = 20;

    /**
     * PARAMETER:
     * Maximum total number of cheap cut rounds
     * before requiring a new pricing round.
     */
    std::size_t m_max_cheap_cut_rounds = 100;

    /**
     * PARAMETER:
     * The fraction of already added vertices
     * to target when determining how many vertices
     * to include during pricing.
     */
    double m_pricing_goal_fraction_existing = 0.05;

    /**
     * PARAMETER:
     * The minimum number of vertices to aim for during pricing.
     */
    std::size_t m_pricing_goal_min = 50;

    /**
     * PARAMETER:
     * The maximum number of vertices to aim for during pricing.
     */
    std::size_t m_pricing_goal_max = 1000;

    /**
     * PARAMETER:
     * Minimum absolute improvement in the relative gap
     * between the best clique and the relaxation objective
     * that is required every m_cheap_cut_rounds_per_gap_check
     * rounds to continue cheap cut rounds.
     */
    double m_max_cheap_cut_required_improvement = 0.02;

    /**
     * PARAMETER:
     * The method to use for solving the Gurobi LP models.
     */
    int m_grb_method = GRB_METHOD_DUAL;

    /**
     * PARAMETER:
     * The presolve level to use for solving the Gurobi LP models.
     */
    int m_grb_presolve = 1;

    /**
     * PARAMETER:
     * After how many rounds the gap is checked for
     * a sufficient improvement to continue cheap cut rounds.
     */
    std::size_t m_cheap_cut_rounds_per_gap_check = 10;

    /**
     * Buffer for the constraint index and extension score,
     * used during constraint extension for prohibiting violated non-edges.
     */
    std::vector<std::pair<std::size_t, int>>
        m_constraint_index_and_score_buffer;

    struct SwitchBack {
        explicit SwitchBack(IncrementalGurobiG2CliqueSolver* that) noexcept
            : that(that) {}

        ~SwitchBack() { that->p_switch_to_lp(); }

        IncrementalGurobiG2CliqueSolver* that;
    };

    /**
     * Add a new vertex.
     */
    void p_add_new_vertex(FullVertexIndex index);

    /**
     * Compute the initial column of the given vertex,
     * placing it into m_column_buffer.
     */
    void p_add_vertex_make_column(Vertex vertex);

    /**
     * Solve the relaxation of the clique model.
     * @throws InterruptError if the computation was interrupted.
     */
    void p_solve_relaxation();

    /**
     * Update the buffered relaxation solution.
     */
    void p_update_relaxed_solution(bool primal_only = false);

    /**
     * Set model parameters & attributes and callback.
     */
    void p_setup_model_params();

    /**
     * Run greedy rounding on the current linear relaxation;
     * returns true if the clique was extended and that led
     * to new vertices being added to the model.
     */
    bool p_round_greedy();

    /**
     * Make a list of violated non-edge constraints.
     * Return true if anything was found.
     */
    bool p_identify_violated_nonedges();

    /**
     * Generate or strengthen existing constraints
     * to prohibit non-edges violated by the current relaxed solution.
     */
    void p_prohibit_violated_nonedges();

    /**
     * Subroutine of p_prohibit_violated_nonedges:
     * search for partial constraints that can be extended
     * to cover the indicated violated non-edge.
     */
    bool p_prohibit_extend_existing(OrderedSolutionIter vit,
                                    OrderedSolutionIter wit);

    /**
     * Subrouting of p_prohibit_violated_nonedges:
     * search for explicit independent set constraints that
     * can be extended to cover the indicated violated non-edge.
     */
    bool p_prohibit_extend_explicit(OrderedSolutionIter vit,
                                    OrderedSolutionIter wit);

    /**
     * Called to update datastructure and constraints
     * after extending the partial constraint indexed by constraint_index
     * from its given old trail size.
     */
    void p_extended_partial(std::size_t constraint_index,
                            std::size_t old_trail_size);

    /**
     * Extend the partial assignment in m_temp_propagator
     * to greedily cover nonedges from m_violated_nonedges,
     * starting at index ibegin and ending at index iend.
     */
    void p_extend_nonedge_cover_greedy(std::size_t ibegin, std::size_t iend);

    /**
     * Collect good and possible vertices by pricing.
     */
    template <typename Range>
    bool p_price_range_collect_good_and_possible(Range&& range);

    /**
     * Total weight of dual LHS constraint of the given
     * vertex v (need not even be in all vertices).
     */
    double p_dual_lhs(Vertex v);

    /**
     * Select vertices to be added, targeting goal added vertices,
     * from the collected good and possible vertices.
     */
    void p_select_priced_vertices(std::size_t goal);

    /**
     * Find cutting planes for the current relaxation
     * using cheap heuristics.
     * Return true if a cut was found; raise
     * InterruptError if the computation was interrupted.
     * If a cut was found,the model is updated and the new
     * relaxation is solved.
     */
    bool p_cheap_cut_round();

    /**
     * Run cheap cut rounds until either no more cuts
     * are found (returning false) or until some termination
     * criterion is met (returning true).
     * The termination criteria are:
     *  - the best clique is optimal for the current subgraph,
     *  - the gap between the best clique and our bound did not improve by
     * m_cheap_cut_required_improvement in the last
     * m_cheap_cut_rounds_per_gap_check rounds.
     *  - a total of m_max_cheap_cut_rounds rounds was performed.
     * Raise an InterruptError if the computation was interrupted or timed out.
     */
    bool p_cheap_cut_rounds();

    /**
     * Return the relative gap of the current solution,
     * i.e., (REL-BEST)/BEST, where BEST is the current
     * best clique size and REL is the relaxation objective;
     * returns 0 if the current clique is optimal.
     */
    double p_relative_gap() const;

    /**
     * Greedily attempt to extend partial constraints
     * in order to achieve a violated cutting plane for
     * the current relaxation.
     * Return true if a cut was found; raise
     * InterruptError if the computation was interrupted.
     */
    bool p_greedy_cut_extend_constraints();

    /**
     * Attempt to expand the partial constraint given by index
     * so that it becomes a cutting plane for the current solution.
     * Returns true if a cutting plane was generated.
     */
    bool p_greedy_cut_attempt_expansion(std::size_t constraint_index);

    /**
     * Extend the partial assignment in m_temp_propagator
     * by pushing the vertices in subrange from [begin, end)
     * if possible.
     */
    double p_greedy_extend_cut(OrderedSolutionIter begin,
                               OrderedSolutionIter end);

    /**
     * Greedily attempt to generate a new cutting plane.
     */
    bool p_greedy_generate_cuts();

    /**
     * Randomly extend the current best clique.
     */
    void p_init_extend_local_clique();

    /**
     * Generate random local cliques to start with.
     * Expands the initial vertex set and potentially
     * leads to a better initial solution.
     */
    void p_initial_random_local_cliques();

    /**
     * Initialize the gurobi model.
     */
    void p_init_model_params();

    /**
     * Switch the problem from LP to IP.
     */
    SwitchBack p_switch_to_mip();

    /**
     * Switch the problem back to LP.
     */
    void p_switch_to_lp();

    /**
     * Greedily extend explicit independent set constraints
     * to generate cutting planes.
     */
    bool p_greedy_cut_extend_explicit();

    /**
     * Get a parameter from params, or a given default.
     */
    template <typename T, typename ParamMap>
    T p_get_param(const ParamMap& params, const std::string& name,
                  T default_value) const {
        auto found = params.find(name);
        if (found == params.end()) {
            return default_value;
        }
        std::istringstream data(found->second);
        T result;
        data >> result;
        return result;
    }
};

// -------- IMPLEMENTATION ----------
void IncrementalGurobiG2CliqueSolver::p_add_new_vertex(FullVertexIndex index) {
    Vertex vertex = m_adj->vertex(index);
    p_add_vertex_make_column(vertex);
    GRBVar var = m_model.addVar(0.0, 1.0, 1.0, GRB_CONTINUOUS, m_column_buffer);
    CliqueModelIndex clique_model_index = m_full_vertex_indices.size();
    m_full_vertex_indices.push_back(index);
    m_existing_clique_vars.push_back(var);
    m_all_clique_vars[index] = var;
    m_present_vertices_with_concrete_literal[vertex.first].push_back(
        clique_model_index);
    m_present_vertices_with_concrete_literal[vertex.second].push_back(
        clique_model_index);
}

void IncrementalGurobiG2CliqueSolver::p_add_vertex_make_column(Vertex vertex) {
    m_column_buffer.clear();
    for (const IncrementalGurobiCliqueConstraint& constraint : m_constraints) {
        if (!constraint.is_explicit() && constraint.contains(vertex)) {
            m_column_buffer.addTerm(1.0, constraint.constraint());
        }
    }
}

void IncrementalGurobiG2CliqueSolver::add_complete_assignment(
    const DynamicBitset& bitset) {
    m_expr_buffer.clear();
    m_constraints.emplace_back(bitset);
    IncrementalGurobiCliqueConstraint& constraint = m_constraints.back();
    for (CliqueModelIndex clique_index = 0, n = m_existing_clique_vars.size();
         clique_index < n; ++clique_index)
    {
        FullVertexIndex full_index = m_full_vertex_indices[clique_index];
        Vertex vertex = m_adj->vertex(full_index);
        if (lit::is_true_in(vertex.first, bitset) &&
            lit::is_true_in(vertex.second, bitset))
        {
            m_expr_buffer += m_existing_clique_vars[clique_index];
        }
    }
    constraint.constraint() = m_model.addConstr(m_expr_buffer <= 1.0);
    m_raw_constraints.push_back(constraint.constraint());
}

void IncrementalGurobiG2CliqueSolver::add_partial_assignment(
    SharedDBPropagator propagator_) {
    m_expr_buffer.clear();
    m_constraints.emplace_back(std::move(propagator_));
    IncrementalGurobiCliqueConstraint& constraint = m_constraints.back();
    const SharedDBPropagator& prop = constraint.partial();
    for (CliqueModelIndex clique_index = 0, n = m_existing_clique_vars.size();
         clique_index < n; ++clique_index)
    {
        FullVertexIndex full_index = m_full_vertex_indices[clique_index];
        Vertex vertex = m_adj->vertex(full_index);
        if (prop.is_true(vertex.first) && prop.is_true(vertex.second)) {
            m_expr_buffer += m_existing_clique_vars[clique_index];
        }
    }
    constraint.constraint() = m_model.addConstr(m_expr_buffer <= 1.0);
    m_raw_constraints.push_back(constraint.constraint());
}

void IncrementalGurobiG2CliqueSolver::p_update_relaxed_solution(
    bool primal_only) {
    m_relaxed_solution.update(m_model, m_existing_clique_vars,
                              m_raw_constraints, primal_only);
    m_last_objective = m_model.get(GRB_DoubleAttr_ObjVal);
    m_last_rounded_objective =
        static_cast<std::size_t>(std::floor(m_last_objective + 0.001));
}

void IncrementalGurobiG2CliqueSolver::p_setup_model_params() {
    m_model.set(GRB_IntParam_Threads, 1);
    m_model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    m_model.set(GRB_IntParam_Method, GRB_METHOD_DUAL);
    m_model.set(GRB_IntParam_Presolve, 1);
    m_model.setCallback(&m_abort_callback);
}

void IncrementalGurobiG2CliqueSolver::p_solve_relaxation() {
    bool repeat = true;
    while (repeat) {
        throw_if_interrupted();
        repeat = false;
        m_model.optimize();
        int status = m_model.get(GRB_IntAttr_Status);
        if (status == GRB_INTERRUPTED || status == GRB_TIME_LIMIT ||
            status == GRB_MEM_LIMIT || status == GRB_NODE_LIMIT ||
            status == GRB_NODE_LIMIT || status == GRB_ITERATION_LIMIT)
        {
            throw InterruptError();
        }
        if (status != GRB_OPTIMAL) {
            throw std::logic_error(
                "Unexpected status after solving relaxation!");
        }
        p_update_relaxed_solution();
        if (p_round_greedy()) {
            repeat = true;
        }
    }
}

bool IncrementalGurobiG2CliqueSolver::p_round_greedy() {
    m_greedy_rounder.round_single_iteration(m_relaxed_solution,
                                            m_full_vertex_indices);
    const auto& clique = m_greedy_rounder.clique();
    if (clique.size() <= m_best_clique.size()) {
        return false;
    }
    m_best_clique = clique;
    bool result = false;
    for (FullVertexIndex vi : clique) {
        if (m_is_present[vi])
            continue;
        result = true;
        add_vertex(vi);
    }
    return result;
}

void IncrementalGurobiG2CliqueSolver::solve_full_relaxation() {
    for (;;) {
        p_solve_relaxation();
        if (m_last_rounded_objective <= m_best_clique.size())
            return;
        if (!p_identify_violated_nonedges())
            return;
        p_prohibit_violated_nonedges();
    }
}

bool IncrementalGurobiG2CliqueSolver::p_identify_violated_nonedges() {
    m_violated_nonedges.clear();
    const auto& ordered_values = m_relaxed_solution.nonzeros();
    for (auto i = ordered_values.begin(), e = ordered_values.end(); i != e; ++i)
    {
        double weight = i->value;
        double thresh = 1.01 - weight;
        if (weight < thresh)
            break;
        CliqueModelIndex clique_index = i->clique_model_index;
        FullVertexIndex vertex_index = m_full_vertex_indices[clique_index];
        for (OrderedSolutionIter j = std::next(i); j != e; ++j) {
            double w_j = j->value;
            if (w_j < thresh)
                break;
            FullVertexIndex vertex_index2 =
                m_full_vertex_indices[j->clique_model_index];
            if (m_adj->is_definitive_nonedge(vertex_index, vertex_index2)) {
                m_violated_nonedges.emplace_back(i, j);
            }
        }
    }
    return !m_violated_nonedges.empty();
}

void IncrementalGurobiG2CliqueSolver::p_prohibit_violated_nonedges() {
    m_covered_violated_nonedges.assign(m_violated_nonedges.size(), false);
    auto last_first = m_relaxed_solution.nonzeros().cend();
    std::size_t num_new_constraints = 0;
    for (std::size_t i = 0, k = m_violated_nonedges.size(); i < k; ++i) {
        auto [vit, wit] = m_violated_nonedges[i];
        if (vit == last_first || m_covered_violated_nonedges[i]) {
            // do not generate or strengthen constraints if this non-edge
            continue;
        }
        if (p_prohibit_extend_existing(vit, wit) ||
            p_prohibit_extend_explicit(vit, wit))
        {
            // we have found a constraint we could strengthen to include this
            // non-edge
            m_covered_violated_nonedges[i].set();
            continue;
        }
        last_first = vit;
        if (++num_new_constraints > m_prohibit_vne_max_new_constraints) {
            // only generate constraints up to a certain number of new
            // constraints
            continue;
        }
        FullVertexIndex vind = m_full_vertex_indices[vit->clique_model_index];
        FullVertexIndex wind = m_full_vertex_indices[wit->clique_model_index];
        if (!push_vertex_pair(m_temp_propagator, m_adj->vertex(vind),
                              m_adj->vertex(wind)))
        {
            // this can rarely happen if we have new clauses in our clause DB
            // from conflict resolution; update the adjacency matrix.
            // as we have reached this point, we do not yet have a broken
            // constraint regarding this non-edge in our system.
            m_adj->nonedge_to_edge(vind, wind);
            m_covered_violated_nonedges[i].set();
            continue;
        }
        m_covered_violated_nonedges[i].set();
        p_extend_nonedge_cover_greedy(i + 1, k);
        p_extend_nonedge_cover_greedy(0, i);
        add_partial_assignment(m_temp_propagator);
        m_temp_propagator.reset_to_zero();
    }
}

bool IncrementalGurobiG2CliqueSolver::p_prohibit_extend_existing(
    OrderedSolutionIter vit, OrderedSolutionIter wit) {
    CliqueModelIndex ci1 = vit->clique_model_index;
    CliqueModelIndex ci2 = wit->clique_model_index;
    FullVertexIndex fi1 = m_full_vertex_indices[ci1];
    FullVertexIndex fi2 = m_full_vertex_indices[ci2];
    Vertex v1 = m_adj->vertex(fi1);
    Vertex v2 = m_adj->vertex(fi2);
    m_constraint_index_and_score_buffer.clear();
    for (std::size_t i = 0, ncons = m_constraints.size(); i < ncons; ++i) {
        if (!m_constraints[i].is_partial())
            continue;
        auto& prop = m_constraints[i].partial();
        if (prop.is_false(v1.first) || prop.is_false(v1.second) ||
            prop.is_false(v2.first) || prop.is_false(v2.second))
            continue;
        int score = int(prop.is_open(v1.first)) + prop.is_open(v1.second);
        if (v2.first != v1.first && v2.first != v1.second) {
            score += prop.is_open(v2.first);
        }
        if (v2.second != v1.first && v2.second != v1.second) {
            score += prop.is_open(v2.second);
        }
        if (score == 1) {
            std::size_t trail_length = prop.get_trail().size();
            if (push_vertex_pair(prop, v1, v2)) {
                p_extended_partial(i, trail_length);
                return true;
            } else {
                continue;
            }
        }
        m_constraint_index_and_score_buffer.emplace_back(i, score);
    }
    for (int goal_score = 2; goal_score < 5; ++goal_score) {
        for (auto entry : m_constraint_index_and_score_buffer) {
            if (entry.second != goal_score)
                continue;
            auto& prop = m_constraints[entry.first].partial();
            std::size_t trail_length = prop.get_trail().size();
            if (push_vertex_pair(prop, v1, v2)) {
                p_extended_partial(entry.first, trail_length);
                return true;
            }
        }
    }
    return false;
}

void IncrementalGurobiG2CliqueSolver::p_extended_partial(
    std::size_t constraint_index, std::size_t old_trail_size) {
    assert(m_constraints.size() > constraint_index);
    assert(m_constraints[constraint_index].is_partial());
    assert(m_constraints[constraint_index].partial().get_trail().size() >=
           old_trail_size);
    const SharedDBPropagator& prop = m_constraints[constraint_index].partial();
    const auto& trail = prop.get_trail();
    GRBConstr constr = m_constraints[constraint_index].constraint();
    const Lit nclit = 2 * m_adj->num_concrete();
    for (Lit lnew : IteratorRange{trail.begin() + old_trail_size, trail.end()})
    {
        if (lnew >= nclit)
            continue;
        for (CliqueModelIndex clique_model_index :
             m_present_vertices_with_concrete_literal[lnew])
        {
            FullVertexIndex full_index =
                m_full_vertex_indices[clique_model_index];
            Vertex v = m_adj->vertex(full_index);
            if (v.first == lnew) {
                if (prop.is_true(v.second)) {
                    m_model.chgCoeff(constr,
                                     m_existing_clique_vars[clique_model_index],
                                     1.0);
                }
            } else {
                if (prop.is_true(v.first)) {
                    m_model.chgCoeff(constr,
                                     m_existing_clique_vars[clique_model_index],
                                     1.0);
                }
            }
        }
    }
}

void IncrementalGurobiG2CliqueSolver::p_extend_nonedge_cover_greedy(
    std::size_t ibegin, std::size_t iend) {
    for (std::size_t j = ibegin; j < iend; ++j) {
        if (m_covered_violated_nonedges[j])
            continue;
        auto [vit, wit] = m_violated_nonedges[j];
        CliqueModelIndex ci1 = vit->clique_model_index;
        CliqueModelIndex ci2 = wit->clique_model_index;
        Vertex v1 = m_adj->vertex(m_full_vertex_indices[ci1]);
        Vertex v2 = m_adj->vertex(m_full_vertex_indices[ci2]);
        if (m_temp_propagator.is_false(v1.first) ||
            m_temp_propagator.is_false(v1.second) ||
            m_temp_propagator.is_false(v2.first) ||
            m_temp_propagator.is_false(v2.second))
        {
            continue;
        }
        if (push_vertex_pair(m_temp_propagator, v1, v2)) {
            m_covered_violated_nonedges[j].set();
        }
    }
}

template <typename Range>
bool IncrementalGurobiG2CliqueSolver::price_vertex_range(Range&& range,
                                                         std::size_t goal) {
    if (!p_price_range_collect_good_and_possible(std::forward<Range>(range))) {
        return false;
    }
    p_select_priced_vertices(goal);
    add_vertices(m_pricing_selected);
    solve_full_relaxation();
    return true;
}

template <typename Range>
bool IncrementalGurobiG2CliqueSolver::p_price_range_collect_good_and_possible(
    Range&& range) {
    m_pricing_good.clear();
    m_pricing_possible.clear();
    std::size_t abort_check = 0;
    for (FullVertexIndex fi : range) {
        if (m_is_present[fi])
            continue;
        double dual_weight = p_dual_lhs(m_adj->vertex(fi));
        if (dual_weight <= 0.99) {
            m_pricing_good.push_back({dual_weight, fi});
        } else if (!m_pricing_ignore_possible && dual_weight <= 1.01) {
            m_pricing_possible.push_back({dual_weight, fi});
        }
        if (++abort_check == 16384) {
            abort_check = 0;
            throw_if_interrupted();
        }
    }
    return !m_pricing_good.empty() || !m_pricing_possible.empty();
}

double IncrementalGurobiG2CliqueSolver::p_dual_lhs(Vertex v) {
    double value = 0.0;
    const std::size_t ncons = m_constraints.size();
    const double* dual_values = m_relaxed_solution.dual_values();
    for (std::size_t ci = 0; ci < ncons; ++ci) {
        if (m_constraints[ci].contains(v)) {
            value += dual_values[ci];
        }
    }
    return value;
}

void IncrementalGurobiG2CliqueSolver::p_select_priced_vertices(
    std::size_t goal) {
    m_pricing_selected.clear();

    auto select_from = [&](std::vector<PricingEntry>& entries) {
        if (goal >= entries.size()) {
            std::transform(
                entries.begin(), entries.end(),
                std::back_inserter(m_pricing_selected),
                [](const PricingEntry& entry) { return entry.full_index; });
        } else {
            auto relevant_end = entries.begin() + goal;
            std::nth_element(entries.begin(), relevant_end, entries.end());
            std::transform(
                entries.begin(), relevant_end,
                std::back_inserter(m_pricing_selected),
                [](const PricingEntry& entry) { return entry.full_index; });
        }
    };

    if (!m_pricing_good.empty()) {
        select_from(m_pricing_good);
    } else {
        select_from(m_pricing_possible);
    }
}

bool IncrementalGurobiG2CliqueSolver::p_cheap_cut_round() {
    if (p_greedy_cut_extend_explicit() || p_greedy_cut_extend_constraints() ||
        p_greedy_generate_cuts())
    {
        solve_full_relaxation();
        return true;
    }
    return false;
}

bool IncrementalGurobiG2CliqueSolver::p_greedy_cut_extend_explicit() {
    struct AchievableValue {
        double value;
        double existing_value;
        std::size_t constraint_index;
        bool operator<(const AchievableValue& other) const noexcept {
            return value > other.value;
        }
    };

    std::vector<AchievableValue> achievable_values;
    for (std::size_t i = 0, ncons = m_constraints.size(); i < ncons; ++i) {
        if (!m_constraints[i].is_explicit())
            continue;
        auto& eset = m_constraints[i].explicit_set();
        double total_value = 0.0;
        double existing_value;
        for (CliqueModelIndex i : eset.clique_model_indices()) {
            total_value += m_relaxed_solution[i];
        }
        existing_value = total_value;
        for (const auto& ov : m_relaxed_solution.nonzeros()) {
            std::size_t fi = m_full_vertex_indices[ov.clique_model_index];
            if (!eset.extendable(fi))
                continue;
            total_value += ov.value;
        }
        if (total_value <= 1.01) {
            continue;
        }
        achievable_values.push_back({total_value, existing_value, i});
    }

    std::sort(achievable_values.begin(), achievable_values.end());
    std::size_t num_cuts = 0;
    std::vector<std::pair<FullVertexIndex, CliqueModelIndex>> added_buffer;
    auto is_independent = [&](FullVertexIndex f) {
        for (auto entry : added_buffer) {
            if (!m_adj->is_definitive_nonedge(f, entry.first)) {
                return false;
            }
        }
        return true;
    };
    for (const AchievableValue& av : achievable_values) {
        if (num_cuts >= m_greedy_cut_extension_max_num_cuts) {
            break;
        }
        auto& constr = m_constraints[av.constraint_index];
        auto& eset = constr.explicit_set();
        added_buffer.clear();
        double value = av.existing_value;
        for (const auto& ov : m_relaxed_solution.nonzeros()) {
            std::size_t fi = m_full_vertex_indices[ov.clique_model_index];
            if (!eset.extendable(fi))
                continue;
            if (is_independent(fi)) {
                added_buffer.emplace_back(fi, ov.clique_model_index);
                value += ov.value;
            }
        }
        if (value >= 1.01) {
            for (auto entry : added_buffer) {
                if (!eset.extend_if_possible(entry.first, entry.second)) {
                    throw std::logic_error("Could not extend explicit set even "
                                           "though it should be possible!");
                }
                m_model.chgCoeff(constr.constraint(),
                                 m_existing_clique_vars[entry.second], 1.0);
            }
            num_cuts += 1;
        }
    }
    return num_cuts > 0;
}

bool IncrementalGurobiG2CliqueSolver::p_greedy_cut_extend_constraints() {
    struct AchievableValue {
        double value;
        std::size_t constraint_index;
        bool operator<(const AchievableValue& other) const noexcept {
            return value > other.value;
        }
    };

    std::vector<AchievableValue> achievable_values;
    for (std::size_t i = 0, ncons = m_constraints.size(); i < ncons; ++i) {
        if (!m_constraints[i].is_partial())
            continue;
        auto& prop = m_constraints[i].partial();
        double total_value = 0.0;
        for (const auto& ov : m_relaxed_solution.nonzeros()) {
            CliqueModelIndex ci = ov.clique_model_index;
            FullVertexIndex fi = m_full_vertex_indices[ci];
            Vertex v = m_adj->vertex(fi);
            if (!can_push(prop, v))
                continue;
            total_value += ov.value;
        }
        if (total_value >= 1.01) {
            achievable_values.push_back({total_value, i});
        }
    }
    std::sort(achievable_values.begin(), achievable_values.end());
    std::size_t num_cuts = 0;
    for (const AchievableValue& av : achievable_values) {
        if (num_cuts >= m_greedy_cut_extension_max_num_cuts) {
            break;
        }
        num_cuts += p_greedy_cut_attempt_expansion(av.constraint_index);
    }
    return num_cuts > 0;
}

bool IncrementalGurobiG2CliqueSolver::p_greedy_cut_attempt_expansion(
    std::size_t constraint_index) {
    assert(constraint_index < m_constraints.size());
    assert(m_constraints[constraint_index].is_partial());
    SharedDBPropagator& prop = m_constraints[constraint_index].partial();
    double total_value = 0.0;
    std::size_t old_trail_length = prop.get_trail().size();
    auto old_level = prop.get_current_level();
    for (const auto& ov : m_relaxed_solution.nonzeros()) {
        CliqueModelIndex ci = ov.clique_model_index;
        FullVertexIndex vi = m_full_vertex_indices[ci];
        Vertex v = m_adj->vertex(vi);
        if (push_vertex(prop, v) >= 0)
            total_value += ov.value;
    }
    if (total_value >= 1.01) {
        p_extended_partial(constraint_index, old_trail_length);
        return true;
    }
    if (m_greedy_cut_extension_undo_failed_extensions) {
        while (prop.get_current_level() > old_level) {
            prop.pop_level();
        }
    } else {
        p_extended_partial(constraint_index, old_trail_length);
    }
    return false;
}

double
IncrementalGurobiG2CliqueSolver::p_greedy_extend_cut(OrderedSolutionIter begin,
                                                     OrderedSolutionIter end) {
    double total = 0.0;
    for (OrderedSolutionIter i = begin; i != end; ++i) {
        CliqueModelIndex ci = i->clique_model_index;
        FullVertexIndex fi = m_full_vertex_indices[ci];
        Vertex v = m_adj->vertex(fi);
        if (push_vertex(m_temp_propagator, v) < 0) {
            continue;
        }
        m_greedy_cut_buffer.push_back(fi);
        total += i->value;
    }
    return total;
}

bool IncrementalGurobiG2CliqueSolver::p_greedy_generate_cuts() {
    const auto& nonzeros = m_relaxed_solution.nonzeros();
    m_greedy_cut_candidates.clear();
    for (OrderedSolutionIter b = nonzeros.begin(), i = b, e = nonzeros.end();
         i != e; ++i)
    {
        assert(m_temp_propagator.get_current_level() == 0);
        m_greedy_cut_buffer.clear();
        double v = p_greedy_extend_cut(i, e);
        v += p_greedy_extend_cut(b, i);
        if (v >= 1.01 &&
            m_greedy_cut_candidates.would_push(GreedyCutCandidate{v}))
        {
            m_greedy_cut_candidates.emplace(v, m_greedy_cut_buffer);
        }
        m_temp_propagator.reset_to_zero();
    }
    if (m_greedy_cut_candidates.elements().empty()) {
        return false;
    }
    for (const GreedyCutCandidate& candidate :
         m_greedy_cut_candidates.elements())
    {
        for (std::size_t vertex : candidate.inducing_vertices) {
            Vertex v = m_adj->vertex(vertex);
            if (push_vertex(m_temp_propagator, v) < 0) {
                throw std::logic_error(
                    "Incorrect greedy cut candidate: produced conflict!");
            }
        }
        add_partial_assignment(m_temp_propagator);
        m_temp_propagator.reset_to_zero();
    }
    return true;
}

double IncrementalGurobiG2CliqueSolver::p_relative_gap() const {
    if (optimal_on_current_vertices())
        return 0.0;
    double best = static_cast<double>(m_best_clique.size());
    return (m_last_objective - best) / best;
}

bool IncrementalGurobiG2CliqueSolver::p_cheap_cut_rounds() {
    if (optimal_on_current_vertices()) {
        return false;
    }
    std::size_t last_checked_at = 0;
    double last_gap = p_relative_gap();
    for (std::size_t num_rounds_total = 0;
         num_rounds_total < m_max_cheap_cut_rounds; ++num_rounds_total)
    {
        if (num_rounds_total - last_checked_at ==
            m_cheap_cut_rounds_per_gap_check)
        {
            double gap = p_relative_gap();
            if (gap > last_gap - m_max_cheap_cut_required_improvement) {
                break;
            }
            last_gap = gap;
            last_checked_at = num_rounds_total;
        }
        if (!p_cheap_cut_round()) {
            return false;
        }
        if (optimal_on_current_vertices()) {
            break;
        }
    }
    return true;
}

bool IncrementalGurobiG2CliqueSolver::run_iterations(
    std::size_t num_iterations, const std::vector<FullVertexIndex>& vertices) {
    auto pricing_goal = [&]() -> std::size_t {
        std::size_t fraction_of_existing =
            m_full_vertex_indices.size() * m_pricing_goal_fraction_existing;
        std::size_t hard_min = m_pricing_goal_min;
        std::size_t hard_max = m_pricing_goal_max;
        std::size_t value = (std::max)(fraction_of_existing, hard_min);
        value = (std::min)(value, hard_max);
        value = (std::min)(value, vertices.size());
        return value;
    };

    std::size_t old_clique_size = m_best_clique.size();
    if (m_relaxed_solution.nonzeros().empty()) {
        solve_full_relaxation();
    }
    bool cuts_are_stuck = false;
    for (std::size_t i = 0; i < num_iterations; ++i) {
        if (optimal_on_current_vertices() || cuts_are_stuck) {
            if (optimal_on_current_vertices()) {
            }
            if (!price_vertex_range(vertices, pricing_goal())) {
                break;
            }
            cuts_are_stuck =
                !optimal_on_current_vertices() && !p_cheap_cut_rounds();
        } else {
            cuts_are_stuck = !p_cheap_cut_rounds();
            if (!price_vertex_range(vertices, pricing_goal())) {
                if (cuts_are_stuck) {
                    break;
                }
            } else {
                cuts_are_stuck = false;
            }
        }
    }
    return m_best_clique.size() > old_clique_size;
}

void IncrementalGurobiG2CliqueSolver::p_init_extend_local_clique() {
    m_greedy_rounder.randomly_extend_given_clique(m_best_clique);
    add_vertices(m_greedy_rounder.clique());
}

void IncrementalGurobiG2CliqueSolver::p_initial_random_local_cliques() {
    auto& rng = sammy::rng();
    std::uniform_int_distribution<std::size_t> dist(0, m_adj->n() - 1);
    for (std::size_t i = 0; i < m_initial_random_local_cliques; ++i) {
        m_greedy_rounder.randomly_extend_given_clique({dist(rng)});
        const auto& clique = m_greedy_rounder.clique();
        add_vertices(clique);
        if (clique.size() > m_best_clique.size()) {
            m_best_clique = clique;
        }
    }
}

void IncrementalGurobiG2CliqueSolver::p_init_model_params() {
    m_model.set(GRB_IntParam_Threads, 1);
    m_model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    m_model.set(GRB_IntParam_Method, m_grb_method);
    m_model.set(GRB_IntParam_Presolve, m_grb_presolve);
    m_model.setCallback(&m_abort_callback);
}

bool IncrementalGurobiG2CliqueSolver::run_as_mip() {
    SwitchBack x = p_switch_to_mip();
    m_model.optimize();
    int status = m_model.get(GRB_IntAttr_Status);
    bool result = false;
    if (m_model.get(GRB_IntAttr_SolCount) != 0) {
        p_update_relaxed_solution(/*primal_only=*/true);
        if (p_round_greedy()) {
            result = true;
        }
    }
    if (status == GRB_INTERRUPTED || status == GRB_TIME_LIMIT ||
        status == GRB_MEM_LIMIT || status == GRB_NODE_LIMIT ||
        status == GRB_NODE_LIMIT || status == GRB_ITERATION_LIMIT)
    {
        throw InterruptError();
    }
    if (status != GRB_OPTIMAL) {
        throw std::logic_error("Unexpected status after solving relaxation!");
    }
    return result;
}

IncrementalGurobiG2CliqueSolver::SwitchBack
IncrementalGurobiG2CliqueSolver::p_switch_to_mip() {
    for (GRBVar var : m_existing_clique_vars) {
        var.set(GRB_CharAttr_VType, GRB_BINARY);
        var.set(GRB_DoubleAttr_Start, 0.0);
    }
    for (FullVertexIndex fi : m_best_clique) {
        m_all_clique_vars[fi].set(GRB_DoubleAttr_Start, 1.0);
    }
    m_model.set(GRB_IntParam_LazyConstraints, 1);
    m_model.setCallback(&m_lazy_nonedge_callback);
    return SwitchBack(this);
}

bool IncrementalGurobiG2CliqueSolver::p_prohibit_extend_explicit(
    OrderedSolutionIter vit, OrderedSolutionIter wit) {
    CliqueModelIndex ci1 = vit->clique_model_index;
    CliqueModelIndex ci2 = wit->clique_model_index;
    FullVertexIndex fi1 = m_full_vertex_indices[ci1];
    FullVertexIndex fi2 = m_full_vertex_indices[ci2];
    GRBVar var1 = m_existing_clique_vars[ci1];
    GRBVar var2 = m_existing_clique_vars[ci2];
    m_constraint_index_and_score_buffer.clear();
    for (std::size_t i = 0, ncons = m_constraints.size(); i < ncons; ++i) {
        if (!m_constraints[i].is_explicit())
            continue;
        auto& expl = m_constraints[i].explicit_set();
        if (expl.contains(fi1)) {
            if (expl.extend_if_possible(fi2, ci2)) {
                m_model.chgCoeff(m_constraints[i].constraint(), var2, 1.0);
                return true;
            }
        }
        if (expl.contains(fi2)) {
            if (expl.extend_if_possible(fi1, ci1)) {
                m_model.chgCoeff(m_constraints[i].constraint(), var1, 1.0);
                return true;
            }
        }
    }
    for (std::size_t i = 0, ncons = m_constraints.size(); i < ncons; ++i) {
        if (!m_constraints[i].is_explicit())
            continue;
        auto& expl = m_constraints[i].explicit_set();
        if (expl.extend_if_possible(fi1, ci1, fi2, ci2)) {
            m_model.chgCoeff(m_constraints[i].constraint(), var1, 1.0);
            m_model.chgCoeff(m_constraints[i].constraint(), var2, 1.0);
            return true;
        }
    }
    return false;
}

void IncrementalGurobiG2CliqueSolver::p_switch_to_lp() {
    for (GRBVar var : m_existing_clique_vars) {
        var.set(GRB_CharAttr_VType, GRB_CONTINUOUS);
    }
    m_model.set(GRB_IntParam_LazyConstraints, 0);
    m_model.setCallback(&m_abort_callback);
}

} // namespace sammy

#endif
