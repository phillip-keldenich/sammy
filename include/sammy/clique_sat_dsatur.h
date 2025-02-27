#ifndef SAMMY_CLIQUE_SAT_DSATUR_H_INCLUDED_
#define SAMMY_CLIQUE_SAT_DSATUR_H_INCLUDED_

#include "algorithm_ex.h"
#include "best_k.h"
#include "class_completer.h"
#include "dynamic_bitset.h"
#include "gurobi.h"
#include "literals.h"
#include "output.h"
#include "pair_infeasibility_map.h"
#include "partial_solution.h"
#include "thread_interrupt.h"
#include "vertex_operations.h"
#include <boost/iterator/transform_iterator.hpp>
#include <variant>

namespace sammy {

/**
 * The maximum number of new constraints added in
 * each call of p_prohibit_violated_nonedges.
 */
static constexpr std::size_t CAP_PROHIBIT_VIOLATED_NUM_NEW_CONSTRAINTS = 20;

/**
 * The number of cheap cut rounds to do before
 * performing a check to see whether the gap
 * has shrunk enough to continue cheap cuts.
 */
static constexpr std::size_t CHEAP_CUT_ROUNDS_PER_GAP_CHECK = 10;

/**
 * Gap: (objective - m_best_clique.size()) / m_best_clique.size()
 * The gap reduction required between two groups of
 * CHEAP_CUT_ROUNDS_PER_GAP_CHECK cheap cut rounds to
 * continue cheap cuts.
 */
static constexpr double CHEAP_CUT_GAP_REDUCTION_REQUIRED = 0.02;

/**
 * The maximum number of cheap cut rounds to do
 * before giving up on cheap cuts for the time being.
 */
static constexpr std::size_t CHEAP_CUT_ROUNDS_HARD_CAP = 100;

/**
 * The number of greedy cuts to add in each round of
 * greedy cut generation.
 */
static constexpr std::size_t GREEDY_CUT_CANDIDATE_SET_SIZE = 7;

template <typename IncrementalSATSolver> class CliqueSatDSaturSolver {
  public:
    /**
     * @brief Constructs a new CliqueSatDSaturSolver object.
     * @param considered_vertices All vertices to consider.
     * @param infeasibility_map The pair infeasibility map.
     *                          Is not changed. learn_infeasibilities should be
     * called before.
     * @param clause_db The clause database. May be changed by learning conflict
     * clauses.
     * @param best_local_mes The best MES involving vertices from
     * considered_vertices.
     * @param best_global_mes The best MES involving all vertices.
     * @param best_global_lb The best lower bound on the number of configuration
     * needed.
     * @param covering_assignments A set of complete configurations covering all
     * vertices in considered_vertices.
     */
    CliqueSatDSaturSolver(std::vector<Vertex> considered_vertices,
                          PairInfeasibilityMap* infeasibility_map,
                          ClauseDB& clause_db,
                          const std::vector<Vertex>& best_local_mes,
                          std::size_t best_global_mes,
                          std::size_t best_global_lb,
                          std::vector<DynamicBitset> covering_assignments,
                          bool enable_incremental_clique = true);

    /**
     * Abort the solve (on timeout, optimality, LNS with better solution, ...).
     */
    void abort();

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
    SolveResult solve();

    /**
     * Get the best solution found. Only valid if
     * the SolveResult was IMPROVED_SOLUTION.
     *
     * @return The best solution found.
     */
    PartialSolution get_partial_solution() const;

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
     * on the global MES size.
     */
    bool improved_mes() const noexcept {
        return m_best_clique.size() > m_initial_global_clique_size;
    }

    /**
     * Get a copy of the best clique found.
     */
    std::vector<Vertex> get_best_mes() const {
        std::vector<Vertex> result;
        result.reserve(m_best_clique.size());
        std::transform(m_best_clique.begin(), m_best_clique.end(),
                       std::back_inserter(result),
                       [&](std::size_t vi) { return m_all_vertices[vi]; });
        return result;
    }

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
        std::vector<Vertex> result;
        result.reserve(m_lower_bound_subgraph.size());
        std::transform(m_lower_bound_subgraph.begin(),
                       m_lower_bound_subgraph.end(), std::back_inserter(result),
                       [&](std::size_t vi) { return m_all_vertices[vi]; });
        return result;
    }

    /**
     * Set a callback to be invoked each time we try to search for a new clique.
     * It can return a clique (containing any vertices); the solver will filter
     * out all irrelevant vertices and only consider the ones in m_all_vertices.
     */
    void set_clique_candidate_callback(
        std::function<std::vector<Vertex>()> callback) {
        m_clique_candidate_callback = std::move(callback);
    }

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
     * The literal type of the SAT solver.
     */
    using SatLit = typename IncrementalSATSolver::Lit;

    /**
     * A SAT literal or a fixed boolean value.
     */
    using LitOrFixed = std::variant<SatLit, bool>;

    /**
     * The pair infeasibility map.
     */
    PairInfeasibilityMap* m_inf_map;

    /**
     * The number of features.
     */
    Index m_n_all;

    /**
     * The number of concrete features.
     */
    Index m_n_concrete;

    /**
     * A propagator we can reuse.
     * After use, it must be emptied again.
     */
    SharedDBPropagator m_empty_propagator;

    /**
     * Optional callback that is called to generate clique
     * candidates (may contain vertices not in m_all_vertices,
     * which are then filtered out).
     */
    std::function<std::vector<Vertex>()> m_clique_candidate_callback;

    /**
     * Optional callback that is called each time we find
     * a new lower bound.
     */
    std::function<void(std::size_t, const std::vector<Vertex>&)>
        m_lower_bound_callback;

    /**
     * The set of all vertices we need to consider.
     * These can be reduced by removing implied vertices.
     */
    std::vector<Vertex> m_all_vertices;

    /**
     * For each concrete literal, the list of vertices containing that literal.
     */
    std::vector<std::vector<std::size_t>>
        m_vertices_containing_concrete_literal;

    /**
     * For each concrete literal, the list of vertex indices
     * into m_existing_clique_vars containing that literal
     * (that are present in the clique).
     */
    std::vector<std::vector<std::size_t>>
        m_clique_model_vertices_containing_concrete_literal;

    /**
     * For each (concrete or non-concrete) literal, the set of vertices
     * _implying_ that literal by propagation.
     */
    std::vector<DynamicBitset> m_vertices_implying_literal;

    /**
     * The adjacency matrix of the graph.
     * A 'true' means that there is an edge (i.e., conflict)
     * between two vertices (i.e., interactions).
     * However, a 'false' only means the edge is not there
     * if the corresponding entry in m_definitive_nonedges is true.
     * Initialized from conflicts according to m_vertices_implying_literal,
     * i.e., initially corresponding to graph G_1.
     */
    std::vector<DynamicBitset> m_adjacency_matrix;

    /**
     * A 'true' means that there is definitely no edge,
     * i.e., propagation of interactions m_all_vertices[i]
     * and m_all_vertices[j] at the same time found no conflict.
     */
    std::vector<DynamicBitset> m_definitive_nonedges;

    /**
     * Whether the incremental clique solver should be used.
     * If false (not the default), the solver will not be
     * run and the initial clique is used throughout the entire solve.
     */
    bool m_incremental_clique_turned_off;

    /**
     * The Gurobi clique solver model.
     */
    GRBModel m_clique_model;

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
     * Buffer for variable values from the last solution
     * of the clique solver model.
     */
    std::unique_ptr<double[]> m_existing_clique_var_values;
    std::size_t m_existing_clique_var_value_count = 0;
    std::string m_last_vertex_addition_source = "not added any";

    /**
     * The last relaxation objective value.
     */
    double m_last_clique_objective = -std::numeric_limits<double>::infinity();

    /**
     * The last relaxation objective value rounded to
     * an integer upper bound (up to some numeric slack).
     */
    std::size_t m_last_rounded_objective = 0;

    /**
     * Struct containing clique model index and value
     * for sorting vertices by their value in the last solution.
     */
    struct OrderedSolutionValue {
        std::size_t clique_model_index;
        double value;

        bool operator<(const OrderedSolutionValue& other) const noexcept {
            return value > other.value;
        }
    };

    /**
     * The non-zero elements of the last solution of the clique solver model,
     * sorted by the variable value in the solution.
     */
    std::vector<OrderedSolutionValue> m_clq_ordered_solution;

    using OrderedSolutionIter =
        typename std::vector<OrderedSolutionValue>::const_iterator;

    /**
     * Collection/buffer for violated non-edges
     * in the relaxation of the clique model.
     */
    std::vector<std::pair<OrderedSolutionIter, OrderedSolutionIter>>
        m_clq_violated_nonedges;

    /**
     * Index i has the vertex index of the vertex
     * referred to by m_existing_clique_vars[i].
     */
    std::vector<std::size_t> m_existing_clique_var_indices;

    /**
     * The set of all vertices present in the clique solver model.
     * Indices into m_all_vertices.
     */
    DynamicBitset m_present_in_clique_model;

    /**
     * The vertices that are currently colored.
     */
    std::vector<std::size_t> m_currently_colored;

    /**
     * The set of vertices that are in
     * m_currently_colored, i.e., are currently
     * explicitly colored.
     */
    DynamicBitset m_in_currently_colored;

    /**
     * The order of vertices in the SAT model.
     * I.e., m_sat_vertex_order[i] is the vertex index
     * in m_all_vertices corresponding to the i-th vertex
     * in the SAT model.
     * Reverse mapping via m_vertex_info[vi].sat_index.
     * Over time, this set only grows unless the model
     * is reset on finding a new clique.
     */
    std::vector<std::size_t> m_sat_vertex_order;

    /**
     * Information on all vertices for the DSatur-like
     * configuration assignment algorithm.
     */
    struct VertexInfo {
        DynamicBitset
            open_classes; //< a superset of the open classes for this vertex
        std::size_t num_open_classes; //< number of set bits in open_classes
        bool in_some_class;    //< number of classes containing this vertex
        std::size_t degree;    //< degree of the vertex in the graph
        std::size_t sat_index; //< index of the vertex in the SAT solver

        VertexInfo(std::size_t num_initial_classes, std::size_t ub_num_classes,
                   std::size_t degree)
            : num_open_classes{num_initial_classes}, in_some_class{false},
              degree{degree},
              sat_index{std::numeric_limits<std::size_t>::max()} {
            open_classes.reserve(ub_num_classes);
            open_classes.assign(num_initial_classes, true);
        }

        bool close_class(std::size_t class_index) {
            if (open_classes[class_index]) {
                open_classes[class_index].reset();
                --num_open_classes;
                return true;
            }
            return false;
        }
    };

    /**
     * Compare vertex indices by the VertexInfo.
     * Used for finding the best K
     * vertices according to their info.
     */
    struct CompareVertexIndexByInfo {
        inline bool operator()(std::size_t v1, std::size_t v2) const noexcept;
        CliqueSatDSaturSolver* that;
    };

    /**
     * The vertex info for each vertex.
     */
    std::vector<VertexInfo> m_vertex_info;

    /**
     * The best K vertex indices according to their info.
     */
    BestK<std::size_t, CompareVertexIndexByInfo> m_best_infos;

    /**
     * Information on a single color class (partial configuration).
     */
    struct ColorClass {
        SharedDBPropagator propagator;
        std::size_t index;

        ColorClass(const SharedDBPropagator& empty, std::size_t index)
            : propagator(empty), index(index) {}

        void add_vertex(CliqueSatDSaturSolver* that, std::size_t vertex);
        bool try_add_vertex(CliqueSatDSaturSolver* that, std::size_t vertex);

      private:
        void p_added_vertex(CliqueSatDSaturSolver* that,
                            std::size_t old_trail_end);
    };

    /**
     * The color classes.
     */
    std::vector<ColorClass> m_color_classes;

    /**
     * If we found any, the set of vertices
     * with zero possible configurations.
     */
    std::vector<std::size_t> m_zero_candidates;

    /**
     * If we found any, the set of vertices with only
     * one possible configuration. Last priority
     * before returning to the DSatur-like choice.
     */
    std::vector<std::size_t> m_one_candidates;

    /**
     * Gurobi callback that only checks for aborts.
     */
    class GurobiCallbackCheckAbort : public GRBCallback {
      public:
        explicit GurobiCallbackCheckAbort(CliqueSatDSaturSolver* that) noexcept
            : that(that) {}

        void callback() override;

      private:
        CliqueSatDSaturSolver* that;
    };
    std::unique_ptr<GurobiCallbackCheckAbort> m_check_abort_cb;

    /**
     * The incremental SAT solver.
     */
    IncrementalSATSolver m_incremental_solver;

    /**
     * SAT model variables.
     * If present, m_class_literals[c][i]
     * is the i-th variable of class c.
     */
    std::vector<std::vector<LitOrFixed>> m_class_literals;

    /**
     * Local recorder of events.
     */
    EventRecorder* m_local_recorder{nullptr};

    /**
     * SAT model variables.
     * If present, m_vertex_in_class[c][si]
     * is the variable indicating that the vertex
     * with vertex index m_sat_vertex_order[si] is in class c.
     */
    std::vector<std::vector<LitOrFixed>> m_vertex_in_class;

    /**
     * SAT model variable.
     * A variable used to disable the constraint that
     * a vertex has to be in one of classes 0, ... c - 1.
     * Fixed as soon as c + 1 colors are available.
     */
    std::optional<SatLit> m_not_enough_colors;

    /**
     * The number of colors m_not_enough_colors is for.
     */
    std::size_t m_not_enough_colors_is_for = 0;

    /**
     * Buffer for building SAT clauses.
     */
    std::vector<SatLit> m_clause_buffer;

    /**
     * Buffer for building vertex sets.
     */
    std::vector<std::size_t> m_vertex_set_buffer;

    /**
     * Buffer for building vertex bitsets.
     */
    DynamicBitset m_vertex_bitset_buffer;

    /**
     * Buffer for GRBLinExprs.
     */
    GRBLinExpr m_buffer_expr;

    /**
     * Buffer for GRBColumns.
     */
    GRBColumn m_column_buffer;

    /**
     * Buffer for building vertex index sets.
     */
    std::vector<std::size_t> m_vertices_buffer;

    /**
     * The best clique found on the vertex set.
     * Indices into m_all_vertices.
     */
    std::vector<std::size_t> m_best_clique;

    /**
     * The clique size last used by SAT.
     * If the clique size increases, we have to
     * rebuild the SAT model.
     */
    std::size_t m_last_sat_clique_size;

    /**
     * If the clique model contains vertices that are
     * not yet colored, there may be vertices in the
     * best clique found so far that are not colored yet.
     * This is the set of those vertices.
     *
     * These vertices are added to the SAT calls (to be used
     * to break symmetries).
     * An UNSAT result leaves these vertices uncolored.
     * At any time, these vertices have the highest priority,
     * disregarding the DSatur-like choice.
     */
    std::vector<std::size_t> m_in_clique_uncolored;

    /**
     * Vertices with zero leftover choices.
     */
    std::vector<std::size_t> m_zero_choice_vertices;

    /**
     * The initial size of the clique on m_all_vertices.
     */
    std::size_t m_initial_clique_size;

    /**
     * The size of the best-known clique (overall)
     * when we started the search.
     */
    std::size_t m_initial_global_clique_size;

    /**
     * The initial lower bound (overall) when
     * we started the search.
     */
    std::size_t m_initial_global_lower_bound;

    /**
     * The current lower bound.
     */
    std::size_t m_lower_bound;

    /**
     * The vertices inducing the subgraph on which
     * we established the lower bound.
     */
    std::vector<std::size_t> m_lower_bound_subgraph;

    /**
     * The upper bound (configurations covering all vertices
     * in m_all_vertices).
     */
    std::vector<DynamicBitset> m_covering_assignments;

    /**
     * Flag to abort the search.
     */
    std::atomic<bool> m_aborted{false};

    /**
     * Constraints in the Gurobi clique model:
     * complete assignments give us independent
     * sets that we use as constraints.
     */
    std::vector<DynamicBitset> m_clq_constr_complete_assignments;

    /**
     * Constraints in the Gurobi clique model
     * from the assignments in m_clq_constr_complete_assignments.
     */
    std::vector<GRBConstr> m_clq_constr_complete_constraints;

    /**
     * Dual values of the constraints in the Gurobi clique model
     * from the assignments in m_clq_constr_complete_assignments.
     */
    std::unique_ptr<double[]> m_clq_constr_complete_dual_values;

    /**
     * Constraints in the Gurobi clique model:
     * partial assignments give us independent
     * sets that we use as constraints.
     */
    std::vector<SharedDBPropagator> m_clq_constr_partial_assignments;

    /**
     * Constraints in the Gurobi clique model
     * from the assignments in m_clq_constr_partial_assignments.
     */
    std::vector<GRBConstr> m_clq_constr_partial_constraints;

    /**
     * Dual values of the constraints in the Gurobi clique model
     * from the assignments in m_clq_constr_partial_assignments.
     */
    std::unique_ptr<double[]> m_clq_constr_partial_dual_values;

    /**
     * Buffer for tracking the set of violated non-edges that we already
     * covered.
     */
    DynamicBitset m_covered_violated_nonedges;

    /**
     * Candidate for a greedily generated cutting plane.
     */
    struct GreedyCutCandidate {
        std::vector<std::size_t> vertex_indices;
        double cut_value;

        // is-better-than operator
        bool operator<(const GreedyCutCandidate& other) const noexcept {
            return cut_value > other.cut_value;
        }

        // is-better-than-score operator
        bool operator<(double cut_value) const noexcept {
            return this->cut_value > cut_value;
        }

        template <typename... VArgs>
        GreedyCutCandidate(double cut_value, VArgs&&... args)
            : vertex_indices(std::forward<VArgs>(args)...),
              cut_value(cut_value) {}
    };

    /**
     * Buffer for the best greedy cut candidates during greedy cut generation.
     */
    BestK<GreedyCutCandidate> m_greedy_cut_candidates;

    /**
     * PricingEntry: vertex and dual price + 1.0.
     */
    struct PricingEntry {
        PricingEntry(std::size_t v, double dw) noexcept
            : vertex_index(v), dual_weight(dw) {}
        std::size_t vertex_index;
        double dual_weight;
    };

    /**
     * The subset of 'good' vertices identified during pricing.
     * Good vertices are vertices for which the dual constraint
     * is violated by a margin of at least 0.01, i.e.,
     * the corresponding dual constraint's LHS is at most 0.99.
     */
    std::vector<PricingEntry> m_pricing_good_vertices;

    /**
     * The subset of 'possible' vertices identified during pricing.
     * Possible vertices are vertices for which the dual constraint
     * is close to being violated, i.e., the dual constraint's LHS
     * is between 0.99 and 1.01.
     */
    std::vector<PricingEntry> m_pricing_possible_vertices;

    // ---------- PRIVATE METHODS ----------
    /**
     * Call the lower bound callback with the current LB.
     */
    void p_call_lb_callback();

    /**
     * Initialize the adjacency matrices.
     * Called by the constructor.
     */
    void p_initialize_graph();

    /**
     * Initialize m_vertices_containing_concrete_literal.
     */
    void p_initialize_vertices_containing_concrete_literal();

    /**
     * Initialize m_vertices_implying_literal.
     */
    void p_initialize_vertices_implying_literal();

    /**
     * Initialize m_adjacency_matrix from
     * m_vertices_impliying_literal.
     */
    void p_initialize_matrix_from_implied_literals();

    /**
     * Transform a set of vertices to a set of indices.
     */
    std::vector<std::size_t> p_clique_to_indices(const std::vector<Vertex>& cvs,
                                                 bool must_match = true);

    /**
     * Extract the vertex indices into m_all_vertices that are
     * covered the given assignment into the vector vertices.
     */
    void
    p_extract_vertices_from_configuration(const DynamicBitset& assignment,
                                          std::vector<std::size_t>& vertices);

    /**
     * Extract the vertex indices into m_all_vertices that are
     * covered by the given assignment into m_vertices_buffer,
     * clearing m_vertices_buffer first.
     */
    void p_buffer_vertices_from_configuration(const DynamicBitset& assignment) {
        m_vertices_buffer.clear();
        p_extract_vertices_from_configuration(assignment, m_vertices_buffer);
    }

    enum class VertexChoiceType {
        ALL_COVERED,       //< all vertices are colored!
        FROM_CLIQUE,       //< from m_in_clique_uncolored
        WITH_ZERO_CHOICES, //< an element with zero choices
        WITH_ONE_CHOICE,   //< an element with one choice
        REGULAR_CHOICE     //< a regular DSatur-like choice
    };

    /**
     * A struct describing the next vertex to color.
     */
    struct VertexChoice {
        /**
         * The vertex to color.
         */
        std::size_t vertex;

        /**
         * Characterises the situation/choice.
         */
        VertexChoiceType type;
    };

    /**
     * Find the next vertex to color.
     */
    VertexChoice p_choose_next_vertex();

    /**
     * Check the clique for vertices to color.
     * Keeps the invariant that any vertex in m_best_clique
     * is either in m_currently_colored or in m_in_clique_uncolored.
     *
     * @param vc the vertex choice to fill in.
     * @return true if a vertex was found, false otherwise.
     */
    bool p_choose_next_vertex_color_clique(VertexChoice& vc);

    /**
     * Check m_one_candidates for vertices to color.
     *
     * @param vc the vertex choice to fill in.
     * @return true if a vertex was found, false otherwise.
     */
    bool p_choose_next_vertex_one_color_left(VertexChoice& vc);

    /**
     * Check all vertices for the best candidates
     * according to m_vertex_info.
     *
     * @return true if any uncolored vertex was found, false otherwise.
     */
    bool p_identify_best_candidates();

    /**
     * Check the best-looking candidates from m_best_infos
     * for their precise number of open classes.
     * @return the best vertex index according to our priority.
     */
    std::size_t p_precise_best_candidates();

    /**
     * Initialize the vertex info (m_vertex_info).
     */
    void p_initialize_vertex_info();

    /**
     * Initialize the color classes (m_color_classes).
     * Assigns the clique vertices to a color class each.
     */
    void p_initialize_color_classes();

    /**
     * Complete all color classes to full assignments.
     * @return true if all classes could be extended
     *         false if not all classes could be extended.
     *
     * If not all classes could be extended, the
     * SAT solver is called to find a complete assignment.
     * This may result in the lower bound increasing or
     * in a complete assignment that however does not
     * incorporate all vertices that were implicitly covered before.
     */
    bool p_complete_all_classes();

    /**
     * Handle the situation where the given vertex
     * has no possible choices.
     * @return true if continuing is sensible (in this case, the given vertex
     * will be added to m_currently_colored), false if we have proved that no
     * improvement is possible or were aborted (check m_aborted after return).
     */
    bool p_handle_zero_choice(std::size_t vertex);

    /**
     * Handle a given vertex.
     * Color it (and add it to m_currently_colored),
     * if that is possible.
     *
     * @return true if the vertex could be colored
     *         false if it turned out to have zero choices.
     */
    bool p_handle_vertex(std::size_t vertex);

    /**
     * Remove the given vertex from the stack
     * of vertices to be covered next, or the
     * stack of vertices with one candidate.
     */
    void p_colored_vertex(std::size_t vertex);

    /**
     * When we have a new clique, we completely
     * reset the SAT model.
     */
    void p_reset_sat_model();

    /**
     * Add a new color class, initially containing
     * the given vertex.
     */
    void p_add_color_class(std::size_t vertex);

    /**
     * When we need to run a new SAT/Clique model solve,
     * we use this method to compute the current set.
     */
    std::vector<std::size_t>
    p_compute_current_vertex_set(std::size_t new_vertex);

    /**
     * Compute/update the set of vertices that are
     * in the clique but not in m_currently_colored.
     */
    void p_compute_in_clique_uncolored();

    /**
     * Extend the given sorted unique set by the given elements;
     * the result is sorted & unique again.
     */
    void p_extend_sorted_unique_set(std::vector<std::size_t>& sorted_unique,
                                    const std::vector<std::size_t>& to_add);

    /**
     * Rebuild the SAT model.
     */
    void p_rebuild_sat_model(const std::vector<std::size_t>& current_vertices);

    /**
     * Rebuild a color class of the SAT model.
     */
    void p_rebuild_sat_color_class(std::optional<std::size_t> clique_vertex);

    /**
     * Append a vector of variables to the sat model to
     * obtain  the values of the features in a new color class.
     */
    std::vector<LitOrFixed>& p_append_unconstrained_class_lits();

    /**
     * Append a vector of variables to the sat model
     * to obtain the values of the features in a new color class,
     * initialized with the given vertex.
     */
    void p_append_constrained_class_lits(std::size_t clique_vertex);

    /**
     * Append a vector of variables to the sat model
     * to obtain the presence of vertices in a new color class.
     */
    std::vector<LitOrFixed>& p_append_vertex_in_class();

    /**
     * Lookup a sat literal or fixed value in a color class.
     */
    LitOrFixed p_lookup_sat_literal_in_class(Lit l, std::size_t class_index);

    /**
     * Negate a LitOrFixed.
     */
    static LitOrFixed p_negate(LitOrFixed l);

    /**
     * Return a LitOrFixed corresponding
     * to l1 & l2 (may add a new variable).
     */
    LitOrFixed p_and(LitOrFixed l1, LitOrFixed l2);

    /**
     * Add clauses to enforce that the given
     * class variables corresponds to a satisfying
     * assignment/valid configuration.
     */
    void p_enforce_formula_on(std::size_t class_index);

    /**
     * Add a clause for the given class variables.
     */
    void p_add_clause(const Lit* begin, const Lit* end,
                      std::size_t class_index);

    /**
     * Run SAT model.
     * @return true on SAT, false on UNSAT, nullopt on timeout/abort.
     */
    std::optional<bool> p_run_sat_model();

    /**
     * Extend the SAT model to all vertices in
     * current_vertices that are not yet part of the SAT model.
     */
    void p_extend_sat_model(const std::vector<std::size_t>& current_vertices);

    /**
     * Extend the SAT model to m_lower_bound color classes.
     */
    void p_extend_sat_model_to_classes();

    /**
     * Extract a satisfying assignment and
     * turn it into a partial coloring.
     * This may increase the number of members
     * in m_color_class to m_lower_bound.
     */
    template <typename ModelMapType>
    void p_extract_sat_model(const ModelMapType& model);

    /**
     * Add the constraint that the given vertex
     * must be in some color class.
     */
    void p_add_vertex_color_constraint(std::size_t sat_index,
                                       SatLit not_enough);

    /**
     * Check that m_not_enough_colors (and the vertex color constraints)
     * are set up for all vertices in m_sat_vertex_order.
     */
    void p_check_not_enough_colors();

    /**
     * Setup the clique model if necessary,
     * and solve it to hopefully find a larger clique.
     * Update m_lower_bound and m_best_clique if a better clique is found.
     *
     * @param current_vertices the vertices that are under consideration for the
     * SAT model
     * @return true if a larger clique was found, false otherwise
     */
    bool
    p_zero_choice_run_clique(const std::vector<std::size_t>& current_vertices);

    /**
     * Run the SAT model, adding colors until it is satisfiable.
     * @return true if satisfiable before hitting upper bound, false if
     * unsatisfiable or if we were aborted/interrupted.
     */
    bool p_run_until_satisfiable();

    /**
     * Reset the color classes to an empty state,
     * making sure that num_classes exist.
     * This is called when we replace the current coloring
     * with one from the SAT solver.
     */
    void p_reset_color_classes(std::size_t num_classes);

    /**
     * Look up a LitOrFixed value in a satisfying assignment of the SAT model.
     */
    template <typename ModelMap>
    bool p_lookup_in_model(const ModelMap& map, LitOrFixed l) const;

    /**
     * Completely reset the coloring in the rare case we failed to
     * extend all color classes to valid configurations.
     */
    void p_reset();

    /**
     * Create a LitOrFixed from a SatLit.
     */
    static inline LitOrFixed p_lit(SatLit sl) noexcept {
        return LitOrFixed(std::in_place_type<SatLit>, sl);
    }

    /**
     * Create a LitOrFixed from a fixed boolean value.
     */
    static inline LitOrFixed p_fixed(bool b) noexcept {
        return LitOrFixed(std::in_place_type<bool>, b);
    }

    /**
     * Called to initialize the clique model
     * on the first attempt to improve our clique.
     */
    bool
    p_first_run_clique_model(const std::vector<std::size_t>& current_vertices);

    /**
     * Called to add a complete assignment as
     * a <=1-constraint to the clique model.
     * The constraint RHS contains all vertices
     * that are present in the clique model and
     * covered by the complete assignment.
     * When new vertices are added to the model, they
     * are automatically added to all complete
     * assignment constraints that cover them.
     */
    void
    p_clq_add_complete_constraint(const DynamicBitset& complete_assignment);

    /**
     * Called to add a partial assignment as
     * a <=1-constraint to the clique model.
     * The constraint RHS contains all vertices
     * that are present in the clique model and
     * covered by the partial assignment.
     * When new vertices are added to the model, they
     * are automatically added to all partial
     * assignment constraints that cover them.
     * When the partial assignment is extended,
     * vertices present in the clique model
     * that are newly covered are added to the
     * constraint.
     */
    template <typename Propagator>
    void p_clq_add_partial_constraint(Propagator&& partial_assignment,
                                      bool expr_in_buffer = false);

    /**
     * Called to add a vertex as 0-1-variable to the clique model.
     */
    void p_clq_add_vertex(std::size_t vertex);

    /**
     * Notify that the partial assignment indicated by partial_index
     * has been extended (previous to the extension, the trail should
     * have had length old_trail_size).
     */
    void p_clq_extended_partial(std::size_t partial_index,
                                std::size_t old_trail_size);

    /**
     * Extend a partial assignment by adding a vertex.
     */
    bool p_clq_extend_partial_by_vertex(std::size_t partial_index,
                                        std::size_t vertex_index);

    /**
     * Setup model parameters.
     */
    void p_setup_clique_model();

    /**
     * Solve the current clique model relaxation.
     * @return true if the relaxation was solved, false on abort.
     */
    bool p_solve_relaxation();

    /**
     * Identify non-edges that violated the <=1-constraint
     * in the current clique model relaxation solution.
     * Actually checks for at least a violation by 0.01.
     * The violated nonedges are placed into m_clq_violated_nonedges.
     *
     * @return true if some violated non-edges were found, false otherwise.
     */
    bool p_identify_violated_nonedges();

    /**
     * Solve the 'full' relaxation, i.e.,
     * the relaxation including all lazily generated
     * non-edge constraints.
     * @return true if the full relaxation was solve, false on abort.
     */
    bool p_solve_full_relaxation();

    /**
     * Check if the indicated edge is actually a non-edge.
     * @return true if the edge is a non-edge, false otherwise.
     */
    bool p_is_nonedge(std::size_t vertex1, std::size_t vertex2);

    /**
     * Check if the indicated edge is actually an edge.
     * @return true if the edge is present in G_2, false otherwise.
     */
    bool p_is_edge(std::size_t vertex1, std::size_t vertex2) {
        return !p_is_nonedge(vertex1, vertex2);
    }

    /**
     * We lazily prohibit violated non-edges in the clique model.
     * The violated nonedges are stored in m_clq_violated_nonedges.
     * Each call to p_prohibit_violated_nonedges generates at most
     * one new constraint per vertex incident to a violated non-edge;
     * in practice, much fewer constraints are generated.
     * The method prefers strengthening existing partial constraints to
     * cover new violated non-edges over generating new constraints.
     */
    void p_prohibit_violated_nonedges();

    /**
     * Greedily extend the current partial assignment
     * in m_empty_propagator by adding non-edges from
     * m_clq_violated_nonedges[begin_index, end_index).
     */
    void p_extend_nonedge_cover_greedy(std::size_t begin_index,
                                       std::size_t end_index);

    /**
     * Try to handle the given violated non-edge by
     * extending an existing partial assignment constraint.
     */
    bool p_prohibit_extend_existing(OrderedSolutionIter vit,
                                    OrderedSolutionIter wit);

    /**
     * Use the clique model to try and improve the best known clique.
     * @return false if the search was aborted or the best clique was not
     * improved, true otherwise.
     */
    bool p_use_clique_model(const std::vector<std::size_t>& current_vertices);

    /**
     * Greedily round the current clique relaxation
     * to try and find an improved clique.
     */
    bool p_round_greedy();

    /**
     * Randomly turn the clique in m_vertex_set_buffer
     * with compatible vertices marked in m_vertex_bitset_buffer
     * into a maximal clique.
     */
    void p_randomly_extend_clique();

    /**
     * Ensure that the best clique found so far
     * is contained in the clique model.
     */
    void p_ensure_best_clique_in_model();

    /**
     * Try to add cuts to the clique model
     * coming from the current color classes.
     * @return true if cuts were added, false otherwise.
     */
    bool p_clique_model_cclass_cuts();

    /**
     * Run cheap cutting plane rounds.
     * @return false on abortion, true otherwise.
     */
    bool p_cheap_cut_rounds();

    /**
     * Run a single cheap cut round.
     * @return false on abortion or failure, true otherwise.
     */
    bool p_cheap_cut_round();

    /**
     * Run expensive cutting plane rounds.
     * @return false on abortion, true otherwise.
     */
    bool p_expensive_cut_rounds() {
        // TODO: implement expensive cuts (might not be very successful)
        return true;
    }

    /**
     * The gap in the current clique model.
     */
    double p_clique_gap() const noexcept {
        return (m_last_clique_objective - m_best_clique.size()) /
               m_best_clique.size();
    }

    /**
     * Greedily add vertices from the current relaxation
     * solution to existing partial-assignment cutting planes.
     * @return true if a constraint became violated, false otherwise.
     */
    bool p_greedy_add_to_cuts();

    /**
     * Greedily generate cutting planes from the current relaxation solution
     * and add them to the clique model.
     * @return true if a violated constraint was added, false otherwise.
     */
    bool p_greedy_generate_cuts();

    /**
     * Extend a greedy cut being built (partial assignment in
     * m_empty_propagator, vertices in m_vertex_set_buffer) by adding vertices
     * from the current relaxation solution.
     * @return the cut value of the extended cut.
     */
    double p_extend_greedy_cut(OrderedSolutionIter begin,
                               OrderedSolutionIter end);

    /**
     * Do a limited pricing round (limited means that there is a limit
     * on the number of vertices that can be added to the clique model).
     */
    bool p_limited_price_vertices(const std::vector<std::size_t>& vertices);

    /**
     * Do an unlimited pricing round (that adds all vertices
     * that could be an improvement to the clique model).
     */
    bool p_unlimited_price_vertices(const std::vector<std::size_t>& vertices);

    /**
     * Part of the pricing routine: Go through the given set of vertices,
     * and add good and possible vertices to the vectors m_pricing_good_vertices
     * and m_pricing_possible_vertices.
     * @return true if vertices were added, and false if there either were
     *         no possible vertices or if we were aborted.
     */
    bool
    p_price_collect_good_and_possible(const std::vector<std::size_t>& vertices);

    /**
     * For a vertex that is not present in the clique model,
     * compute the value of the LHS of its corresponding dual constraint.
     */
    double p_weight_in_dual(Vertex v);

    /**
     * Select a number of vertices from the set of
     * good vertices (or possible vertices if no good vertices exist).
     * The vertices are collected in m_vertex_set_buffer.
     */
    void p_pricing_select_with_goal(std::size_t goal_vertices);

    /**
     * Add the vertices in m_vertex_set_buffer to the clique model.
     */
    void p_pricing_add_vertices();

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
};

template <typename IncrementalSATSolver>
CliqueSatDSaturSolver<IncrementalSATSolver>::CliqueSatDSaturSolver(
    std::vector<Vertex> considered_vertices,
    PairInfeasibilityMap* infeasibility_map, ClauseDB& clause_db,
    const std::vector<Vertex>& best_local_mes, std::size_t best_global_mes,
    std::size_t best_global_lb, std::vector<DynamicBitset> covering_assignments,
    bool enable_incremental_clique)
    : m_inf_map(infeasibility_map), m_n_all(clause_db.num_vars()),
      m_n_concrete(infeasibility_map->get_n_concrete()),
      m_empty_propagator(&clause_db),
      m_all_vertices(std::move(considered_vertices)),
      m_vertices_containing_concrete_literal(2 * m_n_concrete),
      m_clique_model_vertices_containing_concrete_literal(2 * m_n_concrete),
      m_vertices_implying_literal(2 * m_n_all,
                                  DynamicBitset(m_all_vertices.size())),
      m_adjacency_matrix(m_all_vertices.size(),
                         DynamicBitset(m_all_vertices.size())),
      m_definitive_nonedges(m_all_vertices.size(),
                            DynamicBitset(m_all_vertices.size())),
      m_incremental_clique_turned_off(!enable_incremental_clique),
      m_clique_model(sammy::gurobi_environment(true)),
      m_all_clique_vars(m_all_vertices.size()),
      m_present_in_clique_model(m_all_vertices.size()),
      m_in_currently_colored(m_all_vertices.size(), false),
      m_best_infos(20, CompareVertexIndexByInfo{this}),
      m_best_clique(p_clique_to_indices(best_local_mes)),
      m_last_sat_clique_size(0), m_initial_clique_size(best_local_mes.size()),
      m_initial_global_clique_size(best_global_mes),
      m_initial_global_lower_bound(best_global_lb),
      m_lower_bound(m_best_clique.size()),
      m_lower_bound_subgraph(m_best_clique),
      m_covering_assignments(std::move(covering_assignments)),
      m_greedy_cut_candidates(GREEDY_CUT_CANDIDATE_SET_SIZE) {
    p_initialize_graph();
    p_initialize_vertex_info();
    p_initialize_color_classes();
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_initialize_graph() {
    p_initialize_vertices_containing_concrete_literal();
    p_initialize_vertices_implying_literal();
    // TODO: the next function can take way too long.
    p_initialize_matrix_from_implied_literals();
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_initialize_vertices_containing_concrete_literal() {
    assert(m_vertices_containing_concrete_literal.size() == 2 * m_n_concrete);
    for (std::size_t vi = 0, vn = m_all_vertices.size(); vi < vn; ++vi) {
        Vertex v = m_all_vertices[vi];
        assert(lit::var(v.first) < m_n_concrete);
        assert(lit::var(v.second) < m_n_concrete);
        m_vertices_containing_concrete_literal[v.first].push_back(vi);
        m_vertices_containing_concrete_literal[v.second].push_back(vi);
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_initialize_vertices_implying_literal() {
    assert(m_vertices_implying_literal.size() == 2 * m_n_all);
    for (std::size_t vi = 0, n = m_all_vertices.size(); vi < n; ++vi) {
        assert(m_empty_propagator.get_current_level() == 0);
        Vertex v = m_all_vertices[vi];
        if (push_vertex(m_empty_propagator, v) < 0)
            throw std::logic_error("Infeasible interaction in m_all_vertices!");
        for (Lit l : m_empty_propagator.get_trail()) {
            m_vertices_implying_literal[l][vi].set();
        }
        m_empty_propagator.reset_to_zero();
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_initialize_matrix_from_implied_literals() {
    assert(m_adjacency_matrix.size() == m_all_vertices.size());
    assert(m_adjacency_matrix.empty() ||
           m_adjacency_matrix[0].size() == m_all_vertices.size());
    for (std::size_t vi = 0, n = m_all_vertices.size(); vi < n; ++vi) {
        assert(m_empty_propagator.get_current_level() == 0);
        Vertex v = m_all_vertices[vi];
        DynamicBitset& row = m_adjacency_matrix[vi];
        push_vertex(m_empty_propagator, v);
        for (Lit lpos : m_empty_propagator.get_trail()) {
            row |= m_vertices_implying_literal[lit::negate(lpos)];
        }
        m_empty_propagator.reset_to_zero();
        if (vi % 16 == 15)
            throw_if_interrupted();
    }
}

template <typename IncrementalSATSolver>
std::vector<std::size_t>
CliqueSatDSaturSolver<IncrementalSATSolver>::p_clique_to_indices(
    const std::vector<Vertex>& cvs, bool must_match) {
    PairHashSet<Vertex> in_clique(cvs.begin(), cvs.end());
    std::vector<std::size_t> result;
    result.reserve(cvs.size());
    for (std::size_t vi = 0, n = m_all_vertices.size(); vi < n; ++vi) {
        Vertex v = m_all_vertices[vi];
        if (in_clique.count(v)) {
            result.push_back(vi);
        }
    }
    if (must_match && (result.size() != cvs.size())) {
        throw std::logic_error("Clique contains missing vertex!");
    }
    return result;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::
    p_extract_vertices_from_configuration(const DynamicBitset& assignment,
                                          std::vector<std::size_t>& vertices) {
    Lit block = NIL;
    for (std::size_t i = 0, n = m_all_vertices.size(); i < n; ++i) {
        Vertex v = m_all_vertices[i];
        if (v.first == block)
            continue;
        if (!lit::is_true_in(v.first, assignment)) {
            block = v.first;
            continue;
        }
        if (lit::is_true_in(v.second, assignment)) {
            vertices.push_back(i);
        }
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::abort() {
    m_aborted.store(true);
    m_incremental_solver.terminate();
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::CompareVertexIndexByInfo::
operator()(std::size_t v1, std::size_t v2) const noexcept {
    // return true if 'v1' has higher priority than 'v2'
    const VertexInfo& vinfo1 = that->m_vertex_info[v1];
    const VertexInfo& vinfo2 = that->m_vertex_info[v2];
    if (vinfo1.in_some_class != vinfo2.in_some_class) {
        return !vinfo1.in_some_class;
    }
    if (vinfo1.num_open_classes != vinfo2.num_open_classes) {
        return vinfo1.num_open_classes < vinfo2.num_open_classes;
    }
    return vinfo1.degree > vinfo2.degree;
}

template <typename IncrementalSATSolver>
typename CliqueSatDSaturSolver<IncrementalSATSolver>::VertexChoice
CliqueSatDSaturSolver<IncrementalSATSolver>::p_choose_next_vertex() {
    VertexChoice result;
    // first priority: color leftover vertices from the clique
    if (p_choose_next_vertex_color_clique(result))
        return result;
    // second priority: color the vertices that only have one color left
    if (p_choose_next_vertex_one_color_left(result))
        return result;
    // third priority: use the DSatur-style choice
    if (!p_identify_best_candidates()) {
        result.vertex = std::numeric_limits<std::size_t>::max();
        result.type = VertexChoiceType::ALL_COVERED;
        return result;
    }
    std::size_t v = p_precise_best_candidates();
    result.vertex = v;
    if (!m_zero_candidates.empty()) {
        result.type = VertexChoiceType::WITH_ZERO_CHOICES;
    } else if (m_vertex_info[v].num_open_classes == 1) {
        result.type = VertexChoiceType::WITH_ONE_CHOICE;
    } else {
        result.type = VertexChoiceType::REGULAR_CHOICE;
    }
    return result;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_choose_next_vertex_color_clique(VertexChoice& vc) {
    while (!m_in_clique_uncolored.empty()) {
        std::size_t v = m_in_clique_uncolored.back();
        if (m_vertex_info[v].in_some_class) {
            m_in_clique_uncolored.pop_back();
            m_currently_colored.push_back(v);
            m_in_currently_colored[v].set();
        } else {
            vc = VertexChoice{v, VertexChoiceType::FROM_CLIQUE};
            return true;
        }
    }
    return false;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::
    p_choose_next_vertex_one_color_left(VertexChoice& vc) {
    while (!m_one_candidates.empty()) {
        std::size_t v = m_one_candidates.back();
        m_one_candidates.pop_back();
        if (!m_vertex_info[v].in_some_class) {
            vc = VertexChoice{v, VertexChoiceType::WITH_ONE_CHOICE};
            return true;
        }
    }
    return false;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_identify_best_candidates() {
    m_zero_candidates.clear();
    m_one_candidates.clear();
    m_best_infos.clear();
    for (std::size_t i = 0, n = m_all_vertices.size(); i < n; ++i) {
        if (m_vertex_info[i].in_some_class)
            continue;
        if (m_vertex_info[i].num_open_classes == 0) {
            m_zero_candidates.push_back(i);
        } else if (m_vertex_info[i].num_open_classes == 1) {
            m_one_candidates.push_back(i);
        }
        m_best_infos.push(i);
    }
    return !m_best_infos.elements().empty();
}

template <typename IncrementalSATSolver>
std::size_t
CliqueSatDSaturSolver<IncrementalSATSolver>::p_precise_best_candidates() {
    for (std::size_t candidate : m_best_infos.elements()) {
        Vertex vc = m_all_vertices[candidate];
        auto& vi = m_vertex_info[candidate];
        std::size_t open_before = vi.num_open_classes;
        if (open_before == 0)
            continue;
        for (std::size_t candidate_class : vi.open_classes.ones()) {
            auto& cclass = m_color_classes[candidate_class];
            SharedDBPropagator& prop = cclass.propagator;
            if (!can_push(prop, vc)) {
                vi.num_open_classes -= 1;
                vi.open_classes[candidate_class].reset();
            }
        }
        if (open_before != vi.num_open_classes) {
            if (vi.num_open_classes == 0) {
                m_zero_candidates.push_back(candidate);
            } else if (vi.num_open_classes == 1) {
                m_one_candidates.push_back(candidate);
            }
        }
    }
    const auto& elms = m_best_infos.elements();
    std::size_t best_v = *(std::min_element(elms.begin(), elms.end(),
                                            CompareVertexIndexByInfo{this}));
    return best_v;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_initialize_vertex_info() {
    const auto n = m_all_vertices.size();
    m_vertex_info.reserve(n);
    std::size_t ub_num_classes = m_covering_assignments.size();
    std::size_t initial_num_classes = m_best_clique.size();
    for (std::size_t vi = 0; vi < n; ++vi) {
        const DynamicBitset& row = m_adjacency_matrix[vi];
        std::size_t degree = row.count();
        m_vertex_info.emplace_back(initial_num_classes, ub_num_classes, degree);
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_initialize_color_classes() {
    if (m_empty_propagator.get_current_level() != 0)
        throw std::logic_error("Non-empty propagator!");
    m_color_classes.reserve(m_covering_assignments.size());
    for (std::size_t i = 0; i < m_best_clique.size(); ++i) {
        m_color_classes.emplace_back(m_empty_propagator, i);
        m_color_classes.back().add_vertex(this, m_best_clique[i]);
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::ColorClass::add_vertex(
    CliqueSatDSaturSolver* that, std::size_t vertex) {
    Vertex v = that->m_all_vertices[vertex];
    const auto& trail = propagator.get_trail();
    std::size_t old_trail_length = trail.size();
    int p = push_vertex(propagator, v);
    if (p < 0) {
        throw std::logic_error("Called add_vertex on incompatible class!");
    }
    if (p > 0) {
        p_added_vertex(that, old_trail_length);
    }
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::ColorClass::try_add_vertex(
    CliqueSatDSaturSolver* that, std::size_t vertex) {
    Vertex v = that->m_all_vertices[vertex];
    const auto& trail = propagator.get_trail();
    std::size_t old_trail_length = trail.size();
    int p = push_vertex(propagator, v);
    if (p < 0) {
        return false;
    }
    if (p > 0) {
        p_added_vertex(that, old_trail_length);
    }
    return true;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::ColorClass::p_added_vertex(
    CliqueSatDSaturSolver* that, std::size_t old_trail_end) {
    // TODO: unnecessary if we do the thing below (but that may be expensive)
    // except for maybe 'extra edges' (from extended graph G_2)?
    //    for(std::size_t vother : that->m_adjacency_matrix[vertex].ones()) {
    //        that->m_vertex_info[vother].close_class(index);
    //    }
    const auto& trail = propagator.get_trail();
    Lit nclit = that->m_n_concrete * 2;
    for (Lit lnew : IteratorRange{trail.begin() + old_trail_end, trail.end()}) {
        Lit lneg = lit::negate(lnew);
        if (lnew < nclit) {
            for (std::size_t vpot :
                 that->m_vertices_containing_concrete_literal[lnew])
            {
                Vertex vo = that->m_all_vertices[vpot];
                if (propagator.is_true(vo.first) &&
                    propagator.is_true(vo.second))
                {
                    that->m_vertex_info[vpot].in_some_class = true;
                }
            }
        }
        // TODO the thing below:
        for (std::size_t vout : that->m_vertices_implying_literal[lneg].ones())
        {
            that->m_vertex_info[vout].close_class(index);
        }
    }
}

template <typename IncrementalSATSolver>
auto CliqueSatDSaturSolver<IncrementalSATSolver>::solve() -> SolveResult {
    for (;;) {
        if (m_aborted.load()) {
            return SolveResult::ABORTED;
        }
        VertexChoice vc = p_choose_next_vertex();
        switch (vc.type) {
        case VertexChoiceType::ALL_COVERED:
            if (!p_complete_all_classes()) {
                if (m_aborted.load())
                    return SolveResult::ABORTED;
                if (m_lower_bound == m_covering_assignments.size())
                    return SolveResult::SOLUTION_WAS_OPTIMAL;
                continue;
            }
            if (m_covering_assignments.size() <= m_color_classes.size())
                return SolveResult::SOLUTION_WAS_OPTIMAL;
            return SolveResult::IMPROVED_SOLUTION;

        case VertexChoiceType::WITH_ZERO_CHOICES: {
            p_report_event("BEGIN_HANDLE_ZERO_CHOICE",
                           {{"vertices_handled", m_currently_colored.size()}},
                           "vertices_handled");
            bool handled = p_handle_zero_choice(vc.vertex);
            p_report_event("DONE_HANDLE_ZERO_CHOICE",
                           {{"vertices_handled", m_currently_colored.size()}},
                           "vertices_handled");
            if (handled) {
                p_colored_vertex(vc.vertex);
                continue;
            }
            return m_aborted.load() ? SolveResult::ABORTED
                                    : SolveResult::SOLUTION_WAS_OPTIMAL;
        }

        default: {
            if (!p_handle_vertex(vc.vertex)) {
                p_report_event(
                    "BEGIN_HANDLE_ZERO_CHOICE",
                    {{"vertices_handled", m_currently_colored.size()}},
                    "vertices_handled");
                bool handled = p_handle_zero_choice(vc.vertex);
                p_report_event(
                    "DONE_HANDLE_ZERO_CHOICE",
                    {{"vertices_handled", m_currently_colored.size()}},
                    "vertices_handled");
                if (!handled) {
                    return m_aborted.load() ? SolveResult::ABORTED
                                            : SolveResult::SOLUTION_WAS_OPTIMAL;
                }
            }
            p_colored_vertex(vc.vertex);
            break;
        }
        }
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_colored_vertex(
    std::size_t vertex) {
    if (!m_in_clique_uncolored.empty() &&
        m_in_clique_uncolored.back() == vertex)
    {
        m_in_clique_uncolored.pop_back();
    }
    if (!m_one_candidates.empty() && m_one_candidates.back() == vertex) {
        m_one_candidates.pop_back();
    }
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_handle_vertex(
    std::size_t vertex) {
    auto& vinfo = m_vertex_info[vertex];
    if (vinfo.in_some_class) {
        m_currently_colored.push_back(vertex);
        m_in_currently_colored[vertex].set();
        return true;
    }
    for (std::size_t candidate : vinfo.open_classes.ones()) {
        auto& cclass = m_color_classes[candidate];
        if (cclass.try_add_vertex(this, vertex)) {
            m_currently_colored.push_back(vertex);
            m_in_currently_colored[vertex].set();
            return true;
        } else {
            vinfo.close_class(candidate);
        }
    }
    return false;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_handle_zero_choice(
    std::size_t vertex) {
    auto current_vertices = p_compute_current_vertex_set(vertex);
    bool rclique = false;
    if (!m_incremental_clique_turned_off) {
        p_report_event("BEGIN_ZERO_CHOICE_CLIQUE", {});
        rclique = p_zero_choice_run_clique(current_vertices);
        p_report_event("DONE_ZERO_CHOICE_CLIQUE",
                       {{"best_clique_size", m_best_clique.size()}},
                       "best_clique_size");
    }
    if (rclique) {
        // found a better clique (m_best_clique, m_lower_bound updated)!
        if (m_best_clique.size() > m_color_classes.size()) {
            if (m_best_clique.size() == m_covering_assignments.size()) {
                return false;
            }
            p_extend_sorted_unique_set(current_vertices, m_best_clique);
            p_reset_sat_model();
            p_rebuild_sat_model(current_vertices);
            return p_run_until_satisfiable();
        }
    }
    if (m_aborted.load())
        return false;
    p_report_event("BEGIN_ZERO_CHOICE_SAT_MODELING", {});
    if (m_best_clique.size() != m_last_sat_clique_size) {
        p_reset_sat_model();
        p_rebuild_sat_model(current_vertices);
    } else {
        p_extend_sat_model_to_classes();
        p_extend_sat_model(current_vertices);
    }
    p_report_event("DONE_ZERO_CHOICE_SAT_MODELING", {});
    p_report_event("BEGIN_ZERO_CHOICE_SAT_SOLVE", {});
    std::optional<bool> result = p_run_sat_model();
    p_report_event("DONE_ZERO_CHOICE_SAT_SOLVE", {});
    if (!result)
        return false;
    if (*result)
        return true;
    m_lower_bound += 1;
    m_lower_bound_subgraph = std::move(current_vertices);
    p_call_lb_callback();
    if (m_lower_bound == m_covering_assignments.size())
        return false;
    p_add_color_class(vertex);
    p_compute_in_clique_uncolored();
    return true;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_reset_sat_model() {
    if (m_class_literals.empty())
        return;
    m_incremental_solver.reset();
    m_class_literals.clear();
    m_vertex_in_class.clear();
    for (std::size_t vi : m_sat_vertex_order) {
        m_vertex_info[vi].sat_index = std::numeric_limits<std::size_t>::max();
    }
    m_sat_vertex_order.clear();
    m_not_enough_colors.reset();
    m_not_enough_colors_is_for = 0;
}

template <typename IncrementalSATSolver>
std::vector<std::size_t>
CliqueSatDSaturSolver<IncrementalSATSolver>::p_compute_current_vertex_set(
    std::size_t new_vertex) {
    std::vector<std::size_t> cover_vertex_set = m_in_clique_uncolored;
    cover_vertex_set.insert(cover_vertex_set.end(), m_zero_candidates.begin(),
                            m_zero_candidates.end());
    cover_vertex_set.push_back(new_vertex);
    cover_vertex_set.insert(cover_vertex_set.end(), m_currently_colored.begin(),
                            m_currently_colored.end());
    std::sort(cover_vertex_set.begin(), cover_vertex_set.end());
    cover_vertex_set.erase(
        std::unique(cover_vertex_set.begin(), cover_vertex_set.end()),
        cover_vertex_set.end());
    return cover_vertex_set;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_add_color_class(
    std::size_t vertex) {
    std::size_t new_class = m_color_classes.size();
    for (std::size_t vi = 0, n = m_all_vertices.size(); vi < n; ++vi) {
        m_vertex_info[vi].open_classes.push_back(true);
        m_vertex_info[vi].num_open_classes += 1;
    }
    m_color_classes.emplace_back(m_empty_propagator, new_class);
    m_color_classes.back().add_vertex(this, vertex);
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_compute_in_clique_uncolored() {
    m_in_clique_uncolored.clear();
    for (std::size_t v : m_best_clique) {
        if (!m_in_currently_colored[v]) {
            m_in_clique_uncolored.push_back(v);
        }
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_extend_sorted_unique_set(
    std::vector<std::size_t>& sorted_unique,
    const std::vector<std::size_t>& to_add) {
    auto ins_begin =
        sorted_unique.insert(sorted_unique.end(), to_add.begin(), to_add.end());
    std::sort(ins_begin, sorted_unique.end());
    std::inplace_merge(sorted_unique.begin(), ins_begin, sorted_unique.end());
    sorted_unique.erase(std::unique(sorted_unique.begin(), sorted_unique.end()),
                        sorted_unique.end());
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_rebuild_sat_model(
    const std::vector<std::size_t>& current_vertices) {
    // on calling this, the model should be reset/empty
    assert(m_sat_vertex_order.empty());
    assert(m_vertex_in_class.empty());
    assert(m_class_literals.empty());

    // rebuild SAT order
    m_sat_vertex_order.reserve(current_vertices.size());
    for (std::size_t vi : current_vertices) {
        m_vertex_info[vi].sat_index = m_sat_vertex_order.size();
        m_sat_vertex_order.push_back(vi);
    }

    // rebuild classes
    m_last_sat_clique_size = m_best_clique.size();
    std::size_t goal_classes = m_lower_bound;
    for (std::size_t ci = 0; ci < goal_classes; ++ci) {
        std::optional<std::size_t> clique_vertex = std::nullopt;
        if (ci < m_best_clique.size()) {
            clique_vertex = m_best_clique[ci];
        }
        p_rebuild_sat_color_class(clique_vertex);
    }

    // rebuild vertex constraints
    SatLit not_enough = m_incremental_solver.new_var();
    m_not_enough_colors = not_enough;
    m_not_enough_colors_is_for = goal_classes;
    for (std::size_t v : range(m_sat_vertex_order.size())) {
        p_add_vertex_color_constraint(v, not_enough);
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_rebuild_sat_color_class(
    std::optional<std::size_t> clique_vertex) {
    if (clique_vertex) {
        p_append_constrained_class_lits(*clique_vertex);
    } else {
        p_append_unconstrained_class_lits();
    }
    p_append_vertex_in_class();
    p_enforce_formula_on(m_class_literals.size() - 1);
}

template <typename IncrementalSATSolver>
auto CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_append_unconstrained_class_lits()
    -> std::vector<LitOrFixed>& {
    m_class_literals.emplace_back();
    auto& cl = m_class_literals.back();
    cl.reserve(m_n_all);
    for (Var v = 0; v < m_n_all; ++v) {
        Lit l = lit::positive_lit(v);
        if (m_empty_propagator.is_open(l)) {
            cl.emplace_back(p_lit(m_incremental_solver.new_var()));
        } else {
            cl.emplace_back(p_fixed(m_empty_propagator.is_true(l)));
        }
    }
    return cl;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::
    p_append_constrained_class_lits(std::size_t clique_vertex) {
    Vertex v = m_all_vertices[clique_vertex];
    push_vertex(m_empty_propagator, v);
    m_class_literals.emplace_back();
    auto& cl = m_class_literals.back();
    cl.reserve(m_n_all);
    for (Var v = 0; v < m_n_all; ++v) {
        Lit l = lit::positive_lit(v);
        if (m_empty_propagator.is_open(l)) {
            cl.emplace_back(m_incremental_solver.new_var());
        } else {
            cl.emplace_back(m_empty_propagator.is_true(l));
        }
    }
    m_empty_propagator.reset_to_zero();
}

template <typename IncrementalSATSolver>
auto CliqueSatDSaturSolver<IncrementalSATSolver>::p_append_vertex_in_class()
    -> std::vector<LitOrFixed>& {
    std::size_t class_index = m_vertex_in_class.size();
    m_vertex_in_class.emplace_back();
    assert(m_class_literals.size() == m_vertex_in_class.size());
    auto& v_in_c = m_vertex_in_class.back();
    v_in_c.reserve(m_sat_vertex_order.size());
    for (std::size_t vi : m_sat_vertex_order) {
        auto [l1, l2] = m_all_vertices[vi];
        LitOrFixed cl1 = p_lookup_sat_literal_in_class(l1, class_index);
        LitOrFixed cl2 = p_lookup_sat_literal_in_class(l2, class_index);
        v_in_c.push_back(p_and(cl1, cl2));
    }
    return v_in_c;
}

template <typename IncrementalSATSolver>
auto CliqueSatDSaturSolver<IncrementalSATSolver>::p_lookup_sat_literal_in_class(
    Lit l, std::size_t class_index) -> LitOrFixed {
    Var v = lit::var(l);
    bool neg = lit::negative(l);
    LitOrFixed v_in_c = m_class_literals[class_index][v];
    return neg ? p_negate(v_in_c) : v_in_c;
}

template <typename IncrementalSATSolver>
auto CliqueSatDSaturSolver<IncrementalSATSolver>::p_negate(LitOrFixed l)
    -> LitOrFixed {
    return std::visit(overloaded{[](bool& b) { return p_fixed(!b); },
                                 [](SatLit& l) { return p_lit(-l); }},
                      l);
}

template <typename IncrementalSATSolver>
auto CliqueSatDSaturSolver<IncrementalSATSolver>::p_and(LitOrFixed l1,
                                                        LitOrFixed l2)
    -> LitOrFixed {
    auto one_fixed = [&](SatLit b1, bool b2) {
        if (!b2)
            return p_fixed(false);
        return p_lit(b1);
    };
    return std::visit(
        overloaded{[&](bool& b1, bool& b2) { return p_fixed(b1 & b2); },
                   [&](bool& b1, SatLit& b2) { return one_fixed(b2, b1); },
                   [&](SatLit& b1, bool& b2) { return one_fixed(b1, b2); },
                   [&](SatLit& b1, SatLit& b2) {
                       SatLit nv = m_incremental_solver.new_var();
                       m_incremental_solver.add_short_clause(-nv, b1);
                       m_incremental_solver.add_short_clause(-nv, b2);
                       m_incremental_solver.add_short_clause(nv, -b1, -b2);
                       return p_lit(nv);
                   }},
        l1, l2);
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_enforce_formula_on(
    std::size_t class_index) {
    const auto& db = m_empty_propagator.db();
    for (auto [l1, l2] : db.binary_clauses()) {
        const Lit l[2] = {l1, l2};
        p_add_clause(l + 0, l + 2, class_index);
    }
    for (CRef c = 1, cn = db.literal_db_size(); c < cn; c = db.next_clause(c)) {
        auto r = db.lits_of(c);
        p_add_clause(r.begin(), r.end(), class_index);
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_add_clause(
    const Lit* begin, const Lit* end, std::size_t class_index) {
    m_clause_buffer.clear();
    for (Lit l : IteratorRange{begin, end}) {
        LitOrFixed lf = p_lookup_sat_literal_in_class(l, class_index);
        bool clause_satisfied =
            std::visit(overloaded{[&](bool& b) -> bool { return b; },
                                  [&](SatLit& s) -> bool {
                                      m_clause_buffer.push_back(s);
                                      return false;
                                  }},
                       lf);
        if (clause_satisfied)
            return;
    }
    if (m_clause_buffer.empty())
        throw std::logic_error("UNSAT class!");
    m_incremental_solver.add_clause(m_clause_buffer.begin(),
                                    m_clause_buffer.end());
}

template <typename IncrementalSATSolver>
std::optional<bool>
CliqueSatDSaturSolver<IncrementalSATSolver>::p_run_sat_model() {
    std::vector<SatLit> assumptions;
    assumptions.push_back(-*m_not_enough_colors);
    std::optional<bool> result = m_incremental_solver.solve(assumptions);
    if (!result || !*result) {
        return result;
    }
    p_extract_sat_model(m_incremental_solver.get_model());
    return true;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_add_vertex_color_constraint(
    std::size_t sat_index, SatLit not_enough) {
    m_clause_buffer.clear();
    for (std::size_t c = 0; c < m_vertex_in_class.size(); ++c) {
        LitOrFixed v_in_c = m_vertex_in_class[c][sat_index];
        bool is_satisfied =
            std::visit(overloaded{[&](bool& b) -> bool { return b; },
                                  [&](SatLit& l) -> bool {
                                      m_clause_buffer.push_back(l);
                                      return false;
                                  }},
                       v_in_c);
        if (is_satisfied)
            return;
    }
    m_clause_buffer.push_back(not_enough);
    m_incremental_solver.add_clause(m_clause_buffer.begin(),
                                    m_clause_buffer.end());
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_extend_sat_model(
    const std::vector<std::size_t>& current_vertices) {
    std::vector<std::size_t> new_vertices;
    std::size_t old_sat_count = m_sat_vertex_order.size();
    for (std::size_t v : current_vertices) {
        if (m_vertex_info[v].sat_index ==
            std::numeric_limits<std::size_t>::max())
        {
            std::size_t sindex = m_sat_vertex_order.size();
            m_sat_vertex_order.push_back(v);
            m_vertex_info[v].sat_index = sindex;
            new_vertices.push_back(v);
        }
    }
    std::size_t num_classes = m_vertex_in_class.size();
    for (std::size_t c = 0; c < num_classes; ++c) {
        auto& v_in_c = m_vertex_in_class[c];
        for (std::size_t new_vertex : new_vertices) {
            auto [l1, l2] = m_all_vertices[new_vertex];
            LitOrFixed sl1 = p_lookup_sat_literal_in_class(l1, c);
            LitOrFixed sl2 = p_lookup_sat_literal_in_class(l2, c);
            v_in_c.push_back(p_and(sl1, sl2));
        }
    }
    for (std::size_t new_vertex :
         range(old_sat_count, old_sat_count + new_vertices.size()))
    {
        p_add_vertex_color_constraint(new_vertex, *m_not_enough_colors);
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_extend_sat_model_to_classes() {
    std::size_t old_classes = m_vertex_in_class.size();
    if (old_classes == m_lower_bound)
        return;
    for (std::size_t c = m_vertex_in_class.size(); c < m_lower_bound; ++c) {
        std::optional<std::size_t> clique_vertex = std::nullopt;
        if (c < m_best_clique.size()) {
            clique_vertex = m_best_clique[c];
        }
        p_rebuild_sat_color_class(clique_vertex);
    }
    p_check_not_enough_colors();
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_check_not_enough_colors() {
    if (!m_not_enough_colors || m_not_enough_colors_is_for != m_lower_bound) {
        if (m_not_enough_colors)
            m_incremental_solver.fix(*m_not_enough_colors);
        m_not_enough_colors = m_incremental_solver.new_var();
        m_not_enough_colors_is_for = m_lower_bound;
        for (std::size_t v : range(m_sat_vertex_order.size())) {
            p_add_vertex_color_constraint(v, *m_not_enough_colors);
        }
    }
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_zero_choice_run_clique(
    const std::vector<std::size_t>& current_vertices) {
    std::size_t old_clique_size = m_best_clique.size();
    if (m_clique_candidate_callback) {
        auto candidate = m_clique_candidate_callback();
        if (candidate.size() > old_clique_size) {
            auto candidate_indices = p_clique_to_indices(candidate, false);
            if (candidate_indices.size() > old_clique_size) {
                m_best_clique = std::move(candidate_indices);
                if (m_best_clique.size() > m_lower_bound) {
                    m_lower_bound = m_best_clique.size();
                    m_lower_bound_subgraph = m_best_clique;
                    p_call_lb_callback();
                }
                if (!m_existing_clique_vars.empty()) {
                    p_ensure_best_clique_in_model();
                }
            }
        }
    }
    if (m_existing_clique_vars.empty()) {
        if (!p_first_run_clique_model(current_vertices))
            return false;
    } else {
        if (m_existing_clique_var_value_count !=
            m_existing_clique_var_indices.size())
        {
            if (!p_solve_full_relaxation())
                return false;
        }
        if (p_clique_model_cclass_cuts()) {
            if (!p_solve_full_relaxation())
                return false;
        }
    }
    if (!p_use_clique_model(current_vertices))
        return false;
    return old_clique_size < m_best_clique.size();
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_run_until_satisfiable() {
    for (;;) {
        auto res = p_run_sat_model();
        if (!res)
            return false;
        if (!*res) {
            m_lower_bound += 1;
            m_lower_bound_subgraph = m_sat_vertex_order;
            p_call_lb_callback();
            if (m_lower_bound == m_covering_assignments.size()) {
                return false;
            }
            p_extend_sat_model_to_classes();
        } else {
            return true;
        }
    }
}

template <typename IncrementalSATSolver>
template <typename ModelMap>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_extract_sat_model(
    const ModelMap& model_map) {
    const std::size_t sn = m_sat_vertex_order.size();
    const std::size_t cn = m_vertex_in_class.size();
    p_reset_color_classes(cn);
    std::vector<std::size_t> class_with_vertex;
    for (std::size_t si = 0; si < sn; ++si) {
        for (std::size_t c = 0; c < cn; ++c) {
            LitOrFixed l = m_vertex_in_class[c][si];
            if (p_lookup_in_model(model_map, l)) {
                class_with_vertex.push_back(c);
                break;
            }
        }
    }
    if (class_with_vertex.size() != sn) {
        std::cerr << "SAT model does not assign a color to each vertex!"
                  << std::endl;
        std::cerr << "SAT solver name: " << m_incremental_solver.name()
                  << std::endl;
        throw std::logic_error(
            "SAT model does not assign a color to each vertex!");
    }
    for (std::size_t si = 0; si < sn; ++si) {
        std::size_t vi = m_sat_vertex_order[si];
        m_color_classes[class_with_vertex[si]].add_vertex(this, vi);
        if (!m_in_currently_colored[vi]) {
            m_in_currently_colored[vi].set();
            m_currently_colored.push_back(vi);
        }
    }
}

template <typename IncrementalSATSolver>
template <typename ModelMap>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_lookup_in_model(
    const ModelMap& map, LitOrFixed l) const {
    return std::visit(overloaded{[&](bool& b) -> bool { return b; },
                                 [&](SatLit& sl) -> bool { return map[sl]; }},
                      l);
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_reset_color_classes(
    std::size_t num_classes) {
    const std::size_t cold = m_color_classes.size();
    if (num_classes > cold) {
        for (std::size_t c = 0; c < cold; ++c) {
            m_color_classes[c].propagator.reset_to_zero();
        }
        for (std::size_t c = cold; c < num_classes; ++c) {
            m_color_classes.emplace_back(m_empty_propagator, c);
        }
        for (VertexInfo& v : m_vertex_info) {
            v.open_classes.resize(num_classes);
            v.open_classes.set();
            v.num_open_classes = num_classes;
            v.in_some_class = false;
        }
    } else {
        for (std::size_t c = 0; c < num_classes; ++c) {
            m_color_classes[c].propagator.reset_to_zero();
        }
        for (VertexInfo& v : m_vertex_info) {
            v.open_classes.set();
            v.num_open_classes = num_classes;
            v.in_some_class = false;
        }
    }
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_complete_all_classes() {
    struct EmptyHandler {};
    EmptyHandler handler;
    ClassCompleter<EmptyHandler> completer{m_n_concrete, m_n_all, &handler};
    bool result = true;
    for (std::size_t cindex : range(m_color_classes.size())) {
        SharedDBPropagator& prop = m_color_classes[cindex].propagator;
        if (!completer.complete_class(prop)) {
            result = false;
        }
    }
    if (result)
        return true;
    p_reset();
    return false;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_reset() {
    for (VertexInfo& v : m_vertex_info) {
        v.open_classes.set();
        v.num_open_classes = m_color_classes.size();
        v.in_some_class = false;
    }
    p_reset_sat_model();
    m_empty_propagator.incorporate_or_throw();
    for (ColorClass& cc : m_color_classes) {
        cc.propagator.reset_or_throw();
    }
    m_currently_colored = m_best_clique;
    m_in_currently_colored.reset();
    m_zero_candidates.clear();
    m_one_candidates.clear();
    m_in_clique_uncolored.clear();
    for (std::size_t vi : m_currently_colored) {
        m_in_currently_colored[vi].set();
    }
    for (std::size_t ci : range(m_best_clique.size())) {
        std::size_t vi = m_best_clique[ci];
        m_color_classes[ci].add_vertex(this, vi);
    }
}

template <typename IncrementalSATSolver>
PartialSolution
CliqueSatDSaturSolver<IncrementalSATSolver>::get_partial_solution() const {
    auto transformer =
        [](const ColorClass& color) -> const SharedDBPropagator& {
        return color.propagator;
    };
    auto begin = boost::iterators::make_transform_iterator(
        m_color_classes.begin(), +transformer);
    auto end = boost::iterators::make_transform_iterator(m_color_classes.end(),
                                                         +transformer);
    return PartialSolution(m_n_all, m_inf_map, begin, end);
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_first_run_clique_model(
    const std::vector<std::size_t>& current_vertices) {
    p_setup_clique_model();
    for (const DynamicBitset& ca : m_covering_assignments) {
        p_clq_add_complete_constraint(ca);
    }
    for (const ColorClass& cc : m_color_classes) {
        p_clq_add_partial_constraint(cc.propagator);
    }
    for (std::size_t vertex : m_best_clique) {
        m_last_vertex_addition_source = "first_run_clique_model 1";
        p_clq_add_vertex(vertex);
    }
    auto added_vertices =
        sample_from_range(current_vertices, 1000, sammy::rng());
    for (std::size_t vertex : added_vertices) {
        m_last_vertex_addition_source = "first_run_clique_model 2";
        p_clq_add_vertex(vertex);
    }
    if (!p_solve_full_relaxation()) {
        return false;
    }
    return true;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_clq_add_complete_constraint(
    const DynamicBitset& complete_assignment) {
    m_buffer_expr.clear();
    for (std::size_t i = 0, ncv = m_existing_clique_vars.size(); i < ncv; ++i) {
        std::size_t vi = m_existing_clique_var_indices[i];
        Vertex v = m_all_vertices[vi];
        if (lit::is_true_in(v.first, complete_assignment) &&
            lit::is_true_in(v.second, complete_assignment))
        {
            m_buffer_expr += m_existing_clique_vars[i];
        }
    }
    m_clq_constr_complete_assignments.push_back(complete_assignment);
    m_clq_constr_complete_constraints.push_back(
        m_clique_model.addConstr(m_buffer_expr <= 1));
}

template <typename IncrementalSATSolver>
template <typename Propagator>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_clq_add_partial_constraint(
    Propagator&& partial_assignment, bool expr_in_buffer) {
    if (!expr_in_buffer) {
        m_buffer_expr.clear();
        for (std::size_t i = 0, ncv = m_existing_clique_vars.size(); i < ncv;
             ++i)
        {
            std::size_t vi = m_existing_clique_var_indices[i];
            Vertex v = m_all_vertices[vi];
            if (partial_assignment.is_true(v.first) &&
                partial_assignment.is_true(v.second))
            {
                m_buffer_expr += m_existing_clique_vars[i];
            }
        }
    }
    m_clq_constr_partial_assignments.emplace_back(
        std::forward<Propagator>(partial_assignment));
    m_clq_constr_partial_constraints.push_back(
        m_clique_model.addConstr(m_buffer_expr <= 1));
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_clq_add_vertex(
    std::size_t vertex) {
    if (m_present_in_clique_model[vertex])
        return;
    m_column_buffer.clear();
    Vertex v = m_all_vertices[vertex];
    std::size_t clique_model_index = m_existing_clique_vars.size();
    for (std::size_t complete_i :
         range(m_clq_constr_complete_assignments.size()))
    {
        const auto& assignment = m_clq_constr_complete_assignments[complete_i];
        if (lit::is_true_in(v.first, assignment) &&
            lit::is_true_in(v.second, assignment))
        {
            m_column_buffer.addTerm(
                1.0, m_clq_constr_complete_constraints[complete_i]);
        }
    }
    for (std::size_t partial_i : range(m_clq_constr_partial_assignments.size()))
    {
        const auto& assignment = m_clq_constr_partial_assignments[partial_i];
        if (assignment.is_true(v.first) && assignment.is_true(v.second)) {
            m_column_buffer.addTerm(
                1.0, m_clq_constr_partial_constraints[partial_i]);
        }
    }
    GRBVar var =
        m_clique_model.addVar(0.0, 1.0, 1.0, GRB_CONTINUOUS, m_column_buffer);
    m_present_in_clique_model[vertex].set();
    m_existing_clique_var_indices.push_back(vertex);
    m_existing_clique_vars.push_back(var);
    m_all_clique_vars[vertex] = var;
    m_clique_model_vertices_containing_concrete_literal[v.first].push_back(
        clique_model_index);
    m_clique_model_vertices_containing_concrete_literal[v.second].push_back(
        clique_model_index);
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_clq_extended_partial(
    std::size_t partial_index, std::size_t old_trail_size) {
    const auto& propagator = m_clq_constr_partial_assignments[partial_index];
    const auto& trail = propagator.get_trail();
    const Lit nclit = 2 * m_n_concrete;
    GRBConstr constr = m_clq_constr_partial_constraints[partial_index];
    for (Lit lnew : IteratorRange{trail.begin() + old_trail_size, trail.end()})
    {
        if (lnew >= nclit)
            continue;
        for (std::size_t clique_model_index :
             m_clique_model_vertices_containing_concrete_literal[lnew])
        {
            std::size_t vertex_index =
                m_existing_clique_var_indices[clique_model_index];
            Vertex v = m_all_vertices[vertex_index];
            if (propagator.is_true(v.first) && propagator.is_true(v.second)) {
                m_clique_model.chgCoeff(
                    constr, m_existing_clique_vars[clique_model_index], 1.0);
            }
        }
    }
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::
    p_clq_extend_partial_by_vertex(std::size_t partial_index,
                                   std::size_t vertex) {
    Vertex v = m_all_vertices[vertex];
    auto& propagator = m_clq_constr_partial_assignments[partial_index];
    const std::size_t old_trail_size = propagator.get_trail().size();
    int p = push_vertex(propagator, v);
    if (p < 0)
        return false;
    if (p == 0)
        return true;
    p_clq_extended_partial(partial_index, old_trail_size);
    return true;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<
    IncrementalSATSolver>::GurobiCallbackCheckAbort::callback() {
    if (that->m_aborted.load())
        this->abort();
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_setup_clique_model() {
    // thread count, sense and experimentally determined method/presolve level
    m_clique_model.set(GRB_IntParam_Threads, 1);
    m_clique_model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    m_clique_model.set(GRB_IntParam_Method, GRB_METHOD_DUAL);
    m_clique_model.set(GRB_IntParam_Presolve, 1);
    m_check_abort_cb = std::make_unique<GurobiCallbackCheckAbort>(this);
    m_clique_model.setCallback(m_check_abort_cb.get());
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_solve_relaxation() {
    bool repeat = true;
    while (repeat) {
        repeat = false;
        m_clique_model.optimize();
        int status = m_clique_model.get(GRB_IntAttr_Status);
        if (status == GRB_INTERRUPTED) {
            m_aborted.store(true);
            return false;
        }
        if (status != GRB_OPTIMAL) {
            throw std::logic_error(
                "Unexpected status after solving relaxation!");
        }
        m_existing_clique_var_value_count = m_existing_clique_vars.size();
        m_last_vertex_addition_source = "cleared after solving";
        m_existing_clique_var_values.reset(
            m_clique_model.get(GRB_DoubleAttr_X, m_existing_clique_vars.data(),
                               int(m_existing_clique_vars.size())));
        m_clq_constr_complete_dual_values.reset(m_clique_model.get(
            GRB_DoubleAttr_Pi, m_clq_constr_complete_constraints.data(),
            int(m_clq_constr_complete_constraints.size())));
        m_clq_constr_partial_dual_values.reset(m_clique_model.get(
            GRB_DoubleAttr_Pi, m_clq_constr_partial_constraints.data(),
            int(m_clq_constr_partial_constraints.size())));
        m_last_clique_objective = m_clique_model.get(GRB_DoubleAttr_ObjVal);
        m_last_rounded_objective = static_cast<std::size_t>(
            std::floor(m_last_clique_objective + 0.001));
        m_clq_ordered_solution.clear();
        for (std::size_t i = 0; i < m_existing_clique_vars.size(); ++i) {
            assert(i < m_existing_clique_var_value_count);
            double x = m_existing_clique_var_values[i];
            if (x >= 1.0e-4) {
                m_clq_ordered_solution.emplace_back(OrderedSolutionValue{i, x});
            }
        }
        std::sort(m_clq_ordered_solution.begin(), m_clq_ordered_solution.end());
        if (m_aborted.load())
            return false;
        if (p_round_greedy()) {
            repeat = true;
        }
    }
    return !m_aborted.load();
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_solve_full_relaxation() {
    for (;;) {
        if (!p_solve_relaxation())
            return false;
        if (m_last_rounded_objective <= m_best_clique.size())
            return true;
        if (!p_identify_violated_nonedges())
            return true;
        p_prohibit_violated_nonedges();
    }
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_identify_violated_nonedges() {
    m_clq_violated_nonedges.clear();
    for (OrderedSolutionIter i = m_clq_ordered_solution.begin(),
                             e = m_clq_ordered_solution.end();
         i != e; ++i)
    {
        double weight = i->value;
        double thresh = 1.01 - weight;
        if (weight < thresh)
            break;
        std::size_t clique_index = i->clique_model_index;
        std::size_t vertex_index = m_existing_clique_var_indices[clique_index];
        for (OrderedSolutionIter j = std::next(i); j != e; ++j) {
            double w_j = j->value;
            if (w_j < thresh)
                break;
            std::size_t vertex_index2 =
                m_existing_clique_var_indices[j->clique_model_index];
            if (p_is_nonedge(vertex_index, vertex_index2)) {
                m_clq_violated_nonedges.emplace_back(i, j);
            }
        }
    }
    return !m_clq_violated_nonedges.empty();
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_is_nonedge(
    std::size_t vertex1, std::size_t vertex2) {
    if (m_adjacency_matrix[vertex1][vertex2])
        return false;
    if (m_definitive_nonedges[vertex1][vertex2])
        return true;
    if (!push_vertex_pair(m_empty_propagator, m_all_vertices[vertex1],
                          m_all_vertices[vertex2]))
    {
        m_adjacency_matrix[vertex1][vertex2] = true;
        m_adjacency_matrix[vertex2][vertex1] = true;
        return false;
    } else {
        m_empty_propagator.reset_to_zero();
        m_definitive_nonedges[vertex1][vertex2] = true;
        m_definitive_nonedges[vertex2][vertex1] = true;
        return true;
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_prohibit_violated_nonedges() {
    m_covered_violated_nonedges.assign(m_clq_violated_nonedges.size(), false);
    auto last_first = m_clq_ordered_solution.cend();
    const auto nve = m_clq_violated_nonedges.size();
    std::size_t num_new_constraints = 0;
    for (std::size_t i = 0; i < nve; ++i) {
        auto [vit, wit] = m_clq_violated_nonedges[i];
        if (vit == last_first || m_covered_violated_nonedges[i])
            continue;
        if (p_prohibit_extend_existing(vit, wit)) {
            m_covered_violated_nonedges[i].set();
            continue;
        }
        last_first = vit;
        if (++num_new_constraints > CAP_PROHIBIT_VIOLATED_NUM_NEW_CONSTRAINTS)
            continue;
        std::size_t vind =
            m_existing_clique_var_indices[vit->clique_model_index];
        std::size_t wind =
            m_existing_clique_var_indices[wit->clique_model_index];
        if (!push_vertex_pair(m_empty_propagator, m_all_vertices[vind],
                              m_all_vertices[wind]))
        {
            throw std::logic_error(
                "Incorrectly identified edge as violated nonedge");
        }
        m_covered_violated_nonedges[i].set();
        p_extend_nonedge_cover_greedy(i + 1, nve);
        p_extend_nonedge_cover_greedy(0, i);
        p_clq_add_partial_constraint(m_empty_propagator);
        m_empty_propagator.reset_to_zero();
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_extend_nonedge_cover_greedy(
    std::size_t begin_index, std::size_t end_index) {
    for (std::size_t j = begin_index; j < end_index; ++j) {
        if (m_covered_violated_nonedges[j])
            continue;
        auto [vit, wit] = m_clq_violated_nonedges[j];
        std::size_t ci1 = vit->clique_model_index;
        std::size_t ci2 = wit->clique_model_index;
        Vertex v1 = m_all_vertices[m_existing_clique_var_indices[ci1]];
        Vertex v2 = m_all_vertices[m_existing_clique_var_indices[ci2]];
        if (push_vertex_pair(m_empty_propagator, v1, v2)) {
            m_covered_violated_nonedges[j].set();
        }
    }
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_prohibit_extend_existing(
    OrderedSolutionIter vit, OrderedSolutionIter wit) {
    std::size_t ci1 = vit->clique_model_index;
    std::size_t ci2 = wit->clique_model_index;
    Vertex v1 = m_all_vertices[m_existing_clique_var_indices[ci1]];
    Vertex v2 = m_all_vertices[m_existing_clique_var_indices[ci2]];
    for (std::size_t i = 0, ninc = m_clq_constr_partial_constraints.size();
         i < ninc; ++i)
    {
        SharedDBPropagator& prop = m_clq_constr_partial_assignments[i];
        if (prop.is_false(v1.first) || prop.is_false(v1.second) ||
            prop.is_false(v2.first) || prop.is_false(v2.second))
            continue;
        std::size_t trail_length = prop.get_trail().size();
        if (push_vertex_pair(prop, v1, v2)) {
            p_clq_extended_partial(i, trail_length);
            return true;
        }
    }
    return false;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_round_greedy() {
    // track the vertices in the clique
    m_vertex_set_buffer.clear();
    // track possible extensions of the clique
    m_vertex_bitset_buffer.assign(m_all_vertices.size(), true);
    // greedy rounding procedure - we use the basic G_1 definition
    // that has often been extended, at least on the
    // vertices in the ordered solution, by violated non-edge detection
    for (const auto& ov : m_clq_ordered_solution) {
        std::size_t ci = ov.clique_model_index;
        std::size_t vi = m_existing_clique_var_indices[ci];
        if (!m_vertex_bitset_buffer[vi])
            continue;
        assert(!m_adjacency_matrix[vi][vi]);
        m_vertex_bitset_buffer &= m_adjacency_matrix[vi];
        m_vertex_set_buffer.push_back(vi);
    }
    p_randomly_extend_clique();
    if (m_vertex_set_buffer.size() > m_best_clique.size()) {
        m_best_clique.swap(m_vertex_set_buffer);
        if (m_best_clique.size() > m_lower_bound) {
            m_lower_bound = m_best_clique.size();
            m_lower_bound_subgraph = m_best_clique;
            p_call_lb_callback();
        }
        p_ensure_best_clique_in_model();
        return true;
    }
    return false;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_randomly_extend_clique() {
    auto& rng = sammy::rng();
    const std::size_t n = m_all_vertices.size();
    std::uniform_int_distribution<std::size_t> indices(0, n - 1);
    for (std::size_t trial_counter = 1; trial_counter <= 100; ++trial_counter) {
        std::size_t vi = indices(rng);
        if (!m_vertex_bitset_buffer[vi])
            continue;
        assert(!m_adjacency_matrix[vi][vi]);
        m_vertex_bitset_buffer &= m_adjacency_matrix[vi];
        m_vertex_set_buffer.push_back(vi);
        trial_counter = 0;
    }
    bool any = m_vertex_bitset_buffer.any();
    while (any) {
        any = false;
        for (std::size_t vi : m_vertex_bitset_buffer.ones()) {
            any = true;
            assert(!m_adjacency_matrix[vi][vi]);
            m_vertex_bitset_buffer &= m_adjacency_matrix[vi];
            m_vertex_set_buffer.push_back(vi);
            break;
        }
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_ensure_best_clique_in_model() {
    for (std::size_t vi : m_best_clique) {
        m_last_vertex_addition_source = "ensure_best_clique_in_model";
        p_clq_add_vertex(vi);
    }
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_clique_model_cclass_cuts() {
    bool found = false;
    for (const auto& cc : m_color_classes) {
        const auto& prop = cc.propagator;
        double sum = 0.0;
        m_buffer_expr.clear();
        for (std::size_t i = 0, nc = m_existing_clique_var_indices.size();
             i != nc; ++i)
        {
            std::size_t vi = m_existing_clique_var_indices[i];
            Vertex v = m_all_vertices[vi];
            if (prop.is_false(v.first) || prop.is_false(v.second))
                continue;
            assert(i < m_existing_clique_var_value_count);
            sum += m_existing_clique_var_values[i];
            m_buffer_expr += m_existing_clique_vars[i];
        }
        if (sum >= 1.01) {
            p_clq_add_partial_constraint(prop, true);
            found = true;
        }
    }
    return found;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_use_clique_model(
    const std::vector<std::size_t>& current_vertices) {
    std::size_t old_clique_size = m_best_clique.size();
    for (std::size_t i = 0; i < 5; ++i) {
        if (m_last_rounded_objective > m_best_clique.size() &&
            !p_cheap_cut_rounds())
            return false;
        if (m_last_rounded_objective > m_best_clique.size() &&
            !p_expensive_cut_rounds())
            return false;
        if (!p_limited_price_vertices(current_vertices))
            return false;
    }
    if (!p_unlimited_price_vertices(current_vertices))
        return false;
    if (m_last_rounded_objective > m_best_clique.size() &&
        !p_cheap_cut_rounds())
        return false;
    if (m_last_rounded_objective > m_best_clique.size() &&
        !p_expensive_cut_rounds())
        return false;
    return !m_aborted.load() && m_best_clique.size() > old_clique_size;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_cheap_cut_rounds() {
    bool repeat = true;
    std::size_t total_it = 0;
    while (repeat && m_last_rounded_objective > m_best_clique.size()) {
        repeat = false;
        double gap_before = p_clique_gap();
        for (std::size_t it = 0; it < CHEAP_CUT_ROUNDS_PER_GAP_CHECK;
             ++it, ++total_it)
        {
            if (!p_cheap_cut_round()) {
                return !m_aborted.load();
            }
            if (m_last_rounded_objective <= m_best_clique.size())
                return true;
        }
        double gap_after = p_clique_gap();
        repeat = (gap_before - gap_after >= CHEAP_CUT_GAP_REDUCTION_REQUIRED);
        repeat &= (total_it < CHEAP_CUT_ROUNDS_HARD_CAP);
    }
    return true;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_cheap_cut_round() {
    if (p_greedy_add_to_cuts() || p_greedy_generate_cuts() ||
        p_clique_model_cclass_cuts())
    {
        return p_solve_full_relaxation();
    }
    return false;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_greedy_add_to_cuts() {
    bool found = false;
    for (std::size_t partial_i : range(m_clq_constr_partial_assignments.size()))
    {
        auto& prop = m_clq_constr_partial_assignments[partial_i];
        double total_value = 0.0;
        for (const OrderedSolutionValue& ov : m_clq_ordered_solution) {
            std::size_t ci = ov.clique_model_index;
            std::size_t vi = m_existing_clique_var_indices[ci];
            Vertex v = m_all_vertices[vi];
            if (!can_push(prop, v))
                continue;
            total_value += ov.value;
        }
        if (total_value < 1.01)
            continue;
        total_value = 0.0;
        std::size_t old_trail_length = prop.get_trail().size();
        for (const OrderedSolutionValue& ov : m_clq_ordered_solution) {
            std::size_t ci = ov.clique_model_index;
            std::size_t vi = m_existing_clique_var_indices[ci];
            Vertex v = m_all_vertices[vi];
            if (push_vertex(prop, v) >= 0)
                total_value += ov.value;
        }
        p_clq_extended_partial(partial_i, old_trail_length);
        if (total_value >= 1.01)
            found = true;
    }
    return found;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_greedy_generate_cuts() {
    m_greedy_cut_candidates.clear();
    for (OrderedSolutionIter b = m_clq_ordered_solution.begin(), i = b,
                             e = m_clq_ordered_solution.end();
         i != e; ++i)
    {
        m_vertex_set_buffer.clear();
        assert(m_empty_propagator.get_current_level() == 0);
        double v = p_extend_greedy_cut(i, e);
        v += p_extend_greedy_cut(b, i);
        if (v >= 1.01 && m_greedy_cut_candidates.would_push(v)) {
            m_greedy_cut_candidates.emplace(v, m_vertex_set_buffer);
        }
        m_empty_propagator.reset_to_zero();
    }
    for (const GreedyCutCandidate& candidate :
         m_greedy_cut_candidates.elements())
    {
        for (std::size_t vertex : candidate.vertex_indices) {
            Vertex v = m_all_vertices[vertex];
            if (push_vertex(m_empty_propagator, v) < 0) {
                throw std::logic_error(
                    "Incorrect greedy cut candidate: produced conflict!");
            }
        }
        p_clq_add_partial_constraint(m_empty_propagator);
        m_empty_propagator.reset_to_zero();
    }
    return !m_greedy_cut_candidates.elements().empty();
}

template <typename IncrementalSATSolver>
double CliqueSatDSaturSolver<IncrementalSATSolver>::p_extend_greedy_cut(
    OrderedSolutionIter begin, OrderedSolutionIter end) {
    double value = 0.0;
    for (OrderedSolutionIter i = begin; i != end; ++i) {
        std::size_t ci = i->clique_model_index;
        std::size_t vi = m_existing_clique_var_indices[ci];
        Vertex v = m_all_vertices[vi];
        if (push_vertex(m_empty_propagator, v) < 0)
            continue;
        value += i->value;
        m_vertex_set_buffer.push_back(vi);
    }
    return value;
}

template <typename IncrementalSATSolver>
double CliqueSatDSaturSolver<IncrementalSATSolver>::p_weight_in_dual(Vertex v) {
    double value = 0.0;
    for (std::size_t complete_i :
         range(m_clq_constr_complete_assignments.size()))
    {
        const auto& assignment = m_clq_constr_complete_assignments[complete_i];
        if (lit::is_true_in(v.first, assignment) &&
            lit::is_true_in(v.second, assignment))
        {
            value += m_clq_constr_complete_dual_values[complete_i];
            if (value >= 1.02)
                return 1.02;
        }
    }
    for (std::size_t partial_i : range(m_clq_constr_partial_assignments.size()))
    {
        const auto& assignment = m_clq_constr_partial_assignments[partial_i];
        if (assignment.is_true(v.first) && assignment.is_true(v.second)) {
            value += m_clq_constr_partial_dual_values[partial_i];
            if (value >= 1.02)
                return 1.02;
        }
    }
    return value;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::
    p_price_collect_good_and_possible(
        const std::vector<std::size_t>& vertices) {
    m_pricing_good_vertices.clear();
    m_pricing_possible_vertices.clear();
    for (std::size_t vertex : vertices) {
        if (m_present_in_clique_model[vertex])
            continue;
        Vertex v = m_all_vertices[vertex];
        double w = p_weight_in_dual(v);
        if (w >= 1.01)
            continue;
        if (w <= 0.99) {
            m_pricing_good_vertices.emplace_back(vertex, w);
        } else {
            m_pricing_possible_vertices.emplace_back(vertex, w);
        }
    }
    return !m_pricing_good_vertices.empty() ||
           !m_pricing_possible_vertices.empty();
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_limited_price_vertices(
    const std::vector<std::size_t>& vertices) {
    std::size_t goal_vertices = (std::min)(
        vertices.size(),
        (std::max)(std::size_t(10),
                   std::size_t(0.2 * m_existing_clique_vars.size())));
    if (m_aborted.load())
        return false;
    if (!p_price_collect_good_and_possible(vertices))
        return false;
    if (m_aborted.load())
        return false;
    p_pricing_select_with_goal(goal_vertices);
    p_pricing_add_vertices();
    return p_solve_full_relaxation();
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_pricing_select_with_goal(
    std::size_t goal_vertices) {
    auto& good = m_pricing_good_vertices;
    auto& possible = m_pricing_possible_vertices;
    m_vertex_set_buffer.clear();
    if (good.empty()) {
        // no good vertices; take possible vertices
        if (goal_vertices > possible.size())
            goal_vertices = possible.size();
        std::transform(possible.begin(), possible.begin() + goal_vertices,
                       std::back_inserter(m_vertex_set_buffer),
                       [](const PricingEntry& e) { return e.vertex_index; });
        return;
    }
    if (good.size() <= goal_vertices) {
        // not enough good vertices; take all
        std::transform(good.begin(), good.end(),
                       std::back_inserter(m_vertex_set_buffer),
                       [](const PricingEntry& e) { return e.vertex_index; });
        return;
    }
    auto goal_iter = good.begin() + goal_vertices;
    std::nth_element(good.begin(), goal_iter, good.end(),
                     [](const PricingEntry& e1, const PricingEntry& e2) {
                         return e1.dual_weight > e2.dual_weight;
                     });
    std::transform(good.begin(), goal_iter,
                   std::back_inserter(m_vertex_set_buffer),
                   [](const PricingEntry& e) { return e.vertex_index; });
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_pricing_add_vertices() {
    if (m_vertex_set_buffer.empty())
        return;
    m_last_vertex_addition_source = "pricing_add_vertices";
    for (std::size_t v : m_vertex_set_buffer)
        p_clq_add_vertex(v);
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_unlimited_price_vertices(
    const std::vector<std::size_t>& vertices) {
    if (m_aborted.load())
        return false;
    if (!p_price_collect_good_and_possible(vertices))
        return false;
    if (m_aborted.load())
        return false;
    auto* source = &m_pricing_good_vertices;
    if (m_pricing_good_vertices.empty())
        source = &m_pricing_possible_vertices;
    m_vertex_set_buffer.clear();
    std::transform(source->begin(), source->end(),
                   std::back_inserter(m_vertex_set_buffer),
                   [](const PricingEntry& e) { return e.vertex_index; });
    p_pricing_add_vertices();
    return p_solve_full_relaxation();
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_call_lb_callback() {
    if (!m_lower_bound_callback)
        return;
    std::vector<Vertex> buffer;
    buffer.reserve(m_lower_bound_subgraph.size());
    std::transform(m_lower_bound_subgraph.begin(), m_lower_bound_subgraph.end(),
                   std::back_inserter(buffer),
                   [this](std::size_t vi) { return m_all_vertices[vi]; });
    m_lower_bound_callback(m_lower_bound, buffer);
}

} // namespace sammy

#endif
