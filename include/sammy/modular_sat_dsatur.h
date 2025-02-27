#ifndef SAMMY_MODULAR_SAT_DSATUR_H_INCLUDED_
#define SAMMY_MODULAR_SAT_DSATUR_H_INCLUDED_

#include "best_k.h"
#include "incremental_gurobi_g2_clique_solver.h"

namespace sammy {

class DSaturStylePartialConfigurations {
    using FullVertexIndex = std::size_t;

    /**
     * Information tracked for each vertex.
     */
    struct VertexInfo {
        std::size_t num_open_classes;
        std::size_t degree;
        DynamicBitset open_classes;
        bool in_some_class;
    };

    struct CompareVertexIndexByInfo {
        explicit CompareVertexIndexByInfo(
            DSaturStylePartialConfigurations* that) noexcept
            : that(that) {}

        bool operator()(FullVertexIndex v1, FullVertexIndex v2) const noexcept {
            // return true if 'v1' has higher priority than 'v2'
            const auto& vi1 = that->m_all_vertex_infos[v1];
            const auto& vi2 = that->m_all_vertex_infos[v2];
            if (vi1.in_some_class != vi2.in_some_class) {
                return !vi1.in_some_class;
            }
            if (vi1.num_open_classes != vi2.num_open_classes) {
                return vi1.num_open_classes < vi2.num_open_classes;
            }
            return vi1.degree > vi2.degree;
        }

        DSaturStylePartialConfigurations* that;
    };

    LazyG2AdjacencyMatrix* m_adjacency_matrix;
    std::vector<VertexInfo> m_all_vertex_infos;
    std::vector<SharedDBPropagator> m_partial_configurations;
    std::vector<FullVertexIndex> m_vertices_with_low_candidate_count[4];
    BestK<FullVertexIndex, CompareVertexIndexByInfo>
        m_best_with_high_candidate_count;
};

template <typename IncrementalSatSolver> class ModularCliqueSatDSatur {
  public:
    using SatModelIndex = std::size_t;
    using FullVertexIndex = std::size_t;
    using SLit = sammy::lit::Lit;
    using Lit = typename IncrementalSatSolver::Lit;
    using LitOrVal = std::variant<bool, Lit>;

    explicit ModularCliqueSatDSatur(
        IncrementalGurobiG2CliqueSolver* clique_solver,
        const std::vector<DynamicBitset>& best_covering_assignment)
        : m_clique_solver(clique_solver),
          m_adjacency_matrix(&clique_solver->adjacency_matrix()),
          m_initial_lb(clique_solver->get_best_clique().size()),
          m_initial_ub(best_covering_assignment.size()),
          m_best_solution(best_covering_assignment),
          m_current_lb(m_initial_lb) {}

  private:
    /**
     * The clique solver part.
     * The vertices in the 'clique' part and
     * the 'sat' part are different, relatively
     * independent collections.
     */
    IncrementalGurobiG2CliqueSolver* m_clique_solver;

    /**
     * Adjacency matrix, useful both for the clique
     * and the SAT part. Also keeps the full vertex set,
     * i.e., all (uncovered) interactions.
     */
    LazyG2AdjacencyMatrix* m_adjacency_matrix;

    /**
     * The best initial lower and upper bound.
     * These are used to determine the initially
     * available as well as the totally available
     * number of configurations.
     */
    std::size_t m_initial_lb, m_initial_ub;

    /**
     * The best solution found so far, initialized
     * to the initial upper bound solution.
     */
    std::vector<DynamicBitset> m_best_solution;

    /**
     * The current lower bound.
     */
    std::size_t m_current_lb;

    /**
     * The incremental SAT solver.
     */
    IncrementalSatSolver m_sat_solver;

    /**
     * Literals used in assumptions to disable
     * configurations that are not yet allowed.
     * Index 0 corresponds to the configuration with index m_initial_lb,
     * index 1 to m_initial_lb + 1, and so on.
     */
    std::vector<Lit> m_configuration_allowed;

    /**
     * The individual configurations that are present in the SAT model.
     * Configurations are added to the SAT model as soon as it becomes
     * clear that they are necessary (initially, MES size, later the best
     * known lower bound upon UNSAT results from the SAT solver.)
     */
    std::vector<std::vector<LitOrVal>> m_sat_configurations;

    /**
     * The partial configurations that are used and filled by
     * the incremental DSatur-like algorithm.
     */
    std::vector<SharedDBPropagator> m_dsatur_partial_configurations;

    /**
     * The vertices that should be explicitly handled,
     * in the order they became explicitly handled.
     * A superset of m_sat_explicit_vertices;
     * m_sat_explicit_vertices are extended to this
     * set when the SAT model is run.
     */
    std::vector<FullVertexIndex> m_explicit_vertices;
    DynamicBitset m_is_explicit;

    /**
     * The vertices explicitly handled and
     * inserted into the current SAT model.
     */
    std::vector<FullVertexIndex> m_sat_explicit_vertices;
    DynamicBitset m_is_sat_explicit;
};

} // namespace sammy

#endif
