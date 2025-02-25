#include <doctest/doctest.h>
#include <sammy/coloring.h>

using namespace sammy;

TEST_CASE("[SATKColoring] SAT coloring validity") {
    EventRecorder recorder;
    ClauseDB clauses{3, std::vector<ExternalClause>{{1,2,3}, {-1,-2,-3}}};
    PairInfeasibilityMap inf_map{clauses, 3};
    ThreadGroup<void> threads;
    ColoringHeuristicSolver solver{&clauses, 3, &threads, &inf_map};
    solver.set_quiet(true);
    solver.initialize_feasibilities();
    solver.color_lazy();
    solver.extract_feasibilities();
    std::vector<Vertex> vertices;
    for(Var i = 0; i < 2; ++i) {
        for(Var j = i + 1; j < 3; ++j) {
            Vertex vs[4] = {{lit::positive_lit(i), lit::positive_lit(j)},
                            {lit::negative_lit(i), lit::positive_lit(j)},
                            {lit::positive_lit(i), lit::negative_lit(j)},
                            {lit::negative_lit(i), lit::negative_lit(j)}};
            for(Vertex v : vs) {
                CHECK(solver.is_pair_known_feasible(v.first, v.second));
            }
            vertices.insert(vertices.end(), std::begin(vs), std::end(vs));
        }
    }
    UniverseSubgraph subgraph{&clauses, &threads, &inf_map, vertices};
    subgraph.extend_matrix_by_propagation();
    std::vector<std::size_t> clique{
        subgraph.vertex_index(Vertex{lit::positive_lit(0), lit::positive_lit(1)}),
        subgraph.vertex_index(Vertex{lit::negative_lit(0), lit::positive_lit(1)}),
        subgraph.vertex_index(Vertex{lit::positive_lit(0), lit::negative_lit(1)}),
        subgraph.vertex_index(Vertex{lit::negative_lit(0), lit::negative_lit(1)})
    };
    std::size_t coloring_number = 6;
    SATKColoringSolver csolver{&subgraph, &recorder, clique, coloring_number};
    auto optres = csolver.solve();
    CHECK(!!optres);
    CHECK(*optres);
    auto coloring = csolver.get_coloring();
    CHECK(coloring.size() == subgraph.n());
    CHECK(*std::max_element(coloring.begin(), coloring.end()) == coloring_number - 1);
    for(std::size_t i = 0, n = subgraph.n(); i < n; ++i) {
        for(std::size_t j : subgraph.matrix_row(i).ones()) {
            CHECK(i != j);
            CHECK(coloring[i] != coloring[j]);
        }
    }
}
