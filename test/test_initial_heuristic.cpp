#include <doctest/doctest.h>
#include <sammy/initial_coloring_heuristic.h>

using namespace sammy;

TEST_CASE("[Initial Coloring Heuristic] Simple coloring with known universe") {
    Var n_all = 4, n_concrete = 3;
    std::vector<ExternalClause> simple_formula{{4}};
    ClauseDB clauses{n_all, simple_formula};
    PairInfeasibilityMap infeasibility_map{n_concrete};
    ThreadGroup<void> thread_pool;
    ColoringHeuristicSolver solver{&clauses, n_concrete, &thread_pool,
                                   &infeasibility_map};
    solver.set_quiet(true);
    solver.initialize_feasibilities();
    solver.color_lazy();
    solver.extract_feasibilities();
    CHECK(infeasibility_map.is_definitely_feasible(0, 2));
    CHECK(infeasibility_map.is_definitely_feasible(0, 3));
    CHECK(infeasibility_map.is_definitely_feasible(0, 4));
    CHECK(infeasibility_map.is_definitely_feasible(0, 5));
    CHECK(infeasibility_map.is_definitely_feasible(1, 2));
    CHECK(infeasibility_map.is_definitely_feasible(1, 3));
    CHECK(infeasibility_map.is_definitely_feasible(1, 4));
    CHECK(infeasibility_map.is_definitely_feasible(1, 5));
    CHECK(infeasibility_map.is_definitely_feasible(2, 4));
    CHECK(infeasibility_map.is_definitely_feasible(2, 5));
    CHECK(infeasibility_map.is_definitely_feasible(3, 4));
    CHECK(infeasibility_map.is_definitely_feasible(3, 5));
    CHECK(solver.all_classes().size() >= 4);
    std::vector<Vertex> spawners = solver.class_spawners();
    solver.reset_coloring();
    solver.color_lazy(spawners);
    CHECK(solver.all_classes().size() >= 4);
    CHECK(solver.all_classes().size() <= 6);
}
