#include "test_instances.h"
#include <doctest/doctest.h>
#include <sammy/clique_sat_dsatur.h>
#include <sammy/implied_vertices.h>
#include <sammy/initial_coloring_heuristic.h>
#include <sammy/learn_infeasibilities.h>
#include <sammy/lingeling_solver.h>
#include <sammy/run_initial.h>
#include <sammy/simplification.h>
#include <sammy/verify.h>

using namespace sammy;

TEST_CASE("[CliqueSATDSatur] soletta optimal solve (initialized with "
          "non-optimal clique)") {
    VarCount count = soletta_2017_03_01_15_25_44_variables();
    ClauseDB clauses{count.n_all, soletta_2017_03_01_15_25_44_clauses()};
    SimplifyDatastructure simplifier{clauses, count.n_concrete};
    SimplificationStats simp_stats;
    SimplifiedInstance simplified = run_simplifier(simplifier, simp_stats);
    PairInfeasibilityMap inf_map{simplified.num_concrete};
    ThreadGroup<void> threads;
    ColoringHeuristicSolver solver{&simplified.formula, simplified.num_concrete,
                                   &threads, &inf_map};
    solver.initialize_feasibilities();
    solver.set_quiet(true);
    solver.color_lazy();
    solver.extract_feasibilities();
    learn_infeasibilities(simplified.formula, &inf_map);
    PartialSolution best_solution = solver.get_partial_solution();
    std::vector<Vertex> spawners = solver.class_spawners();
    solver.reset_coloring();
    solver.color_lazy(spawners);
    if (solver.all_classes().size() <= best_solution.size()) {
        best_solution = solver.get_partial_solution();
        spawners = solver.class_spawners();
    }
    ParallelFastCliqueBuilder clq_builder{
        SharedDBPropagator(&simplified.formula), &threads};
    std::vector<Vertex> best_clique;
    for (std::size_t i = 0; i < 5; ++i) {
        auto clique = clq_builder.random_multistart_best_clique(20, spawners);
        if (clique.size() > best_clique.size()) {
            best_clique = clique;
        }
        solver.reset_coloring();
        solver.color_lazy(clique);
        if (solver.all_classes().size() <= best_solution.size()) {
            best_solution = solver.get_partial_solution();
        }
        spawners = solver.class_spawners();
    }
    CHECK(!mutually_exclusive_set_has_error(
        simplified.formula, simplified.num_concrete, best_clique));
    for (const auto& assignment : best_solution.assignments()) {
        CHECK(!solution_has_error(simplified.formula, assignment));
    }
    simplified = remove_subsumed(simplified);
    CHECK(best_clique.size() <= 37);
    std::vector<Vertex> unfiltered_universe = inf_map.collect_vertices();
    CHECK(unfiltered_universe.size() == 124'383);
    std::vector<Vertex> filtered_universe = eliminate_implied_vertices(
        unfiltered_universe, SharedDBPropagator{&simplified.formula});
    filtered_universe.insert(filtered_universe.end(), best_clique.begin(),
                             best_clique.end());
    std::sort(filtered_universe.begin(), filtered_universe.end());
    filtered_universe.erase(
        std::unique(filtered_universe.begin(), filtered_universe.end()),
        filtered_universe.end());
    CliqueSatDSaturSolver<LingelingSolver> csds{
        filtered_universe,          &inf_map,
        simplified.formula,         best_clique,
        best_clique.size(),         best_clique.size(),
        best_solution.assignments()};
    using SR = CliqueSatDSaturSolver<LingelingSolver>::SolveResult;
    SR result = csds.solve();
    CHECK(result != SR::ABORTED);
    if (result != SR::SOLUTION_WAS_OPTIMAL) {
        CHECK(result == SR::IMPROVED_SOLUTION);
        best_solution = csds.get_partial_solution();
    }
    CHECK(best_solution.size() == 37);
    for (const auto& assignment : best_solution.assignments()) {
        CHECK(!solution_has_error(simplified.formula, assignment));
    }
    best_solution.iterate_all_uncovered([&](Lit lmin, Lit lmax) {
        CHECK(lmin == lmax);
        CHECK(false);
    });
}
