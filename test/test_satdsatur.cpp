#include "test_instances.h"
#include <doctest/doctest.h>
#include <sammy/initial_coloring_heuristic.h>
#include <sammy/learn_infeasibilities.h>
#include <sammy/lingeling_solver.h>
#include <sammy/satdsatur_colorer.h>
#include <sammy/verify.h>

using namespace sammy;

TEST_CASE("[SATDsatur] Test with known universe") {
    ClauseDB clauses{3, std::vector<ExternalClause>{{1, 2, 3}, {-1, -2, -3}}};
    std::vector<Vertex> clique{{lit::internalize(1), lit::internalize(2)},
                               {lit::internalize(1), lit::internalize(-2)},
                               {lit::internalize(-1), lit::internalize(2)},
                               {lit::internalize(-1), lit::internalize(-2)}};
    PairInfeasibilityMap inf_map{3};
    ThreadGroup<void> threads{0};
    ColoringHeuristicSolver solver{&clauses, 3, &threads, &inf_map};
    solver.set_quiet(true);
    solver.initialize_feasibilities();
    solver.color_lazy();
    solver.extract_feasibilities();
    learn_infeasibilities(clauses, &inf_map);
    auto whole_universe = inf_map.collect_vertices();
    UniverseSubgraph subgraph{&clauses, &threads, &inf_map, whole_universe};
    subgraph.extend_matrix_by_propagation();
    std::vector<std::size_t> clique_inds;
    std::transform(clique.begin(), clique.end(),
                   std::back_inserter(clique_inds),
                   [&](Vertex v) { return subgraph.vertex_index(v); });
    SATDSaturCoverer<LingelingSolver> coverer{&subgraph, std::move(clique_inds),
                                              10};
    coverer.run_coloring();
    CHECK(coverer.was_successful());
    CHECK(coverer.did_improve());
    PartialSolution psoln = coverer.get_partial_solution();
    CHECK(psoln.size() == 6);
    psoln.iterate_all_uncovered([](Lit, Lit) { CHECK(false); });
    for (std::size_t i = 0; i < psoln.size(); ++i) {
        CHECK(!solution_has_error(clauses, psoln.get_assignment(i)));
    }
}

TEST_CASE("[SATDSatur] Test with soletta instance") {
    VarCount varcount = soletta_2017_03_01_15_25_44_variables();
    ClauseDB clauses{varcount.n_all, soletta_2017_03_01_15_25_44_clauses()};
    PairInfeasibilityMap inf_map{varcount.n_concrete};
    ThreadGroup<void> threads{0};
    ColoringHeuristicSolver solver{&clauses, varcount.n_concrete, &threads,
                                   &inf_map};
    solver.set_quiet(true);
    solver.initialize_feasibilities();
    solver.color_lazy();
    solver.extract_feasibilities();
    learn_infeasibilities(clauses, &inf_map);
    PartialSolution best_solution = solver.get_partial_solution();
    auto whole_universe = inf_map.collect_vertices();
    CHECK(whole_universe.size() == 272307);
    CHECK(best_solution.size() == solver.all_classes().size());
    SharedDBPropagator propagator{&clauses};
    std::vector<std::pair<int, int>> known_mes = {
        {74, 131},   {131, 402}, {74, 143},   {143, 402}, {2, 436},
        {74, 436},   {2, 131},   {2, 143},    {402, 436}, {279, 402},
        {31, 436},   {303, 350}, {303, -359}, {74, 350},  {2, -359},
        {74, 351},   {2, 279},   {2, 351},    {31, 351},  {-359, 402},
        {131, 303},  {84, -389}, {-70, 381},  {-79, 345}, {-378, 451},
        {-142, 370}, {279, 303}};
    std::vector<Vertex> best_mes;
    std::transform(known_mes.begin(), known_mes.end(),
                   std::back_inserter(best_mes), [](std::pair<int, int> ve) {
                       return Vertex{lit::internalize(ve.first),
                                     lit::internalize(ve.second)};
                   });
    for (std::size_t i = 0; i < 3; ++i) {
        solver.reset_coloring();
        solver.color_lazy(best_mes);
        if (solver.all_classes().size() < best_solution.size()) {
            best_solution = solver.get_partial_solution();
        }
    }
    CHECK(best_solution.size() >= best_mes.size());
    for (std::size_t i : range(best_solution.size())) {
        CHECK(!solution_has_error(clauses, best_solution.get_assignment(i)));
    }
    CHECK(!mutually_exclusive_set_has_error(clauses, varcount.n_concrete,
                                            best_mes));
    std::vector<Vertex> uniquely_covered;
    best_solution.iterate_all_uniquely_covered(
        [&](Lit lmin, Lit lmax) { uniquely_covered.emplace_back(lmin, lmax); });
    uniquely_covered.insert(uniquely_covered.end(), best_mes.begin(),
                            best_mes.end());
    std::sort(uniquely_covered.begin(), uniquely_covered.end());
    uniquely_covered.erase(
        std::unique(uniquely_covered.begin(), uniquely_covered.end()),
        uniquely_covered.end());
    UniverseSubgraph subgraph{&clauses, &threads, &inf_map, uniquely_covered};
    SATDSaturCoverer<LingelingSolver> coverer{
        &subgraph, subgraph.to_indices(best_mes), best_solution.size()};
    coverer.run_coloring(false);
    CHECK(coverer.was_successful());
    CHECK(coverer.did_improve());
    std::vector<SharedDBPropagator> incomplete =
        coverer.get_incomplete_solution();
    solver.reset_coloring();
    for (SharedDBPropagator& p : incomplete) {
        solver.add_color_class(std::move(p));
    }
    solver.color_lazy(best_mes);
    CHECK(incomplete.size() >= best_mes.size());
    for (std::size_t i : range(solver.all_classes().size())) {
        CHECK(!solution_has_error(
            clauses, solver.get_partial_solution().get_assignment(i)));
    }
}
