#include "test_instances.h"
#include <doctest/doctest.h>
#include <sammy/barrage.h>
#include <sammy/implied_vertices.h>
#include <sammy/initial_coloring_heuristic.h>
#include <sammy/learn_infeasibilities.h>
#include <sammy/simplification.h>

using namespace sammy;

TEST_CASE("[eliminate_implied_vertices] Test implied vertices (soletta "
          "unsimplified)") {
    VarCount vc = soletta_2017_03_01_15_25_44_variables();
    std::vector<ExternalClause> raw_clauses =
        soletta_2017_03_01_15_25_44_clauses();
    ClauseDB clauses(vc.n_all, raw_clauses);
    ThreadGroup<void> threads;
    PairInfeasibilityMap inf_map(vc.n_concrete);
    ColoringHeuristicSolver solver(&clauses, vc.n_concrete, &threads, &inf_map);
    solver.set_quiet(true);
    solver.initialize_feasibilities();
    solver.color_lazy();
    solver.extract_feasibilities();
    learn_infeasibilities(clauses, &inf_map);
    std::vector<Vertex> universe = inf_map.collect_vertices();
    CHECK(universe.size() == 272'307);
    std::vector<Vertex> filtered =
        eliminate_implied_vertices(universe, SharedDBPropagator{&clauses});
    CHECK(filtered.size() < universe.size());
    solver.reset_coloring();
    solver.add_vertices_to_queue(filtered.begin(), filtered.end());
    solver.color_vertices_in_queue();
    CHECK(solver.num_uncovered() == 0);
}

TEST_CASE(
    "[eliminate_implied_vertices] Test implied vertices (soletta simplified)") {
    VarCount vc = soletta_2017_03_01_15_25_44_variables();
    std::vector<ExternalClause> raw_clauses =
        soletta_2017_03_01_15_25_44_clauses();
    ClauseDB old_clauses(vc.n_all, raw_clauses);
    SimplifyDatastructure simplifier(old_clauses, vc.n_concrete);
    SimplificationStats stats;
    SimplifiedInstance simplified = run_simplifier(simplifier, stats);
    ThreadGroup<void> threads;
    PairInfeasibilityMap inf_map(simplified.num_concrete);
    ColoringHeuristicSolver solver(&simplified.formula, simplified.num_concrete,
                                   &threads, &inf_map);
    solver.set_quiet(true);
    solver.initialize_feasibilities();
    solver.color_lazy();
    solver.extract_feasibilities();
    learn_infeasibilities(simplified.formula, &inf_map);
    std::vector<Vertex> universe = inf_map.collect_vertices();
    std::vector<Vertex> filtered = eliminate_implied_vertices(
        universe, SharedDBPropagator{&simplified.formula});
    CHECK(filtered.size() < universe.size());
    solver.reset_coloring();
    solver.add_vertices_to_queue(filtered.begin(), filtered.end());
    solver.color_vertices_in_queue();
    CHECK(solver.num_uncovered() == 0);
}

#ifdef SAMMY_TEST_DATA

TEST_CASE("[eliminate_implied_vertices] Test implied correct & unchanged") {
    std::filesystem::path path_data(SAMMY_TEST_DATA);
    path_data /= "soletta_17_03_09_initial.json";
    std::ifstream input(path_data, std::ios::in);
    OutputObject data = OutputObject::parse(input);
    auto initial_phase_data = import_initial_phase_result(data);
    ClauseDB& clauses = initial_phase_data.first;
    InitialPhaseResult& initial_result = initial_phase_data.second;
    SharedDBPropagator propagator{&clauses};
    auto ticket = publish_clauses(clauses);
    EventRecorder recorder;
    PortfolioSolver solver{ticket, &recorder, std::move(initial_result)};
    CHECK(solver.get_universe_size() == 124634);
    solver.reduce_universe();
    CHECK(solver.implied_cache().have_reduced_universe());
    CHECK(solver.implied_cache().reduced_universe_size() == 20978);
    CHECK(solver.implied_cache().original_universe_size() == 124634);
}

#endif
