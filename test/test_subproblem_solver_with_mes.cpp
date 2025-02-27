#include "test_instances.h"
#include <doctest/doctest.h>
#include <sammy/barrage.h>
#include <sammy/fast_clique.h>
#include <sammy/initial_coloring_heuristic.h>
#include <sammy/learn_infeasibilities.h>
#include <sammy/subproblem_solver_with_mes.h>
#include <sammy/thread_clauses.h>

using namespace sammy;

void check_mes(SharedDBPropagator& prop, const std::vector<Vertex>& v) {
    for (std::size_t i = 0, n = v.size(); i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            prop.reset_to_zero();
            CHECK(push_vertex(prop, v[i]) >= 0);
            CHECK(push_vertex(prop, v[j]) < 0);
        }
    }
}

#ifdef SAMMY_TEST_DATA

TEST_CASE("[SubproblemSolverWithMES] MES solver soletta_17_03_09") {
    std::filesystem::path path_data(SAMMY_TEST_DATA);
    path_data /= "soletta_17_03_09_initial.json";
    std::ifstream input(path_data, std::ios::in);
    OutputObject data = OutputObject::parse(input);
    auto initial_phase_data = import_initial_phase_result(data);
    ClauseDB& clauses = initial_phase_data.first;
    InitialPhaseResult& initial_result = initial_phase_data.second;
    auto ticket = publish_clauses(clauses);
    EventRecorder recorder;
    PortfolioSolver solver{ticket, &recorder, std::move(initial_result)};
    SharedDBPropagator prop(&clauses);
    CHECK(solver.get_universe_size() == 124634);
    solver.reduce_universe();
    CHECK(solver.implied_cache().have_reduced_universe());
    CHECK(solver.implied_cache().reduced_universe_size() == 20978);
    CHECK(solver.implied_cache().original_universe_size() == 124634);
    auto solution = solver.get_best_solution();
    CHECK(solution.size() == 53);
    CHECK(solver.get_best_mes_size() == 33);
    std::vector<Index> indices{0, 2, 3, 5, 6, 9, 11, 44};
    auto removed = solution.remove_assignments(indices);
    CHECK(removed.size() == 8);
    LNSSubproblem subproblem;
    solution.iterate_all_uncovered([&](Lit lmin, Lit lmax) {
        subproblem.uncovered_universe.emplace_back(lmin, lmax);
    });
    FastCliqueBuilder clique_builder{SharedDBPropagator(&clauses)};
    subproblem.mutually_exclusive_set =
        clique_builder.random_multistart_best_clique(
            10, subproblem.uncovered_universe);
    check_mes(prop, subproblem.mutually_exclusive_set);
    subproblem.removed_configurations = std::move(removed);
    subproblem.num_concrete = initial_result.inf_map.get_n_concrete();
    SubproblemMESSolver mes_solver{&solver,
                                   std::move(subproblem),
                                   SharedDBPropagator(&clauses),
                                   &recorder,
                                   42,
                                   1000};
    auto result = mes_solver.solve();
    CHECK(result);
    CHECK(*result);
    CHECK(mes_solver.mes_vertices().size() == 6);
    check_mes(prop, mes_solver.mes_vertices());
}

TEST_CASE("[SubproblemSolverWithMES] MES solver soletta_17_03_09 larger") {
    std::filesystem::path path_data(SAMMY_TEST_DATA);
    path_data /= "soletta_17_03_09_initial.json";
    std::ifstream input(path_data, std::ios::in);
    OutputObject data = OutputObject::parse(input);
    auto initial_phase_data = import_initial_phase_result(data);
    ClauseDB& clauses = initial_phase_data.first;
    InitialPhaseResult& initial_result = initial_phase_data.second;
    auto ticket = publish_clauses(clauses);
    EventRecorder recorder;
    PortfolioSolver solver{ticket, &recorder, std::move(initial_result)};
    SharedDBPropagator prop(&clauses);
    CHECK(solver.get_universe_size() == 124634);
    solver.reduce_universe();
    CHECK(solver.implied_cache().have_reduced_universe());
    CHECK(solver.implied_cache().reduced_universe_size() == 20978);
    CHECK(solver.implied_cache().original_universe_size() == 124634);
    auto solution = solver.get_best_solution();
    CHECK(solution.size() == 53);
    CHECK(solver.get_best_mes_size() == 33);
    std::vector<Index> indices{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
                               12, 13, 14, 15, 16, 17, 18, 31, 32, 33, 35, 44};
    auto removed = solution.remove_assignments(indices);
    CHECK(removed.size() == 24);
    LNSSubproblem subproblem;
    solution.iterate_all_uncovered([&](Lit lmin, Lit lmax) {
        subproblem.uncovered_universe.emplace_back(lmin, lmax);
    });
    FastCliqueBuilder clique_builder{SharedDBPropagator(&clauses)};
    subproblem.mutually_exclusive_set =
        clique_builder.random_multistart_best_clique(
            10, subproblem.uncovered_universe);
    check_mes(prop, subproblem.mutually_exclusive_set);
    subproblem.removed_configurations = std::move(removed);
    subproblem.num_concrete = initial_result.inf_map.get_n_concrete();
    SubproblemMESSolver mes_solver{&solver,
                                   std::move(subproblem),
                                   SharedDBPropagator(&clauses),
                                   &recorder,
                                   42,
                                   1000};
    auto result = mes_solver.solve();
    CHECK(result);
    CHECK(*result);
    CHECK(mes_solver.mes_vertices().size() == 18);
    check_mes(prop, mes_solver.mes_vertices());
}

#endif
