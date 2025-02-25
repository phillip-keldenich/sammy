#include <doctest/doctest.h>
#include <sammy/initial_coloring_heuristic.h>
#include <sammy/partial_solution.h>
#include <sammy/learn_infeasibilities.h>
#include "test_instances.h"

using namespace sammy;

TEST_CASE("[PartialSolution] Test with known universe") {
    ClauseDB clauses{3, std::vector<ExternalClause>{{1,2,3}, {-1,-2,-3}}};
    std::vector<Vertex> clique{
        {lit::internalize(1), lit::internalize(2)},
        {lit::internalize(1), lit::internalize(-2)},
        {lit::internalize(-1), lit::internalize(2)},
        {lit::internalize(-1), lit::internalize(-2)}
    };
    PairInfeasibilityMap inf_map{3};
    ThreadGroup<void> threads{0};
    ColoringHeuristicSolver solver{&clauses, 3, &threads, &inf_map};
    solver.set_quiet(true);
    solver.initialize_feasibilities();
    solver.color_lazy();
    solver.extract_feasibilities();
    learn_infeasibilities(clauses, &inf_map);
    auto whole_universe = inf_map.collect_vertices();
    CHECK(whole_universe.size() == 12);
    CHECK(inf_map[lit::internalize(1)][lit::internalize(-1)]);
    CHECK(inf_map[lit::internalize(2)][lit::internalize(-2)]);
    CHECK(inf_map[lit::internalize(3)][lit::internalize(-3)]);
    std::vector<std::vector<int>> fixed_solution{
        {-1, -2, 3}, {1, 2, -3}, {-1, 2, 3},
        {1, -2, -3}, {-1, 2, -3}, {1, -2, 3}
    };
    PartialSolution solution(3, &inf_map);
    for(const auto& c : fixed_solution) {
        DynamicBitset bs(3, false);
        for(std::size_t i : range(3)) {
            bs[i] = (c[i] > 0);
        }
        solution.add_assignment(bs);
    }
    solution.iterate_all_uncovered([] (Lit, Lit) {CHECK(false);});
    CHECK(solution.count_uncovered() == 0);
    std::vector<std::pair<int,int>> uniquely_covered_pairs;
    solution.iterate_all_uniquely_covered([&] (Lit lmin, Lit lmax) {
        uniquely_covered_pairs.emplace_back(lit::externalize(lmin), lit::externalize(lmax));
    });
    std::vector<std::pair<int,int>> uniques{
        {1, 2}, {1, 3}, {-1, -2}, {-1, -3}, {2, 3}, {-2, -3}};
    CHECK(uniques == uniquely_covered_pairs);
    auto removed = solution.remove_assignments({1,4});
    CHECK(removed.size() == 2);
    CHECK(removed[0][0]);
    CHECK(removed[0][1]);
    CHECK(!removed[0][2]);
    CHECK(!removed[1][0]);
    CHECK(removed[1][1]);
    CHECK(!removed[1][2]);
    std::vector<std::pair<int,int>> uncovered{
        {1, 2}, {-1, -3}, {2, -3}
    };
    solution.iterate_all_uncovered([&](Lit lmin, Lit lmax) {
        std::pair<int,int> p{lit::externalize(lmin), lit::externalize(lmax)};
        return std::find(uncovered.begin(), uncovered.end(), p) != uncovered.end();
    });
    CHECK(solution.count_uncovered() == 3);
    uniquely_covered_pairs.clear();
    solution.iterate_all_uniquely_covered([&](Lit lmin, Lit lmax) {
        uniquely_covered_pairs.emplace_back(lit::externalize(lmin), lit::externalize(lmax));
    }, true);
    CHECK(uniquely_covered_pairs == 
          std::vector<std::pair<int,int>>{{1,3}, {1,-3}, {-1, 2}, {-1,-2}, {2, 3}, {-2, -3}});
}


TEST_CASE("[PartialSolution] Test with soletta instance") {
    auto soletta_vars = soletta_2017_03_01_15_25_44_variables();
    ClauseDB clauses{soletta_vars.n_all, soletta_2017_03_01_15_25_44_clauses()};
    PairInfeasibilityMap inf_map{soletta_vars.n_concrete};
    ThreadGroup<void> threads;
    ColoringHeuristicSolver solver{&clauses, soletta_vars.n_concrete, &threads, &inf_map};
    solver.set_quiet(true);
    solver.initialize_feasibilities();
    solver.color_lazy();
    solver.extract_feasibilities();
    learn_infeasibilities(clauses, &inf_map);
    PartialSolution solution = solver.get_partial_solution();
    CHECK(solution.count_uncovered() == 0);
    std::vector<Vertex> uniquely_covered;
    solution.iterate_all_uniquely_covered([&] (Lit lmin, Lit lmax) {uniquely_covered.emplace_back(lmin,lmax);});
    std::vector<std::size_t> unique_count_per_class(solution.size(), 0);
    for(Vertex v : uniquely_covered) {
        unique_count_per_class[solver.find_covering_class(v.first, v.second)] += 1;
    }
    std::size_t s = solution.size();
    CHECK(s >= 20);
    auto removed = solution.remove_assignments({4});
    CHECK(removed.size() == 1);
    CHECK(solution.size() == s-1);
    std::size_t prev_uncovered = solution.count_uncovered();
    CHECK(unique_count_per_class[4] == prev_uncovered);
    uniquely_covered.clear();
    solution.iterate_all_uniquely_covered([&] (Lit lmin, Lit lmax) {uniquely_covered.emplace_back(lmin,lmax);}, true);
    unique_count_per_class.assign(solution.size(), 0);
    for(Vertex v : uniquely_covered) {
        unique_count_per_class[solution.find_covering_class(v.first, v.second)] += 1;
    }
    auto removed2 = solution.remove_assignments({0});
    CHECK(removed2.size() == 1);
    CHECK(solution.size() == s-2);
    CHECK(unique_count_per_class[0] + prev_uncovered == solution.count_uncovered());
    unique_count_per_class.assign(solution.size(), 0);
    std::vector<std::size_t> unique_count_per_class2(solution.size(), 0);
    uniquely_covered.clear();
    solution.iterate_all_uniquely_covered([&] (Lit lmin, Lit lmax) {uniquely_covered.emplace_back(lmin,lmax);}, true);
    for(Vertex v : uniquely_covered) {
        unique_count_per_class[solution.find_covering_class(v.first, v.second)] += 1;
    }
    solution.iterate_all_uniquely_covered_with_class([&] (Index c, Lit, Lit) {
        CHECK(c < solution.size());
        unique_count_per_class2[c] += 1;
    });
    CHECK(unique_count_per_class2 == unique_count_per_class);
}
