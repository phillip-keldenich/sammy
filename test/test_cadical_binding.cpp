#include <sammy/cadical_solver.h>
#include <doctest/doctest.h>

using namespace sammy;

TEST_CASE("[CaDiCaL] Test bindings") {
    CadicalSolver solver;
    solver.add_short_clause(1, 2);
    solver.add_short_clause(1, -2);
    solver.add_short_clause(-1, 2);
    auto res = solver.solve();
    REQUIRE(res.has_value());
    REQUIRE(*res);
    auto model = solver.get_model();
    REQUIRE(model[1]);
    REQUIRE(model[2]);
    solver.add_short_clause(-1, -2);
    res = solver.solve();
    REQUIRE(res.has_value());
    REQUIRE(!*res);
}

TEST_CASE("[CaDiCaL] Test bindings with assertions") {
    CadicalSolver solver;
    solver.add_short_clause(1, 2);
    solver.add_short_clause(1, -2);
    solver.add_short_clause(-1, 2);
    auto res = solver.solve({-1}, 1.0);
    REQUIRE(res.has_value());
    REQUIRE(!*res);
}

TEST_CASE("[CaDiCaL] Test timeout") {
    using Lit = CadicalSolver::Lit;
    CadicalSolver solver;
    std::vector<std::vector<Lit>> pidgeon_holes(32);
    for(int pidgeon = 1; pidgeon <= 31; ++pidgeon) {
        for(int hole = 1; hole < 31; ++hole) {
            pidgeon_holes[pidgeon].push_back(solver.new_var());
        }
        for(Lit l : pidgeon_holes[pidgeon]) {
            solver.add_literal(l);
        }
        solver.finish_clause();
    }
    for(int hole = 1; hole < 31; ++hole) {
        for(int p1 = 1; p1 < 31; ++p1) {
            for(int p2 = p1 + 1; p2 <= 31; ++p2) {
                solver.add_short_clause(-pidgeon_holes[p1][hole - 1], -pidgeon_holes[p2][hole - 1]);
            }
        }
    }
    auto res = solver.solve({}, 2.0);
    if(res.has_value()) {
        REQUIRE(!*res);
    }
}
