#include <sammy/cmsat5_solver.h>
#include <doctest/doctest.h>

using namespace sammy;

using CLit = CMSAT5Solver::Lit;

TEST_CASE("[CMSat5] Test bindings") {
    CMSAT5Solver solver;
    CLit l1 = solver.new_var();
    CLit l2 = solver.new_var();
    solver.add_short_clause(l1, l2);
    solver.add_short_clause(l1, -l2);
    solver.add_short_clause(-l1, l2);
    auto res = solver.solve();
    REQUIRE(res.has_value());
    REQUIRE(*res);
    auto model = solver.get_model();
    REQUIRE(model[l1]);
    REQUIRE(model[l2]);
    solver.add_short_clause(-l1, -l2);
    res = solver.solve();
    REQUIRE(res.has_value());
    REQUIRE(!*res);
}

TEST_CASE("[CMSat5] Test bindings with assertions") {
    CMSAT5Solver solver;
    CLit l1 = solver.new_var();
    CLit l2 = solver.new_var();
    solver.add_short_clause(l1, l2);
    solver.add_short_clause(l1, -l2);
    solver.add_short_clause(-l1, l2);
    auto res = solver.solve({-l1}, 1.0);
    REQUIRE(res.has_value());
    REQUIRE(!*res);
}

TEST_CASE("[CMSat5] Test timeout") {
    CMSAT5Solver solver;
    std::vector<std::vector<CLit>> pidgeon_holes(32);
    for(int pidgeon = 1; pidgeon <= 31; ++pidgeon) {
        for(int hole = 1; hole < 31; ++hole) {
            pidgeon_holes[pidgeon].push_back(solver.new_var());
        }
        for(CLit l : pidgeon_holes[pidgeon]) {
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
