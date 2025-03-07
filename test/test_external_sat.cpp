#include <doctest/doctest.h>
#include <sammy/external_sat_solver.h>

using namespace sammy;

template <ExternalSolverType EST> struct ESTW {
    static constexpr ExternalSolverType value = EST;
};

TEST_CASE_TEMPLATE("[ExternalSAT] Test external SAT solver", ESTWT,
                   ESTW<ExternalSolverType::KISSAT>,
                   ESTW<ExternalSolverType::CADICAL>,
                   ESTW<ExternalSolverType::LINGELING>,
                   ESTW<ExternalSolverType::CRYPTOMINISAT>) {
    ExternalNonIncrementalSAT<ESTWT::value> solver;
    auto v1 = solver.new_var();
    auto v2 = solver.new_var();
    auto v3 = solver.new_var();
    solver.add_short_clause(v1, -v2, -v3);
    solver.add_short_clause(-v1, -v2, -v3);
    solver.add_short_clause(-v2, v3);
    solver.add_short_clause(v2, v3);
    solver.add_short_clause(v1, v2);
    solver.add_short_clause(v1, v3);
    auto res = solver.solve();
    auto model = solver.get_model();
    REQUIRE(res);
    REQUIRE(*res);
    REQUIRE(model.raw().size() == 3);
    REQUIRE(model[v1]);
    REQUIRE(!model[-v1]);
    REQUIRE(!model[v2]);
    REQUIRE(model[-v2]);
    REQUIRE(model[v3]);
    REQUIRE(!model[-v3]);
}
