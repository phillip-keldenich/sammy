#include <doctest/doctest.h>
#include <filesystem>
#include <sammy/barrage.h>
#include <sammy/cadical_solver.h>
#include <sammy/cmsat5_solver.h>
#include <sammy/incremental_sat_lns.h>
#include <sammy/io.h>
#include <sammy/kissat_solver.h>
#include <sammy/lingeling_solver.h>
#include <sammy/sat_dsatur.h>
#include <sammy/sat_lns.h>
#include <sammy/shared_db_propagator.h>

using namespace sammy;

#define INC_SOLVERS CMSAT5Solver, CadicalSolver, LingelingSolver

#define SAT_SOLVERS INC_SOLVERS, KissatSolver

#ifdef SAMMY_TEST_DATA

inline bool is_covered(Vertex v, const DynamicBitset& assignment) {
    auto v1 = lit::var(v.first), v2 = lit::var(v.second);
    bool neg1 = lit::negative(v.first), neg2 = lit::negative(v.second);
    return assignment[v1] != neg1 && assignment[v2] != neg2;
}

static std::unordered_set<std::string> KNOWN_INFEASIBLE = {
    "soletta-2015-07-02_18-48-59-subproblem1.json.xz"};

template <typename SubSolver>
static void test_solver_on(const std::string& subproblem_name) {
    std::filesystem::path test_dir(SAMMY_TEST_DATA);
    std::filesystem::path test_file = test_dir / subproblem_name;
    auto input_object = read_json_path(test_file);
    EventRecorder recorder;
    EventRecorder local_recorder;
    PortfolioSolver portfolio(input_object.at("portfolio_state"), &recorder);
    LNSSubproblem subproblem =
        LNSSubproblem::load(input_object.at("subproblem"));
    std::size_t expected_solution_size =
        subproblem.removed_configurations.size() - 1;
    std::vector<Vertex> universe_copy = subproblem.uncovered_universe;
    ClauseDB& clauses = portfolio.get_clauses();
    sammy::rng().seed(4213);
    SubSolver solver{&portfolio, std::move(subproblem),
                     SharedDBPropagator{&clauses}, &local_recorder, 42};
    auto res = solver.solve();
    REQUIRE(!!res);
    if (*res) {
        CHECK(!KNOWN_INFEASIBLE.count(subproblem_name));
        std::vector<DynamicBitset> solution = solver.get_solution();
        CHECK(solution.size() == expected_solution_size);
        for (const DynamicBitset& assignment : solution) {
            auto new_end = std::remove_if(
                universe_copy.begin(), universe_copy.end(),
                [&assignment](Vertex v) { return is_covered(v, assignment); });
            universe_copy.erase(new_end, universe_copy.end());
        }
        CHECK(universe_copy.empty());
    } else {
        CHECK(KNOWN_INFEASIBLE.count(subproblem_name));
    }
}

TEST_CASE_TEMPLATE(
    "FixedMESSatDSatur on soletta-2015-07-02_18-48-59-subproblem1",
    IncSatSolver, INC_SOLVERS) {
    using SubSolver = FixedMESSatDSaturSolver<IncSatSolver>;
    test_solver_on<SubSolver>(
        "soletta-2015-07-02_18-48-59-subproblem1.json.xz");
}

TEST_CASE_TEMPLATE("FixedMESSAT on soletta-2015-07-02_18-48-59-subproblem1",
                   SatSolver, SAT_SOLVERS) {
    using SubSolver = FixedMESSATImprovementSolver<SatSolver>;
    test_solver_on<SubSolver>(
        "soletta-2015-07-02_18-48-59-subproblem1.json.xz");
}

TEST_CASE_TEMPLATE(
    "FixedMESIncrementalSAT on soletta-2015-07-02_18-48-59-subproblem1",
    IncSatSolver, INC_SOLVERS) {
    using SubSolver = FixedMESIncrementalSATImprovementSolver<IncSatSolver>;
    test_solver_on<SubSolver>(
        "soletta-2015-07-02_18-48-59-subproblem1.json.xz");
}

#endif
