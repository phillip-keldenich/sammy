#include <doctest/doctest.h>
#include <sammy/simplification.h>
#include <sammy/randsat.h>
#include <sammy/verify.h>
#include "test_instances.h"

using namespace sammy;

TEST_CASE("[detect_failed_and_equal_literals] Failed and equal literal detection") {
    ClauseDB formula{7, std::vector<ExternalClause>{
        {1, 2}, {-1, -2}, {3, -4}, {-3, 5}, {-3, -5, 6},
        {-3, -5, -6, 7}, {-4, -7}
    }};
    SimplifyDatastructure simplify_ds{formula, 5};
    CHECK(detect_failed_and_equal_literals(simplify_ds));
    CHECK((simplify_ds.is_eliminated(0) || simplify_ds.is_eliminated(1)));
    CHECK((!simplify_ds.is_eliminated(0) || !simplify_ds.is_eliminated(1)));
    CHECK(simplify_ds.is_eliminated(3));
    CHECK(!simplify_ds.is_eliminated(2));
    CHECK(!simplify_ds.is_eliminated(4));
    CHECK(!simplify_ds.is_eliminated(5));
    CHECK(!simplify_ds.is_eliminated(6));
    SimplifiedInstance instance = simplify_ds.compress();
    CHECK(instance.num_concrete == 3);
    std::vector<bool> simp_sol{false, false, true, true, true};
    std::vector<bool> full_sol = simplify_ds.reconstruct_solution(instance, simp_sol);
    CHECK(full_sol.size() == 7);
    if(simplify_ds.is_eliminated(0)) {
        CHECK(full_sol[0] == true);
        CHECK(full_sol[1] == false);
    } else {
        CHECK(full_sol[0] == false);
        CHECK(full_sol[1] == true);
    }
    CHECK(full_sol[2] == false);
    CHECK(full_sol[3] == false);
    CHECK(full_sol[4] == true);
    CHECK(full_sol[5] == true);
    CHECK(full_sol[6] == true);
}

bool check_error(const std::optional<std::string>& error) {
    if(error) {
        std::cerr << "[ERROR] " << *error << std::endl;
        std::exit(1);
    }
    return !error;
}

TEST_CASE("[detect_failed_and_equal_literals] Level 0 literals & subsumption, soletta") {
    std::vector<ExternalClause> clauses = soletta_2017_03_01_15_25_44_clauses();
    Var n_all = 461, n_concrete = 460;
    ClauseDB clause_db{n_all, clauses};
    SimplifyDatastructure simplify_ds{clause_db, n_concrete};
    CHECK(fix_level0_literals(simplify_ds));
    // with only level-0-literal elimination,
    // the following reduced clause (original: (-296 -267 -96 -49 71 106))
    // should be included in the simplified formula
    std::unordered_set<Lit> reduced_original_clause{
        lit::internalize(-296), lit::internalize(-267),
        lit::internalize(-49), lit::internalize(71),
        lit::internalize(106)};
    CHECK(simplify_ds.has_clause_subsuming(reduced_original_clause));
    eliminate_subsumed(simplify_ds.clauses(), n_all);
    CHECK(simplify_ds.has_clause_subsuming(reduced_original_clause));
    SimplifiedInstance simplified = simplify_ds.compress();
    for(int i = 0; i < 100; ++i) {
        auto sol_simplified1 = randsolve(simplified.formula);
        auto sol_reconstructed = simplify_ds.reconstruct_solution(simplified, sol_simplified1);
        CHECK(check_error(solution_has_error(simplified.formula, sol_simplified1)));
        CHECK(check_error(solution_has_error(clause_db, sol_reconstructed)));
    }
}

TEST_CASE("[detect_failed_and_equal_literals] Failed and equal literals, soletta") {
    std::vector<ExternalClause> clauses = soletta_2017_03_01_15_25_44_clauses();
    Var n_all = 461, n_concrete = 460;
    ClauseDB clause_db{n_all, clauses};
    SimplifyDatastructure simplify_ds{clause_db, n_concrete};
    CHECK(detect_failed_and_equal_literals(simplify_ds));
    eliminate_subsumed(simplify_ds.clauses(), n_all);
    SimplifiedInstance simplified1 = simplify_ds.compress();
    for(int i = 0; i < 5; ++i) {
        auto sol_simplified1 = randsolve(simplified1.formula);
        auto sol_reconstructed = simplify_ds.reconstruct_solution(simplified1, sol_simplified1);
        CHECK(solution_has_error(simplified1.formula, sol_simplified1) == std::nullopt);
        CHECK(solution_has_error(clause_db, sol_reconstructed) == std::nullopt);   
    }
    detect_failed_and_equal_literals(simplify_ds);
    eliminate_subsumed(simplify_ds.clauses(), n_all);
    SimplifiedInstance simplified2 = simplify_ds.compress();
    for(int i = 0; i < 5; ++i) {
        auto sol_simplified2 = randsolve(simplified2.formula);
        auto sol_reconstructed2 = simplify_ds.reconstruct_solution(simplified2, sol_simplified2);
        CHECK(solution_has_error(simplified2.formula, sol_simplified2) == std::nullopt);
        CHECK(solution_has_error(clause_db, sol_reconstructed2) == std::nullopt);   
    }
}
