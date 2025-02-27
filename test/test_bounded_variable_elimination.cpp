#include <doctest/doctest.h>
#include <sammy/simplification.h>
#include <sammy/verify.h>

using namespace sammy;

TEST_CASE("[bounded_variable_elimination] Pure variable elimination case 1") {
    ClauseDB formula{9, std::vector<ExternalClause>{{1, 2, 3},
                                                    {-1, 2, 3},
                                                    {3, 4, 5, -6},
                                                    {6, 7, 8},
                                                    {4, 5, 8},
                                                    {2, 5, 9},
                                                    {2, 3, -4},
                                                    {4, 6, -8},
                                                    {2, 4, 5, -8},
                                                    {2, 5, 8, -9},
                                                    {-3, 9}}};
    SimplifyDatastructure simplifier{formula, 2};
    CHECK(bounded_variable_elimination(simplifier, 0));
    // 1, 2 are concrete
    CHECK(!simplifier.is_eliminated(0));
    CHECK(simplifier.is_concrete(0));
    CHECK(simplifier.is_concrete(1));
    CHECK(!simplifier.is_eliminated(1));
    CHECK(!simplifier.is_concrete(2));
    // 5 and 7 are pure
    CHECK(simplifier.is_eliminated(4));
    CHECK(simplifier.is_eliminated(6));
    // after eliminating 5 and 7, 8 and 9 are pure
    CHECK(simplifier.is_eliminated(7));
    CHECK(simplifier.is_eliminated(8));
    // now, 4 is pure as well
    CHECK(simplifier.is_eliminated(3));
    // now, 3 and 6 are pure as well
    CHECK(!simplifier.is_concrete(2));
    CHECK(simplifier.is_eliminated(2));
    CHECK(simplifier.is_eliminated(5));
    // check that everything was done as pure literal elimination
    for (const SCVec& c : simplifier.get_reconstruction_stack()) {
        CHECK(c.size() == 1);
    }
}

TEST_CASE(
    "[bounded_variable_elimination] Bounded variable elimination case 1") {
    ClauseDB formula{4, std::vector<ExternalClause>{{1, 2, 3, 4},
                                                    {1, -2, -3, -4},
                                                    {-1, 2, -3, -4},
                                                    {-1, -2, -3, -4},
                                                    {-1, -2, 3, -4},
                                                    {1, 2, 3, -4},
                                                    {1, -2, 3},
                                                    {1, 2, -3},
                                                    {1, -2, -3}}};
    SimplifyDatastructure simplifier{formula, 1};
    CHECK(bounded_variable_elimination(simplifier, 0));
    CHECK(simplifier.is_eliminated(3));
    CHECK(simplifier.is_eliminated(2));
    CHECK(simplifier.is_eliminated(1));
    CHECK(simplifier.get_reconstruction_stack().size() >= 5);
    CHECK(simplifier.get_reconstruction_stack()[0].size() == 4);
    CHECK(simplifier.get_reconstruction_stack()[0][0] == lit::positive_lit(3));
    CHECK(simplifier.get_reconstruction_stack()[1][0] == lit::negative_lit(3));
    CHECK(simplifier.get_reconstruction_stack()[2][0] == lit::negative_lit(3));
    CHECK(simplifier.get_reconstruction_stack()[3][0] == lit::negative_lit(3));
    CHECK(simplifier.get_reconstruction_stack()[4][0] == lit::negative_lit(3));
    CHECK(!simplifier.clauses().empty());
    CHECK(simplifier.clauses().size() == 1);
    CHECK(simplifier.clauses()[0].size() == 1);
    CHECK(simplifier.clauses()[0][0] == lit::positive_lit(0));
    SimplifiedInstance simp1 = simplifier.compress();
    CHECK(simp1.num_concrete == 1);
    CHECK(simp1.formula.num_vars() == 1);
    CHECK(simp1.formula.num_unaries() == 1);
    CHECK(simp1.new_to_old.size() == 1);
    CHECK(simp1.new_to_old[0] == 0);
    CHECK(!verify_solution(simp1.formula, std::vector<bool>(1, false)));
    CHECK(verify_solution(simp1.formula, std::vector<bool>(1, true)));
    std::vector<bool> reconstructed =
        simplifier.reconstruct_solution(simp1, std::vector<bool>(1, true));
    CHECK(verify_solution(formula, reconstructed));
    CHECK(detect_failed_and_equal_literals(simplifier));
    CHECK(simplifier.is_eliminated(0));
    SimplifiedInstance simp2 = simplifier.compress();
    reconstructed = simplifier.reconstruct_solution(simp2, std::vector<bool>{});
    CHECK(verify_solution(formula, reconstructed));
}
