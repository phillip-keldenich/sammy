#include <doctest/doctest.h>
#include <sammy/simplification.h>

using namespace sammy;

TEST_CASE("[eliminate_subsumed] Subsumption test 1") {
    auto plit = lit::positive_lit;
    auto nlit = lit::negative_lit;

    std::vector<CVec> subsumption_test1{{plit(0), plit(1), nlit(2)}, // 0
                                        {plit(1), nlit(2), plit(3)}, // 1
                                        {plit(1), nlit(2), plit(3)}, // 2 subsumed by 1
                                        {plit(1), nlit(2), plit(3)}, // 3 subsumed by 1
                                        {plit(0), plit(1), nlit(3)}, // 4
                                        {plit(0), plit(1), nlit(2), nlit(3)}, // 5 subsumed by 0
                                        {nlit(0), nlit(1), plit(2), nlit(3)}}; // 6

    eliminate_subsumed(subsumption_test1, 4);
    CHECK(subsumption_test1.size() == 4);
    CHECK(subsumption_test1[0] == (CVec{plit(0), plit(1), nlit(2)}));
    CHECK(subsumption_test1[1] == (CVec{plit(1), nlit(2), plit(3)}));
    CHECK(subsumption_test1[2] == (CVec{plit(0), plit(1), nlit(3)}));
    CHECK(subsumption_test1[3] == (CVec{nlit(0), nlit(1), plit(2), nlit(3)}));
}


// shorthand for positive/negative literals
static auto p = lit::positive_lit;
static auto n = lit::negative_lit;

TEST_CASE("[eliminate_subsumed] No subsumption at all") {
    // none of these clauses subsume any other
    std::vector<CVec> clauses{
        {p(0), p(1)},
        {p(1), p(2)},
        {p(2), p(3)}
    };

    eliminate_subsumed(clauses, /*n_all=*/4);
    CHECK(clauses.size() == 3);
    CHECK(clauses[0] == (CVec{p(0), p(1)}));
    CHECK(clauses[1] == (CVec{p(1), p(2)}));
    CHECK(clauses[2] == (CVec{p(2), p(3)}));
}

TEST_CASE("[eliminate_subsumed] Duplicate clauses collapse to one") {
    // clause 0 and clause 1 are identical, so one will subsume the other
    std::vector<CVec> clauses{
        {p(0), p(1)},  // 0
        {p(0), p(1)},  // 1  (duplicate)
        {p(1), p(2)}   // 2
    };

    eliminate_subsumed(clauses, /*n_all=*/3);
    // one of the two {0,1} should survive, plus the {1,2}
    CHECK(clauses.size() == 2);
    // order is preserved among survivors
    CHECK(clauses[0] == (CVec{p(0), p(1)}));
    CHECK(clauses[1] == (CVec{p(1), p(2)}));
}

TEST_CASE("[eliminate_subsumed] Throws on initial empty clause") {
    static auto p = lit::positive_lit;
    std::vector<CVec> clauses{
        {},           // 0  <-- now considered invalid
        {p(0)},       // 1
        {p(0), p(1)}  // 2
    };

    // p_init_watches should throw before any subsumption work happens
    CHECK_THROWS_AS(eliminate_subsumed(clauses, /*n_all=*/2),
                    std::logic_error);
}


TEST_CASE("[eliminate_subsumed] Transitive cascade of subsumption") {
    // A={0}, B={0,1}, C={0,1,2}.  A subsumes both B and C in one pass.
    std::vector<CVec> clauses{
        {p(0)},              // 0
        {p(0), p(1)},        // 1
        {p(0), p(1), p(2)}   // 2
    };

    eliminate_subsumed(clauses, /*n_all=*/3);
    // only the smallest (unit) clause should remain
    CHECK(clauses.size() == 1);
    CHECK(clauses[0] == (CVec{p(0)}));
}

TEST_CASE("[eliminate_subsumed] 10 times the same clause") {
    // If we have 10 copies of the same clause, we should end up with one
    // clause after subsumption.
    std::vector<CVec> clauses{
        {p(0), p(1)},  // 0
        {p(0), p(1)},  // 1
        {p(0), p(1)},  // 2
        {p(0), p(1)},  // 3
        {p(0), p(1)},  // 4
        {p(0), p(1)},  // 5
        {p(0), p(1)},  // 6
        {p(0), p(1)},  // 7
        {p(0), p(1)},  // 8
        {p(0), p(1)}   // 9
    };
    eliminate_subsumed(clauses, /*n_all=*/2);
    // only one copy of the clause should remain
    CHECK(clauses.size() == 1);
    CHECK(clauses[0] == (CVec{p(0), p(1)}));
}
