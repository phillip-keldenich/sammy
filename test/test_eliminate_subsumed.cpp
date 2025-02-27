#include <doctest/doctest.h>
#include <sammy/simplification.h>

using namespace sammy;

TEST_CASE("[eliminate_subsumed] Subsumption test 1") {
    auto plit = lit::positive_lit;
    auto nlit = lit::negative_lit;

    std::vector<CVec> subsumption_test1{{plit(0), plit(1), nlit(2)},
                                        {plit(1), nlit(2), plit(3)},
                                        {plit(1), nlit(2), plit(3)},
                                        {plit(1), nlit(2), plit(3)},
                                        {plit(0), plit(1), nlit(3)},
                                        {plit(0), plit(1), nlit(2), nlit(3)},
                                        {nlit(0), nlit(1), plit(2), nlit(3)}};

    eliminate_subsumed(subsumption_test1, 4);
    CHECK(subsumption_test1.size() == 4);
    CHECK(subsumption_test1[0] == (CVec{plit(0), plit(1), nlit(2)}));
    CHECK(subsumption_test1[1] == (CVec{plit(1), nlit(2), plit(3)}));
    CHECK(subsumption_test1[2] == (CVec{plit(0), plit(1), nlit(3)}));
    CHECK(subsumption_test1[3] == (CVec{nlit(0), nlit(1), plit(2), nlit(3)}));
}
