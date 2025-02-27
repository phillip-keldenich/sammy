#include <doctest/doctest.h>
#include <sammy/shared_db_propagator.h>

using namespace sammy;

TEST_CASE("[SharedDBPropagator] Simple propagator test") {
    std::vector<std::vector<int>> clauses_sat;
    std::vector<std::vector<int>> clauses_unsat;
    for (int diff = 1; diff <= 3; ++diff) {
        for (int start = 1; start + 2 * diff <= 8; ++start) {
            clauses_sat.emplace_back(std::initializer_list<int>{
                start, start + diff, start + 2 * diff});
            clauses_sat.emplace_back(std::initializer_list<int>{
                -start, -start - diff, -start - 2 * diff});
        }
    }
    for (int diff = 1; diff <= 4; ++diff) {
        for (int start = 1; start + 2 * diff <= 9; ++start) {
            clauses_unsat.emplace_back(std::initializer_list<int>{
                start, start + diff, start + 2 * diff});
            clauses_unsat.emplace_back(std::initializer_list<int>{
                -start, -start - diff, -start - 2 * diff});
        }
    }

    ClauseDB db_sat{8, clauses_sat};
    ClauseDB db_unsat{9, clauses_unsat};

    SUBCASE("sat propagation") {
        SharedDBPropagator propagator{&db_sat};
        CHECK(propagator.get_trail().empty());
        CHECK(propagator.push_level(lit::internalize(1)));
        CHECK(!propagator.is_conflicting());
        CHECK(propagator.get_trail().size() == 1);
        CHECK(propagator.get_reason(lit::internalize(1)).reason_length == 0);
        CHECK(propagator.push_level(lit::internalize(2)));
        CHECK(!propagator.is_conflicting());
        CHECK(propagator.get_trail().size() == 3);
        CHECK(lit::externalize(propagator.get_trail()[2]) == -3);
        Reason r3 = propagator.get_reason(lit::internalize(-3));
        CHECK(r3.reason_length == 3);
        auto r3lits = r3.lits(db_sat);
        std::vector<Lit> r3lits_expected{
            lit::internalize(-3), lit::internalize(-2), lit::internalize(-1)};
        CHECK(
            std::equal(r3lits.begin(), r3lits.end(), r3lits_expected.begin()));
        auto supporting = propagator.decisions_leading_to(lit::internalize(-3));
        CHECK(supporting.size() == 2);
        CHECK(supporting[0].second != supporting[1].second);
        CHECK(supporting[0].first <= 2);
        CHECK(supporting[0].first >= 1);
        CHECK(supporting[1].first <= 2);
        CHECK(supporting[1].first >= 1);
        CHECK(lit::externalize(supporting[0].second) >= 1);
        CHECK(lit::externalize(supporting[0].second) <= 2);
        CHECK(lit::externalize(supporting[1].second) >= 1);
        CHECK(lit::externalize(supporting[1].second) <= 2);
        CHECK(propagator.is_open(lit::internalize(4)));
        CHECK(!propagator.push_level(lit::internalize(4)));
        CHECK(propagator.is_conflicting());
        propagator.pop_level();
        CHECK(!propagator.is_conflicting());
        CHECK(propagator.push_level(lit::internalize(5)));
        CHECK(propagator.get_trail().size() == 5);
        CHECK(propagator.is_open(lit::internalize(4)));
        CHECK(propagator.is_open(lit::internalize(6)));
        CHECK(propagator.is_open(lit::internalize(7)));
        supporting = propagator.decisions_leading_to(lit::internalize(8));
        CHECK(supporting.size() == 2);
        CHECK(propagator.push_level(lit::internalize(-4)));
        CHECK(!propagator.is_conflicting());
        CHECK(propagator.get_trail().size() == 8);
        supporting = propagator.decisions_leading_to(lit::internalize(-7));
        CHECK(supporting.size() == 3);
        std::vector<Lit> decisions{supporting[0].second, supporting[1].second,
                                   supporting[2].second};
        std::vector<std::int32_t> levels{
            supporting[0].first, supporting[1].first, supporting[2].first};
        std::sort(decisions.begin(), decisions.end());
        std::sort(levels.begin(), levels.end());
        CHECK(decisions == std::vector<Lit>{lit::internalize(2),
                                            lit::internalize(-4),
                                            lit::internalize(5)});
        CHECK(levels == std::vector<std::int32_t>{2, 3, 4});
    }

    SUBCASE("unsat propagation") {
        SharedDBPropagator prop{&db_unsat};
        bool have_result = false;
        bool is_unsat = false;
        while (!have_result) {
            bool all_set = true;
            for (Var i = 1; i <= 9; ++i) {
                Lit x = lit::internalize(i);
                if (prop.is_open(x)) {
                    all_set = false;
                    if (!prop.push_level(x)) {
                        if (!prop.resolve_conflicts()) {
                            is_unsat = true;
                            have_result = true;
                        }
                    }
                    break;
                }
            }
            if (all_set)
                have_result = true;
        }
        CHECK(is_unsat);
    }
}
