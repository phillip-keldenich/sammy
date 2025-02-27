#include <algorithm>
#include <doctest/doctest.h>
#include <sammy/cuda_iteration.h>
#include <sammy/initial_coloring_heuristic.h>

using namespace sammy;

static std::vector<Vertex> filter_universe_if_covered(
    const std::vector<Vertex>& universe,
    const std::vector<DynamicBitset>& literals_in_class) {
    auto is_covered_by = [&](Vertex v, const DynamicBitset& bset) {
        return bset[v.first] && bset[v.second];
    };
    auto is_uncovered = [&](Vertex v) {
        return !std::any_of(
            literals_in_class.begin(), literals_in_class.end(),
            [&](const DynamicBitset& bset) { return is_covered_by(v, bset); });
    };
    std::vector<Vertex> remaining;
    std::copy_if(universe.begin(), universe.end(),
                 std::back_inserter(remaining), is_uncovered);
    return remaining;
}

static std::vector<std::vector<Index>>
make_classes_with_literal(const std::vector<DynamicBitset>& literals_in_class) {
    std::vector<std::vector<Index>> result;
    result.resize(literals_in_class[0].size());
    std::size_t cindex = 0;
    for (const DynamicBitset& bset : literals_in_class) {
        for (Lit l : bset.ones()) {
            result[l].push_back(cindex);
        }
        ++cindex;
    }
    return result;
}

TEST_CASE("[CUDA Iteration] Test with simple instance and partial solutions") {
    std::vector<ExternalClause> external_clauses = {{-1, 2}, {3, 4}, {2, -3}};
    Var num_vars = 4;
    Var n_concrete = 4;
    Lit nclit = 8;
    PairInfeasibilityMap inf_map(n_concrete);
    inf_map.literal_pair_infeasible(lit::internalize(1), lit::internalize(-2));
    inf_map.literal_pair_infeasible(lit::internalize(-3), lit::internalize(-4));
    inf_map.literal_pair_infeasible(lit::internalize(-2), lit::internalize(3));
    std::vector<DynamicBitset> literals_in_class;
    std::vector<std::vector<Index>> classes_with_literal;
    std::vector<Vertex> uncovered_vout;
    auto collect_uncovered = [&](Lit lmin, Lit lmax) {
        uncovered_vout.emplace_back(lmin, lmax);
    };
    std::vector<Vertex> universe;
    std::vector<Vertex> excluded{{lit::internalize(1), lit::internalize(-2)},
                                 {lit::internalize(-3), lit::internalize(-4)},
                                 {lit::internalize(-2), lit::internalize(3)}};
    for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
        for (Lit lmax = lmin + 1; lmax < nclit; ++lmax) {
            if (lit::negate(lmin) == lmax)
                continue;
            if (std::find(excluded.begin(), excluded.end(),
                          std::make_pair(lmin, lmax)) != excluded.end())
                continue;
            universe.emplace_back(lmin, lmax);
        }
    }

    SUBCASE("EMPTY PARTIAL") {
        cuda_iterate_all_uncovered(literals_in_class, classes_with_literal,
                                   &inf_map, collect_uncovered);
        CHECK(uncovered_vout == universe);
    }

    SUBCASE("ONE PARTIAL SAMPLE") {
        CHECK(uncovered_vout.empty());
        DynamicBitset c1(nclit, false);
        c1[lit::internalize(-1)] = true;
        c1[lit::internalize(2)] = true;
        literals_in_class.emplace_back(std::move(c1));
        classes_with_literal = make_classes_with_literal(literals_in_class);
        auto expected = filter_universe_if_covered(universe, literals_in_class);
        cuda_iterate_all_uncovered(literals_in_class, classes_with_literal,
                                   &inf_map, collect_uncovered);
        CHECK(uncovered_vout == expected);
    }
}
