#ifndef CLAUSE_REDUCTION_H_INCLUDED_
#define CLAUSE_REDUCTION_H_INCLUDED_

#include "literals.h"
#include <algorithm>
#include <optional>
#include <unordered_set>
#include <utility>
#include <vector>

namespace sammy {

/**
 * @brief Reduce an external clause by removing duplicates.
 * Detects tautological clauses (which contain x and -x); returns
 * an empty vector for them. Also sorts the literals.
 *
 * @param x
 * @return std::vector<int>
 */
inline ExternalClause reduce_external_clause(const ExternalClause& x) {
    std::vector<int> result;
    result.reserve(x.size());
    if (x.size() > 256) {
        // use an external set for large clauses.
        std::unordered_set<int> entries;
        for (int l : x) {
            if (entries.count(-l)) {
                // x is a tautology; drop it.
                result.clear();
                return result;
            }
            if (entries.insert(l).second) {
                // filter duplicates.
                result.push_back(l);
            }
        }
    } else {
        // use the result vector for small clauses.
        for (int l : x) {
            bool found = false;
            for (int rl : result) {
                if (rl == -l) {
                    // x is a tautology.
                    result.clear();
                    return result;
                }
                if (rl == l) {
                    found = true;
                }
            }
            if (!found) {
                // filter duplicates.
                result.push_back(l);
            }
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

namespace impl {

static inline void
reduced_remove_identical_pairs(std::vector<ExternalClause>& reduced) {
    auto comp_size = [](const auto& v1, const auto& v2) {
        return v1.size() < v2.size();
    };
    std::sort(reduced.begin(), reduced.end(), comp_size);
    ExternalClause l2(2, 0);
    auto begin_binary =
        std::lower_bound(reduced.begin(), reduced.end(), l2, comp_size);
    auto end_binary =
        std::upper_bound(reduced.begin(), reduced.end(), l2, comp_size);
    std::unordered_set<std::pair<int, int>, PairHash> binaries;
    auto out_binary = begin_binary;
    for (; begin_binary != end_binary; ++begin_binary) {
        std::pair<int, int> p{(*begin_binary)[0], (*begin_binary)[1]};
        if (binaries.insert(p).second) {
            if (out_binary != begin_binary) {
                *out_binary = std::move(*begin_binary);
            }
            ++out_binary;
        }
    }
    if (out_binary != begin_binary) {
        reduced.erase(std::move(end_binary, reduced.end(), out_binary),
                      reduced.end());
    }
}

} // namespace impl

/**
 * @brief Reduce external clauses, removing duplicates, tautologies and subsumed
 * clauses.
 *
 * @param clauses
 * @return std::vector<std::vector<int>>
 */
inline std::vector<ExternalClause>
reduce_external_clauses(const std::vector<ExternalClause>& clauses) {
    std::vector<std::vector<int>> reduced;
    reduced.reserve(clauses.size());
    for (const auto& c : clauses) {
        auto res = reduce_external_clause(c);
        if (!res.empty()) {
            reduced.emplace_back(std::move(res));
        }
    }
    impl::reduced_remove_identical_pairs(reduced);
    return reduced;
}

} // namespace hs

#endif
