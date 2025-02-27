#ifndef SAMMY_LITERALS_H_INCLUDED_
#define SAMMY_LITERALS_H_INCLUDED_

#include "time.h"
#include <algorithm>
#include <boost/container/small_vector.hpp>
#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_flat_set.hpp>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace sammy {

/**
 * A literal in the external sense:
 *  - Positive literals are positive numbers,
 *  - Negative literals are negative numbers,
 *  - 0 is not a valid literal.
 */
using ExternalLit = std::int32_t;

/**
 * An external clause, represented as vector of literals.
 */
using ExternalClause = std::vector<ExternalLit>;

/**
 * A literal in the internal sense:
 *  - Even numbers are positive literals,
 *  - Odd numbers are negative literals.
 *  - For variables x, we thus have 2 * x and 2 * x + 1 as literals.
 */
using Lit = std::uint32_t;

/**
 * An internal clause represented as vector of internal literals.
 */
using CVec = std::vector<Lit>;

/**
 * An internal clause represented as small vector of internal literals.
 */
using SCVec = boost::container::small_vector<Lit, 4>;

/**
 * A variable in the internal sense (technically equal type to Lit).
 * Only used to make types clearer in the interfaces.
 */
using Var = Lit;

/**
 * A clause reference to a clause in the database.
 */
using CRef = Lit;

/**
 * A vertex in the universe of pairwise interactions.
 */
using Vertex = std::pair<Lit, Lit>;

/**
 * A vertex in the universe of pairwise interactions
 * in external representation.
 */
using ExternalVertex = std::pair<ExternalLit, ExternalLit>;

/**
 * A value that indicates an invalid variable/literal/clause.
 */
static constexpr Lit NIL = std::numeric_limits<Lit>::max();

namespace lit {

/**
 * @brief External to internal literal conversion.
 *
 * @param l
 * @return Lit
 */
static inline constexpr Lit internalize(ExternalLit l) noexcept {
    return Lit(2) * static_cast<Lit>((std::abs)(l)-1) + Lit(l < 0);
}

static inline constexpr Vertex internalize(ExternalVertex v) noexcept {
    return {internalize(v.first), internalize(v.second)};
}

static inline std::vector<Vertex>
internalize(const std::vector<ExternalVertex>& vertices) {
    std::vector<Vertex> result;
    result.reserve(vertices.size());
    std::transform(vertices.begin(), vertices.end(), std::back_inserter(result),
                   [](ExternalVertex v) { return internalize(v); });
    return result;
}

/**
 * @brief Internal to external literal conversion.
 *
 * @param l
 * @return constexpr ExternalLit
 */
static inline constexpr ExternalLit externalize(Lit l) noexcept {
    ExternalLit result = static_cast<ExternalLit>(l >> 1) + 1;
    if (l & 1)
        result = -result;
    return result;
}

static inline std::vector<std::pair<ExternalLit, ExternalLit>>
externalize(const std::vector<Vertex>& vertices) {
    std::vector<std::pair<ExternalLit, ExternalLit>> result;
    result.reserve(vertices.size());
    std::transform(vertices.begin(), vertices.end(), std::back_inserter(result),
                   [](const auto& v) {
                       return std::pair<ExternalLit, ExternalLit>{
                           externalize(v.first), externalize(v.second)};
                   });
    return result;
}

/**
 * @brief Internal to external literal conversion.
 */
template <
    typename ClauseContainer,
    std::enable_if_t<std::is_integral_v<typename ClauseContainer::value_type>,
                     int> = 0>
static inline ExternalClause externalize(const ClauseContainer& clause) {
    ExternalClause result;
    std::transform(clause.begin(), clause.end(), std::back_inserter(result),
                   [](Lit l) -> ExternalLit { return externalize(l); });
    return result;
}

/**
 * @brief Internal literal negation.
 *
 * @param l
 * @return constexpr Lit
 */
static inline constexpr Lit negate(Lit l) noexcept { return l ^ Lit(1); }

/**
 * @brief Extract the variable from a literal.
 *
 * @param l
 * @return constexpr Lit
 */
static inline constexpr Var var(Lit l) noexcept { return l >> 1; }

/**
 * @brief Turn a variable into its positive literal.
 */
static inline constexpr Lit positive_lit(Var v) noexcept { return v << 1; }

/**
 * @brief Turn a variable into its negative literal.
 */
static inline constexpr Lit negative_lit(Var v) noexcept {
    return (v << 1) + 1;
}

/**
 * @brief Check for negative literal.
 *
 * @param l
 * @return true if the literal is negative.
 * @return false if the literal is positive.
 */
static inline constexpr bool negative(Lit l) noexcept { return l & Lit(1); }

template <typename BitsetType>
static inline bool is_true_in(Lit l, const BitsetType& assignment) noexcept {
    bool value(assignment[var(l)]);
    return negative(l) ? !value : value;
}

template <typename BitsetType>
static inline bool is_false_in(Lit l, const BitsetType& assignment) noexcept {
    return !is_true_in(l, assignment);
}

} // namespace lit

namespace simplify {

static constexpr Lit fixed_positive() noexcept { return NIL - 1; }

static constexpr Lit fixed_negative() noexcept { return NIL; }

static constexpr Lit eliminated() noexcept { return NIL - 2; }

} // namespace simplify

struct PairHash {
    std::size_t operator()(std::pair<int, int> v) const noexcept {
        std::size_t x1 = static_cast<unsigned>(v.first);
        std::size_t x2 = static_cast<unsigned>(v.second);
        x1 = ((x1 << ((CHAR_BIT / 2) * sizeof(std::size_t))) |
              (x1 >> ((CHAR_BIT / 2) * sizeof(std::size_t)))) *
             7;
        x1 += 17 * x2;
        return x1;
    }

    std::size_t
    operator()(std::pair<std::uint32_t, std::uint32_t> v) const noexcept {
        std::size_t x1 = v.first;
        std::size_t x2 = v.second;
        x1 = ((x1 << ((CHAR_BIT / 2) * sizeof(std::size_t))) |
              (x1 >> ((CHAR_BIT / 2) * sizeof(std::size_t)))) *
             7;
        x1 += 17 * x2;
        return x1;
    }

    std::size_t
    operator()(std::pair<std::size_t, std::size_t> v) const noexcept {
        std::size_t x = 7 * v.first;
        x = (x << 25) | (x >> (sizeof(x) * CHAR_BIT - 25));
        std::size_t y = 17 * v.second;
        return x + y;
    }

    using value = std::size_t;
};

template <typename Element> using HashSet = boost::unordered_flat_set<Element>;

template <typename Key, typename Value>
using HashMap = boost::unordered_flat_map<Key, Value>;

using EdgeSet = boost::unordered_flat_set<std::pair<Lit, Lit>>;

template <typename Value>
using VertexMapTo = boost::unordered_flat_map<std::pair<Lit, Lit>, Value>;

template <typename Key> using PairHashSet = boost::unordered_flat_set<Key>;

template <typename Key, typename Value>
using PairMapTo = boost::unordered_flat_map<Key, Value>;

} // namespace sammy

#endif
