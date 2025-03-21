#ifndef SAMMY_IMPLIED_VERTICES_H_INCLUDED_
#define SAMMY_IMPLIED_VERTICES_H_INCLUDED_

#include "algorithm_ex.h"
#include "dynamic_bitset.h"
#include "literals.h"
#include "pair_infeasibility_map.h"
#include "range.h"
#include "shared_db_propagator.h"
#include "vertex_operations.h"

namespace sammy {

/**
 * On demand elimination of implied vertices.
 * Returns a list of vertices such that covering all of them guarantees a full
 * cover, by removing all vertices that are implied by another vertex.
 */
inline std::vector<Vertex>
eliminate_implied_vertices(const std::vector<Vertex>& vertices,
                           SharedDBPropagator& propagator) {
    if (propagator.get_current_level() > 0)
        throw std::logic_error(
            "eliminate_implied_vertices: decision level > 0");
    DynamicBitset implied(vertices.size(), false);
    std::vector<std::vector<std::pair<Lit, std::size_t>>> partners_of(
        2 * propagator.db().num_vars());
    for (std::size_t i = 0, n = vertices.size(); i < n; ++i) {
        Vertex v = vertices[i];
        if (propagator.is_true(v.first) && propagator.is_true(v.second)) {
            implied[i].set();
            continue;
        }
        partners_of[v.first].emplace_back(v.second, i);
        partners_of[v.second].emplace_back(v.first, i);
    }
    for (std::size_t i = 0, n = vertices.size(); i < n; ++i) {
        if (implied[i])
            continue;
        Vertex v = vertices[i];
        reset_and_push_noresolve(propagator, v);
        for (Lit l : propagator.get_trail()) {
            auto& ps = partners_of[l];
            ps.erase(std::remove_if(ps.begin(), ps.end(),
                                    [&](std::pair<Lit, std::size_t> entry) {
                                        if (implied[entry.second])
                                            return true;
                                        if (entry.second == i)
                                            return false;
                                        if (propagator.is_true(entry.first)) {
                                            implied[entry.second].set();
                                            return true;
                                        }
                                        return false;
                                    }),
                     ps.end());
        }
    }
    std::vector<Vertex> result;
    for (std::size_t i = 0, n = vertices.size(); i < n; ++i) {
        if (!implied[i])
            result.push_back(vertices[i]);
    }
    propagator.reset_or_throw();
    return result;
}

inline std::vector<Vertex>
eliminate_implied_vertices(const std::vector<Vertex>& vertices,
                           const SharedDBPropagator& propagator) {
    SharedDBPropagator propagator_(propagator);
    return eliminate_implied_vertices(vertices, propagator_);
}

/**
 * Lookup class for implied vertices.
 * Serves as storage for the entire
 * (reduced) universe of vertices.
 * Offers four main operations:
 *  - Eliminate all implied vertices,
 *    i.e., vertices that are implicitly covered
 *    if another vertex is covered.
 *  - Check if a single vertex is implied, and if so,
 *    by which non-implied vertex.
 *  - Remove all implied vertices from a given list.
 *    Usually, when this is called, the list
 *    to be reduced will contain the implying vertex as well.
 *    Used primarily to reduce the set of uncovered vertices
 *    in a predictable and consistent manner.
 *  - Replace all implied vertices by their implier.
 *    Usually, when this is called, the list
 *    does not contain both a vertex and its implier.
 *    Used primarily to replace vertices in a mutually exclusive set/clique
 *    by canonical implying vertices.
 * Offers a 'partial'/'limited' elimination version that may
 * let some implied vertices slip through, trading off vs. elimination
 * runtime.
 */
class ImpliedVertexCache {
  public:
    explicit ImpliedVertexCache(const PairInfeasibilityMap* inf_map,
                                std::size_t universe_size) noexcept
        : m_pair_inf_map(inf_map), m_universe_size(universe_size) {}

    /**
     * Compute and cache a reduced version of the universe.
     */
    inline void reduce_universe(ClauseDB& clause_db);

    /**
     * Time-limited version of universe reduction.
     */
    inline void limited_reduce_universe(ClauseDB& clause_db, double time_limit);

    /**
     * Check whether we have computed and stored a reduced universe.
     */
    bool have_reduced_universe() const noexcept {
        return !m_reduced_universe.empty();
    }

    /**
     * Return the size of the reduced universe.
     */
    std::size_t reduced_universe_size() const noexcept {
        return m_reduced_universe.size();
    }

    /**
     * Original universe size.
     */
    std::size_t original_universe_size() const noexcept {
        return m_universe_size;
    }

    /**
     * Get a reference to the reduced universe.
     * Returns an empty vector if we have no reduced universe.
     */
    const std::vector<Vertex>& get_reduced_universe() const noexcept {
        return m_reduced_universe;
    }

    /**
     * Check if the given vertex is implied by another vertex.
     */
    bool is_implied(Vertex v) const noexcept { return m_implied_by.count(v); }

    /**
     * Get the vertex implying the given vertex (or an exception).
     */
    Vertex implying_vertex(Vertex v) const noexcept {
        return m_implied_by.at(v);
    }

    /**
     * Get the vertex implying the given vertex,
     * or the given vertex if the vertex is not implied.
     */
    Vertex implying_or_self(Vertex v) const noexcept {
        auto it = m_implied_by.find(v);
        if (it == m_implied_by.end())
            return v;
        return it->second;
    }

    /**
     * Get a copy of the reduced universe, or the entire universe,
     * if we have no reduced universe.
     */
    std::vector<Vertex> get_universe() const {
        if (have_reduced_universe())
            return m_reduced_universe;
        return m_pair_inf_map->collect_vertices(m_universe_size);
    }

    /**
     * Remove the implied vertices from the given list
     * (throws if universe reduction has not been performed).
     */
    void remove_implied(std::vector<Vertex>& vertices) const {
        if (!have_reduced_universe())
            throw std::logic_error("Universe reduction not performed!");
        auto new_end = std::remove_if(vertices.begin(), vertices.end(),
                                      [&](Vertex v) { return is_implied(v); });
        vertices.erase(new_end, vertices.end());
    }

    /**
     * Remove the implied vertices from the given list (either through the
     * ad-hoc method, or using a cached reduced universe).
     */
    void remove_implied(std::vector<Vertex>& vertices,
                        SharedDBPropagator& propagator) const {
        if (!have_reduced_universe()) {
            vertices = eliminate_implied_vertices(vertices, propagator);
        } else {
            remove_implied(vertices);
        }
    }

    /**
     * Remove the implied vertices from the given list (const propagator
     * version).
     */
    void remove_implied(std::vector<Vertex>& vertices,
                        const SharedDBPropagator& propagator) const {
        if (!have_reduced_universe()) {
            SharedDBPropagator copy(propagator);
            remove_implied(vertices, copy);
        } else {
            remove_implied(vertices);
        }
    }

    /**
     * Replace all implied vertices in the given list
     * by their implier; does nothing if we do not actually
     * have a reduced universe.
     */
    void replace_implied(std::vector<Vertex>& vertices) const {
        if (!have_reduced_universe())
            return;
        std::transform(vertices.begin(), vertices.end(), vertices.begin(),
                       [&](Vertex v) { return implying_or_self(v); });
    }

  private:
    /**
     * Member struct that implements the actual elimination algorithm.
     */
    struct EliminationAlgorithm {
        explicit EliminationAlgorithm(ClauseDB& clauses,
                                      ImpliedVertexCache* that)
            : that(that), propagator(&clauses),
              universe(that->m_pair_inf_map->collect_vertices(
                  that->m_universe_size)),
              implier_of(universe.size(),
                         std::numeric_limits<std::size_t>::max()) {}

        bool is_implied(std::size_t vi) const noexcept {
            return implier_of[vi] != std::numeric_limits<std::size_t>::max();
        }

        inline void
        compute_impliers_handle_trail_literal(Lit trail_literal,
                                              std::size_t pushed_index) __attribute__((noinline));
        inline void mark_and_merge_two_sorted_lists(Lit implier, Lit implied) __attribute__((noinline));
        inline void compute_literal_partners_of() __attribute__((noinline));
        inline void compute_impliers() __attribute__((noinline));
        inline void limited_compute_impliers(double time_limit) __attribute__((noinline));
        inline void compute_impliers_single_literal() __attribute__((noinline));
        inline void compress_path(std::size_t v) __attribute__((noinline));

        /**
         * Compress all paths of length > 1 by
         * replacing each intermediate implier
         * by the corresponding ultimate implier,
         * i.e., non-implied vertex.
         */
        void compress_paths() {
            for (std::size_t v : range(universe.size()))
                compress_path(v);
        }

        inline void export_to_cache();

        /**
         * A reference to the ImpliedVertexCache we are doing computation for.
         */
        ImpliedVertexCache* that;

        /**
         * A propagator that we use to detect implications between vertices.
         */
        SharedDBPropagator propagator;

        /**
         * The full, unreduced universe.
         */
        std::vector<Vertex> universe;

        /**
         * For each literal, the list of vertices that contain that literal
         * and the other literal that makes up the vertex.
         */
        std::vector<std::vector<std::pair<Lit, std::size_t>>> partners_of;

        /**
         * For each vertex (indices into universe), the vertex that implies it,
         * or std::numeric_limits<std::size_t>::max() if there is no implier.
         */
        std::vector<std::size_t> implier_of;

        /**
         * Cache of pointers we need to update during path compression.
         */
        std::vector<std::size_t*> pcompress_cache;
    };

    std::vector<Vertex> m_reduced_universe;
    VertexMapTo<Vertex> m_implied_by;
    const PairInfeasibilityMap* m_pair_inf_map;
    std::size_t m_universe_size;
};

void ImpliedVertexCache::EliminationAlgorithm::compute_literal_partners_of() {
    partners_of.resize(2 * propagator.db().num_vars());
    for (std::size_t i = 0, n = universe.size(); i < n; ++i) {
        Vertex v = universe[i];
        if (propagator.is_true(v.first) && propagator.is_true(v.second)) {
            implier_of[i] = i;
            continue;
        }
        partners_of[v.first].emplace_back(v.second, i);
        partners_of[v.second].emplace_back(v.first, i);
    }
    auto sort_partners_list = [](auto& r) {
        const auto compare_1st = [](const auto& e1, const auto& e2) {
            return e1.first < e2.first;
        };
        std::sort(r.begin(), r.end(), compare_1st);
    };
    std::for_each(partners_of.begin(), partners_of.end(), sort_partners_list);
}

void ImpliedVertexCache::EliminationAlgorithm::mark_and_merge_two_sorted_lists(
    Lit implier, Lit implied) {
    auto& left_list = partners_of[implier];
    auto& right_list = partners_of[implied];
    if (left_list.empty() || right_list.empty())
        return;

    auto left_in = left_list.begin(), left_out = left_list.begin(),
         left_end = left_list.end();
    auto right_in = right_list.begin(), right_out = right_list.begin(),
         right_end = right_list.end();

    auto left_skip_implied = [&]() {
        while (is_implied(left_in->second)) {
            if (++left_in == left_end)
                return false;
        }
        return true;
    };

    auto right_skip_implied = [&]() {
        while (is_implied(right_in->second)) {
            if (++right_in == right_end)
                return false;
        }
        return true;
    };

    auto advance_left = [&]() {
        *left_out++ = *left_in;
        if (++left_in == left_end) {
            return false;
        }
        return left_skip_implied();
    };

    auto advance_right = [&]() {
        *right_out++ = *right_in;
        if (++right_in == right_end) {
            return false;
        }
        return right_skip_implied();
    };

    auto advance_both = [&]() {
        *left_out++ = *left_in;
        ++left_in, ++right_in;
        if (left_in != left_end) {
            left_skip_implied();
        }
        if (right_in != right_end) {
            right_skip_implied();
        }
        return left_in != left_end && right_in != right_end;
    };

    auto handle_left_remainder = [&]() {
        if (left_in == left_end)
            return;
        while (advance_left())
            ;
    };

    auto handle_right_remainder = [&]() {
        if (right_in == right_end)
            return;
        while (advance_right())
            ;
    };

    if (!left_skip_implied()) {
        left_list.clear();
        if (right_skip_implied()) {
            handle_right_remainder();
        }
        right_list.erase(right_out, right_end);
        return;
    }

    if (!right_skip_implied()) {
        right_list.clear();
        handle_left_remainder();
        left_list.erase(left_out, left_end);
        return;
    }

    for (;;) {
        Lit left = left_in->first, right = right_in->first;
        if (left < right) {
            if (!advance_left())
                break;
        } else if (right < left) {
            if (!advance_right())
                break;
        } else {
            assert(left_in->second != right_in->second);
            implier_of[right_in->second] = left_in->second;
            if (!advance_both())
                break;
        }
    }
    handle_right_remainder();
    handle_left_remainder();
    left_list.erase(left_out, left_end);
    right_list.erase(right_out, right_end);
}

void ImpliedVertexCache::EliminationAlgorithm::
    compute_impliers_single_literal() {
    for (Lit l = 0, nall = 2 * propagator.db().num_vars(); l < nall; ++l) {
        // find any non-implied interaction involving l;
        // avoid doing propagation if there are no un-implied
        // interactions involving l
        auto& p_of_l = partners_of[l];
        while (!p_of_l.empty() && is_implied(p_of_l.back().second)) {
            p_of_l.pop_back();
        }
        if (p_of_l.empty())
            continue;
        if (!propagator.push_level(l)) {
            throw std::logic_error("Infeasible interaction in universe!");
        }
        for (Lit l2 : propagator.get_trail()) {
            if (l2 == l)
                continue;
            // if l implies l2, (l2, x) is implied by (l, x) for any x.
            // so we walk the sorted partner lists of l and l2 and
            // mark (l2,x) for x for which (l,x) is present,
            // stripping both lists of implied elements as we go
            mark_and_merge_two_sorted_lists(l, l2);
        }
        propagator.pop_level();
    }
}

void ImpliedVertexCache::EliminationAlgorithm::
    compute_impliers_handle_trail_literal(Lit trail_literal,
                                          std::size_t pushed_index) {
    auto& ps = partners_of[trail_literal];
    ps.erase(std::remove_if(ps.begin(), ps.end(),
                            [&](std::pair<Lit, std::size_t> entry) {
                                // simply drop already-implied vertices
                                if (is_implied(entry.second))
                                    return true;
                                // avoid detecting vertices as implying
                                // themselves
                                if (entry.second == pushed_index)
                                    return false;
                                // if the other literal is true, pushed_index
                                // implies entry.first
                                if (propagator.is_true(entry.first)) {
                                    implier_of[entry.second] = pushed_index;
                                    return true;
                                }
                                // not known to be implied
                                return false;
                            }),
             ps.end());
}

void ImpliedVertexCache::EliminationAlgorithm::compute_impliers() {
    for (std::size_t i : range(universe.size())) {
        if (is_implied(i))
            continue;
        Vertex v = universe[i];
        // TODO: it should suffice to check the level 2 part of the trail since
        // we cover implied-by-single-literal earlier;
        // it could also help to only push one literal and pop one literal
        // unless the first literal of our vertex has changed.
        // other datastructure/algorithm ideas?
        // could create a graph-like structure for each level 1 literal,
        // using all second literals of not-yet-implied vertices, but this
        // sounds rather complex
        reset_and_push_noresolve(propagator, v);
        for (Lit l : propagator.get_trail()) {
            compute_impliers_handle_trail_literal(l, i);
        }
    }
}

void ImpliedVertexCache::EliminationAlgorithm::limited_compute_impliers(
    double time_limit) {
    auto begin_time = std::chrono::steady_clock::now();
    std::vector<std::size_t> vertex_indices;
    for (std::size_t i = 0, us = universe.size(); i < us; ++i) {
        if (is_implied(i))
            continue;
        vertex_indices.push_back(i);
    }
    std::shuffle(vertex_indices.begin(), vertex_indices.end(), sammy::rng());
    std::size_t check_count = 0;
    for (std::size_t vertex_index : vertex_indices) {
        if (is_implied(vertex_index))
            continue;
        Vertex vertex = universe[vertex_index];
        reset_and_push_noresolve(propagator, vertex);
        for (Lit l : propagator.get_trail()) {
            compute_impliers_handle_trail_literal(l, vertex_index);
        }
        if (++check_count == 16384) {
            check_count = 0;
            if (seconds_between(begin_time, std::chrono::steady_clock::now()) >=
                time_limit)
            {
                break;
            }
        }
    }
}

void ImpliedVertexCache::EliminationAlgorithm::compress_path(std::size_t vi) {
    // check non-impliedness/self-impliedness
    std::size_t* current = &implier_of[vi];
    std::size_t cval = *current;
    if (cval == std::numeric_limits<std::size_t>::max() || cval == vi)
        return;
    // get ready to compress the path
    pcompress_cache.clear();
    // find entries we need to update
    for (;;) {
        std::size_t* next = &implier_of[cval];
        std::size_t nval = *next;
        if (nval == std::numeric_limits<std::size_t>::max())
            break;
        pcompress_cache.push_back(current);
        cval = nval;
        current = next;
    }
    // update all entries to ultimate implier
    for (std::size_t* ptr : pcompress_cache) {
        *ptr = cval;
    }
}

void ImpliedVertexCache::EliminationAlgorithm::export_to_cache() {
    auto& implied_by = that->m_implied_by;
    auto& unimplied = that->m_reduced_universe;
    for (std::size_t vi : range(universe.size())) {
        Vertex v = universe[vi];
        if (implier_of[vi] != std::numeric_limits<std::size_t>::max()) {
            Vertex implier = universe[implier_of[vi]];
            implied_by.insert_or_assign(v, implier);
        } else {
            unimplied.push_back(v);
        }
    }
}

void ImpliedVertexCache::reduce_universe(ClauseDB& clause_db) {
    EliminationAlgorithm algorithm(clause_db, this);
    algorithm.compute_literal_partners_of();
    algorithm.compute_impliers_single_literal();
    algorithm.compute_impliers();
    algorithm.compress_paths();
    algorithm.export_to_cache();
}

void ImpliedVertexCache::limited_reduce_universe(ClauseDB& clause_db,
                                                 double time_limit) {
    if (time_limit >= 0.0 && !std::isfinite(time_limit)) {
        reduce_universe(clause_db);
        return;
    }
    auto begin_time = std::chrono::steady_clock::now();
    EliminationAlgorithm algorithm(clause_db, this);
    algorithm.compute_literal_partners_of();
    algorithm.compute_impliers_single_literal();
    double trem = time_limit - 
        seconds_between(begin_time, std::chrono::steady_clock::now());
    if (trem > 0.0) {
        algorithm.limited_compute_impliers(trem);
    }
    algorithm.compress_paths();
    algorithm.export_to_cache();
}

} // namespace sammy

#endif
