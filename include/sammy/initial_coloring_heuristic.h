#ifndef SAMMY_INITIAL_COLORING_HEURISTIC_H_INCLUDED_
#define SAMMY_INITIAL_COLORING_HEURISTIC_H_INCLUDED_

#include "algorithm_ex.h"
#include "class_completer.h"
#include "cuda_iteration.h"
#include "error.h"
#include "greedysat.h"
#include "literals.h"
#include "pair_infeasibility_map.h"
#include "parallel_bit_filter.h"
#include "partial_solution.h"
#include <boost/iterator/transform_iterator.hpp>
#include <unordered_set>

namespace sammy {

#ifdef SAMMY_CUDA_SUPPORTED
namespace detail {

class CUDABitFilter {
  public:
    CUDABitFilter(const std::vector<DynamicBitset>& literals_in_class,
                  const std::vector<std::vector<Index>>& classes_with_literal,
                  const PairInfeasibilityMap* inf_map)
        : m_inf_map(inf_map), m_u32_per_bitset((*inf_map)[0].blocks().size() *
                                               sizeof(DynamicBitset::Block) /
                                               sizeof(std::uint32_t)),
          m_device_bit_data((p_fill_prepare_buffer(literals_in_class),
                             m_host_prepare_buffer)),
          m_device_classes_with_literal(
              (p_fill_offset_buffer(classes_with_literal),
               m_host_prepare_buffer)),
          m_device_classes_with_literal_offsets(m_host_prepare_offsets),
          m_device_output_buffer(GOAL_ROWS_PER_CALL() * m_u32_per_bitset) {}

    template <typename Callable /*(Lit lmin, Lit lmax)*/>
    void iterate_all_uncovered(Callable&& callable) {
        Lit n_concrete = m_inf_map->get_n_concrete();
        Lit nclit = 2 * n_concrete;
        Lit lmin = 0;
        std::vector<std::uint32_t> host_output_buffer;
        DynamicBitset row_buffer(nclit, true);
        while (lmin < nclit - 2) {
            Lit num_rows = nclit - lmin - 2;
            if (num_rows > GOAL_ROWS_PER_CALL()) {
                num_rows = GOAL_ROWS_PER_CALL();
            }
            cuda_call_bit_filter_kernel(
                m_device_bit_data.get(), m_u32_per_bitset,
                m_device_classes_with_literal.get(), nclit,
                m_device_classes_with_literal_offsets.get(),
                m_device_output_buffer.get(), lmin, num_rows);
            m_device_output_buffer.copy_to_host(host_output_buffer);
            for (Lit l = 0; l < num_rows; ++l) {
                Lit lm = lmin + l;
                row_buffer.set();
                row_buffer ^= (*m_inf_map)[lm];
                row_buffer.binary_subtract(
                    &host_output_buffer[l * m_u32_per_bitset]);
                for (Lit lmax : row_buffer.ones_from(lm + 1)) {
                    callable(lm, lmax);
                }
            }
            lmin += num_rows;
        }
    }

    std::size_t count_uncovered() {
        Lit n_concrete = m_inf_map->get_n_concrete();
        Lit nclit = 2 * n_concrete;
        Lit lmin = 0;
        std::vector<std::uint32_t> host_output_buffer;
        DynamicBitset row_buffer(nclit, true);
        std::size_t result = 0;
        while (lmin < nclit - 2) {
            Lit num_rows = nclit - lmin - 2;
            if (num_rows > GOAL_ROWS_PER_CALL()) {
                num_rows = GOAL_ROWS_PER_CALL();
            }
            cuda_call_bit_filter_kernel(
                m_device_bit_data.get(), m_u32_per_bitset,
                m_device_classes_with_literal.get(), nclit,
                m_device_classes_with_literal_offsets.get(),
                m_device_output_buffer.get(), lmin, num_rows);
            m_device_output_buffer.copy_to_host(host_output_buffer);
            for (Lit l = 0; l < num_rows; ++l) {
                Lit lm = lmin + l;
                row_buffer.set();
                row_buffer ^= (*m_inf_map)[lm];
                row_buffer.binary_subtract(
                    &host_output_buffer[l * m_u32_per_bitset]);
                result += row_buffer.count_from(lm + 1);
            }
            lmin += num_rows;
        }
        return result;
    }

  private:
    void
    p_fill_prepare_buffer(const std::vector<DynamicBitset>& literals_in_class) {
        if (literals_in_class.empty()) {
            throw std::invalid_argument("literals_in_class must not be empty");
        }
        m_host_prepare_buffer.reserve(literals_in_class.size() *
                                      m_u32_per_bitset);
        for (const DynamicBitset& bs : literals_in_class) {
            for (DynamicBitset::Block b : bs.blocks()) {
                to_prepare_buffer(b, m_host_prepare_buffer);
            }
        }
    }

    void p_fill_offset_buffer(
        const std::vector<std::vector<Index>>& classes_with_literal) {
        m_host_prepare_buffer.clear();
        m_host_prepare_offsets.reserve(classes_with_literal.size() + 1);
        m_host_prepare_offsets.push_back(0);
        for (const std::vector<Index>& cls : classes_with_literal) {
            std::copy(cls.begin(), cls.end(),
                      std::back_inserter(m_host_prepare_buffer));
            m_host_prepare_offsets.push_back(m_host_prepare_buffer.size());
        }
    }

    // pair infeasibility map
    const PairInfeasibilityMap* m_inf_map;
    // words per bitset in m_host_prepare_buffer and on the device
    std::size_t m_u32_per_bitset;
    // buffer on the host to prepare the bitsets for copy,
    // or the classes with literal list
    std::vector<std::uint32_t> m_host_prepare_buffer;
    // host buffer for the offsets in classes with literal
    std::vector<std::uint32_t> m_host_prepare_offsets;

    // pointer to array of concatenated bitsets
    CUDADevicePointer<const std::uint32_t> m_device_bit_data;
    // pointer to array of classes with literal
    CUDADevicePointer<const std::uint32_t> m_device_classes_with_literal;
    // pointer to offsets within m_device_classes_with_literal to the start of
    // each literal
    CUDADevicePointer<const std::uint32_t>
        m_device_classes_with_literal_offsets;
    // pointer to output matrix rows
    CUDADevicePointer<std::uint32_t> m_device_output_buffer;
};

} // namespace detail

template <typename Callable>
static inline void cuda_iterate_all_uncovered(
    const std::vector<DynamicBitset>& literals_in_class,
    const std::vector<std::vector<Index>>& classes_with_literal,
    const PairInfeasibilityMap* inf_map, Callable&& callable) {
    if (literals_in_class.empty()) {
        Lit nclit = 2 * inf_map->get_n_concrete();
        DynamicBitset row_buffer(nclit, true);
        for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
            row_buffer.set();
            row_buffer ^= (*inf_map)[lmin];
            for (Lit lmax : row_buffer.ones_from(lmin + 1)) {
                std::forward<Callable>(callable)(lmin, lmax);
            }
        }
        return;
    }

    detail::CUDABitFilter bit_filter(literals_in_class, classes_with_literal,
                                     inf_map);
    bit_filter.iterate_all_uncovered(std::forward<Callable>(callable));
}

static inline std::size_t cuda_count_uncovered(
    const std::vector<DynamicBitset>& literals_in_class,
    const std::vector<std::vector<Index>>& classes_with_literal,
    const PairInfeasibilityMap* inf_map) {
    if (literals_in_class.empty()) {
        Lit nclit = 2 * inf_map->get_n_concrete();
        DynamicBitset row_buffer(nclit, true);
        std::size_t result = 0;
        for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
            row_buffer.set();
            row_buffer ^= (*inf_map)[lmin];
            result += row_buffer.count_from(lmin + 1);
        }
        return result;
    }

    detail::CUDABitFilter bit_filter(literals_in_class, classes_with_literal,
                                     inf_map);
    return bit_filter.count_uncovered();
}
#endif

/**
 * A class containing a set of propagators,
 * managing the incremental construction of
 * several configurations simultaneously.
 */
template <typename EventListener> class ColorClasses {
    /**
     * @brief A handler for the various events that
     *        can occur during our operations.
     */
    EventListener* event_handler;

    /**
     * @brief The database of all clauses,
     *        shared by all propagators.
     */
    ClauseDB* all_clauses;

    /**
     * @brief A propagator that is kept 'empty',
     *        i.e., with no (permanent) decisions;
     *        it gets all clauses and serves as a source to
     *        copy-init new color classes from
     *        (much faster than re-initialization).
     */
    SharedDBPropagator empty_class;

    /**
     * @brief Each color class is represented by the
     *        literals in the trail of a propagator.
     */
    std::vector<SharedDBPropagator> color_classes;

    /**
     * @brief For each literal, a list of class indices
     *        in which that literal is forced to true.
     */
    std::vector<std::vector<Index>> classes_with_literal;

    /**
     * @brief The list of vertices that spawned the color classes,
     *        i.e., that were first given that color.
     */
    std::vector<Vertex> class_spawners;

    /**
     * @brief An internal buffer of the phases
     *        of variables for p_try_completion.
     */
    std::vector<bool> m_phases;

    /**
     * @brief An internal buffer of the decisions
     *        for p_try_completion.
     */
    std::vector<Lit> m_decisions;

    /**
     * @brief Add a new color class.
     *
     * @param prop The propagator for the new class.
     * @param spawner The vertex (feature literal pair) that spawned the class.
     */
    void p_add_class(SharedDBPropagator&& prop, Vertex spawner) {
        class_spawners.emplace_back((std::min)(spawner.first, spawner.second),
                                    (std::max)(spawner.first, spawner.second));
        color_classes.emplace_back(std::move(prop));
        Index idx = num_classes() - 1;
        for (Lit l : color_classes.back().get_trail()) {
            classes_with_literal[l].push_back(idx);
        }
        event_handler->new_color_class(idx);
    }

    bool p_spawn_class_one_is_false(Vertex spawner) const noexcept {
        if (empty_class.is_false(spawner.first) ||
            empty_class.is_false(spawner.second))
        {
            if (empty_class.is_false(spawner.first)) {
                event_handler->feature_literal_infeasible(spawner.first);
            }
            if (empty_class.is_false(spawner.second)) {
                event_handler->feature_literal_infeasible(spawner.second);
            }
            return true;
        }
        return false;
    }

    void p_init_phases(SharedDBPropagator& propagator) {
        const Lit nv = all_clauses->num_vars();
        m_phases.assign(nv, false);
        for (Lit v = 0; v < nv; ++v) {
            if (propagator.is_true(lit::positive_lit(v))) {
                m_phases[v] = true;
            }
        }
    }

    /**
     * @brief Get a spawner from a propagator by finding the first concrete
     * literals on different levels.
     */
    std::optional<Vertex> p_spawner_from_levels(const SharedDBPropagator& prop,
                                                Lit n_concrete) {
        const std::uint32_t nlvl = prop.get_current_level();
        if (nlvl < 2)
            return std::nullopt;
        const Lit nclit = 2 * n_concrete;
        auto&& is_concrete = [&](Lit l) { return l < nclit; };
        for (std::uint32_t level = 1; level <= nlvl; ++level) {
            auto lbeg = prop.level_begin(level);
            auto lend = prop.level_end(level);
            auto lf = std::find_if(lbeg, lend, is_concrete);
            if (lf != lend) {
                for (++level; level <= nlvl; ++level) {
                    auto lbeg2 = prop.level_begin(level);
                    auto lend2 = prop.level_end(level);
                    auto lf2 = std::find_if(lbeg2, lend2, is_concrete);
                    if (lf2 != lend2) {
                        Lit l1 = *lf, l2 = *lf2;
                        return std::optional<Vertex>{
                            std::in_place, std::min(l1, l2), std::max(l1, l2)};
                    }
                }
                break;
            }
        }
        return std::nullopt;
    }

    /**
     * @brief Called if conflict resolution took us below
     *        the level introduced by the last assumption (decision
     *        that we wanted to keep) in p_try_completion.
     *
     * @param propagator
     * @return true If we could fix the situation, i.e., reintroduce a decision,
     *              or a decision turned into an implied literal.
     * @return false If the set of assumptions is infeasible.
     */
    bool p_try_completion_jumped_decision(SharedDBPropagator& propagator) {
        for (Lit l : m_decisions) {
            // if any of our assumptions became forced to false, return false.
            if (propagator.is_false(l))
                return false;
            // if any of our assumptions became open, try closing it.
            if (propagator.is_open(l)) {
                // if that yields a conflict, resolve and return false.
                if (!propagator.push_level(l)) {
                    if (!propagator.resolve_conflicts())
                        throw UNSATError();
                    return false;
                }
            }
        }
        // remove all decisions that are decisions no more
        m_decisions.erase(
            std::remove_if(m_decisions.begin(), m_decisions.end(),
                           [&](Lit d) { return !propagator.is_decision(d); }),
            m_decisions.end());
        return true;
    }

    bool p_try_completion_push(SharedDBPropagator& propagator, Lit v) {
        Lit d = m_phases[v] ? lit::positive_lit(v) : lit::negative_lit(v);
        if (!propagator.push_level(d)) {
            m_phases[v].flip();
            if (!propagator.resolve_conflicts()) {
                throw UNSATError();
            }
            if (std::size_t(propagator.get_current_level()) <
                m_decisions.size())
            {
                if (!p_try_completion_jumped_decision(propagator))
                    return false;
            }
            if (propagator.is_open(d)) {
                return p_try_completion_push(propagator, v);
            }
        }
        return true;
    }

    /**
     * @brief Try to complete the given partial assignment
     *        (actually solving SAT with a very simple variable selection
     * strategy).
     *
     * @param propagator
     * @return true If the solve was successful (the propagator is returned to a
     * state with the same assumptions as it had initially).
     * @return false If the solve was unsuccessful (the propagator is in some
     * unspecified state); the given assumptions are infeasible, and clauses
     * have been learned to reflect that.
     */
    bool p_try_completion(SharedDBPropagator& propagator) {
        p_init_phases(propagator);
        const Lit nv = all_clauses->num_vars();
        m_decisions = propagator.get_decisions();
        bool any_open = true;
        while (any_open) {
            any_open = false;
            for (Lit v = 0; v < nv; ++v) {
                Lit p = lit::positive_lit(v);
                if (propagator.is_open(p)) {
                    any_open = true;
                    if (!p_try_completion_push(propagator, v))
                        return false;
                }
            }
        }
        while (std::size_t(propagator.get_current_level()) > m_decisions.size())
        {
            propagator.pop_level();
        }
        return true;
    }

    struct ConflictResolutionHandler {
        void assignment_forced(Lit l) const {
            that->classes_with_literal[l].push_back(cls_index);
            that->event_handler->literal_added_to_class(cls_index, l);
        }

        void assignment_undone(Lit l) const {
            auto& list = that->classes_with_literal[l];
            list.erase(std::remove(list.begin(), list.end(), cls_index),
                       list.end());
            that->event_handler->literal_removed_from_class(cls_index, l);
        }

        ColorClasses* that;
        Index cls_index;
    };

  public:
    std::vector<SharedDBPropagator>
    remove_classes(const std::vector<std::size_t>& sorted_classes) {
        if (sorted_classes.empty())
            return {};
        std::vector<SharedDBPropagator> result;
        result.reserve(sorted_classes.size());
        for (std::size_t idx : sorted_classes) {
            result.emplace_back(std::move(color_classes[idx]));
        }
        class_spawners.erase(
            remove_indices(class_spawners.begin(), class_spawners.end(),
                           sorted_classes.begin(), sorted_classes.end()),
            class_spawners.end());
        const auto n = color_classes.size();
        std::vector<Index> old_to_new(n, std::numeric_limits<Index>::max());
        auto iter = sorted_classes.begin();
        Index new_out = 0;
        for (Index i = 0; i < n; ++i) {
            if (iter != sorted_classes.end() && *iter == i) {
                ++iter;
            } else {
                old_to_new[i] = new_out++;
            }
        }
        color_classes.erase(
            remove_indices(color_classes.begin(), color_classes.end(),
                           sorted_classes.begin(), sorted_classes.end()),
            color_classes.end());
        for (auto& cwl : classes_with_literal) {
            std::transform(cwl.begin(), cwl.end(), cwl.begin(),
                           [&](Index c) { return old_to_new[c]; });
            cwl.erase(std::remove(cwl.begin(), cwl.end(),
                                  std::numeric_limits<Index>::max()),
                      cwl.end());
        }
        return result;
    }

    const std::vector<Vertex>& spawners() const noexcept {
        return class_spawners;
    }

    void class_changed(Index cls_index, const Bitset& old_true,
                       const Bitset& new_true, Bitset& buffer) {
        buffer = old_true;
        buffer -= new_true;
        for (Lit l : buffer.ones()) {
            auto& list = classes_with_literal[l];
            list.erase(std::remove(list.begin(), list.end(), cls_index),
                       list.end());
            event_handler->literal_removed_from_class(cls_index, l);
        }
        buffer = new_true;
        buffer -= old_true;
        for (Lit l : buffer.ones()) {
            classes_with_literal[l].push_back(cls_index);
            event_handler->literal_added_to_class(cls_index, l);
        }
    }

    const std::vector<std::vector<Index>>&
    get_classes_with_literal() const noexcept {
        return classes_with_literal;
    }

    const std::vector<Index>& with_literal(Lit l) const noexcept {
        return classes_with_literal[l];
    }

    IteratorRange<std::vector<Lit>::const_iterator>
    pretend_push(Lit l, Index cls_index) {
        auto& cls = color_classes[cls_index];
        if (!cls.push_level(l)) {
            return {cls.get_trail().end(), cls.get_trail().end()};
        } else {
            return {cls.current_level_begin(), cls.get_trail().end()};
        }
    }

    void pretend_pop(Index cls_index) {
        auto& cls = color_classes[cls_index];
        cls.pop_level();
    }

    explicit ColorClasses(EventListener* event_handler, ClauseDB* all_clauses)
        : event_handler(event_handler), all_clauses(all_clauses),
          empty_class(all_clauses),
          classes_with_literal(2 * all_clauses->num_vars()) {}

    bool push_vertex(Index cindex, Lit lmin, Lit lmax) {
        SharedDBPropagator& prop = color_classes[cindex];
        bool min_true = prop.is_true(lmin);
        bool max_true = prop.is_true(lmax);
        bool pushed_min = false, pushed_max = false;
        if (min_true & max_true) {
            return true;
        }
        if (prop.is_false(lmin) || prop.is_false(lmax))
            return false;
        if (!min_true) {
            if (!prop.push_level(lmin)) {
                ConflictResolutionHandler handler{this, cindex};
                if (!prop.resolve_conflicts(handler))
                    throw UNSATError();
                return false;
            }
            if (prop.is_false(lmax)) {
                prop.pop_level();
                return false;
            }
            max_true = prop.is_true(lmax);
            pushed_min = true;
        }
        if (!max_true) {
            if (!prop.push_level(lmax)) {
                if (pushed_min) {
                    prop.pop_level();
                    prop.pop_level();
                } else {
                    ConflictResolutionHandler handler{this, cindex};
                    if (!prop.resolve_conflicts(handler))
                        throw UNSATError();
                }
                return false;
            }
            pushed_max = true;
        }
        auto new_level = prop.get_current_level() - pushed_min - pushed_max + 1;
        auto new_lits_begin = prop.level_begin(new_level);
        auto new_lits_end = prop.get_trail().end();
        for (; new_lits_begin != new_lits_end; ++new_lits_begin) {
            Lit new_lit = *new_lits_begin;
            classes_with_literal[new_lit].push_back(cindex);
            event_handler->literal_added_to_class(cindex, new_lit);
        }
        return true;
    }

    void add_class(SharedDBPropagator&& propagator, Lit n_concrete) {
        Vertex spawner;
        auto lspawner = p_spawner_from_levels(propagator, n_concrete);
        if (lspawner)
            spawner = *lspawner;
        else {
            const Lit nclit = 2 * n_concrete;
            const auto& t = propagator.get_trail();
            for (auto beg = t.begin(), end = t.end(); beg != end; ++beg) {
                Lit l1 = *beg;
                if (l1 < nclit) {
                    for (++beg; beg != end; ++beg) {
                        Lit l2 = *beg;
                        if (l2 < nclit) {
                            spawner = std::make_pair(l1, l2);
                            beg = end;
                            break;
                        }
                    }
                    --beg;
                }
            }
            if (spawner.first > spawner.second)
                std::swap(spawner.first, spawner.second);
        }
        p_add_class(std::move(propagator), spawner);
    }

    void add_class(SharedDBPropagator&& propagator, Lit, Vertex spawner) {
        p_add_class(std::move(propagator), spawner);
    }

    SharedDBPropagator& operator[](Index idx) noexcept {
        return color_classes[idx];
    }

    const SharedDBPropagator& operator[](Index idx) const noexcept {
        return color_classes[idx];
    }

    const std::vector<SharedDBPropagator>& all() const noexcept {
        return color_classes;
    }

    ClauseDB& clauses() noexcept { return *all_clauses; }
    const ClauseDB& clauses() const noexcept { return *all_clauses; }

    Index num_classes() const noexcept { return color_classes.size(); }

    /**
     * @brief Create a new propagator at level 0
     *        that includes all clauses learned so far.
     *        More efficient than creating a propagator
     *        via its constructor (copy-constructs instead).
     * @return SharedDBPropagator
     */
    SharedDBPropagator make_empty_class() {
        assert(empty_class.get_current_level() == 0);
        empty_class.incorporate_new_clauses_at_level_0();
        return empty_class;
    }

    /**
     * @brief Spawn a new color class from the given pair.
     *
     * @param spawner
     * @return true If a new color class was successfully spawned.
     *              This guarantees that the given pair is actually feasible
     *              in the strong sense, i.e., it can be extended to a complete
     * assignment.
     * @return false The given feature literal pair is infeasible.
     *               This is also reported to the event handler.
     */
    bool spawn_class(Vertex spawner) {
        assert(empty_class.get_current_level() == 0);
        empty_class.incorporate_new_clauses_at_level_0();
        if (p_spawn_class_one_is_false(spawner))
            return false;
        SharedDBPropagator new_class(empty_class);
        bool ft = new_class.is_true(spawner.first);
        bool st = new_class.is_true(spawner.second);
        if (ft | st) {
            if (ft & st) {
                event_handler->pair_definitely_feasible(spawner.first,
                                                        spawner.second);
                p_add_class(std::move(new_class), spawner);
                return true;
            }
            if (!new_class.push_level(ft ? spawner.second : spawner.first)) {
                event_handler->feature_literal_pair_infeasible(spawner.first,
                                                               spawner.second);
                new_class.resolve_or_throw();
                return false;
            }
        } else {
            if (!new_class.push_level(spawner.first)) {
                event_handler->feature_literal_infeasible(spawner.first);
                new_class.resolve_or_throw();
                return false;
            }
            if (new_class.is_false(spawner.second)) {
                event_handler->feature_literal_pair_infeasible(spawner.first,
                                                               spawner.second);
                return false;
            }
            if (!new_class.is_true(spawner.second) &&
                !new_class.push_level(spawner.second))
            {
                event_handler->feature_literal_pair_infeasible(spawner.first,
                                                               spawner.second);
                new_class.resolve_or_throw();
                return false;
            }
        }
        if (!event_handler->is_pair_known_feasible(spawner.first,
                                                   spawner.second))
        {
            if (!p_try_completion(new_class)) {
                event_handler->feature_literal_pair_infeasible(spawner.first,
                                                               spawner.second);
                return false;
            }
            event_handler->pair_definitely_feasible(spawner.first,
                                                    spawner.second);
        }
        p_add_class(std::move(new_class), spawner);
        return true;
    }

    void initialize_feasibilities() {
        empty_class.incorporate_or_throw();
        const auto nl = 2 * all_clauses->num_vars();
        for (Lit l = 0; l < nl; ++l) {
            if (empty_class.is_true(l))
                continue;
            if (empty_class.is_false(l)) {
                event_handler->feature_literal_infeasible(l);
                continue;
            }
            if (!empty_class.push_level(l)) {
                event_handler->feature_literal_infeasible(l);
                empty_class.resolve_or_throw();
                continue;
            }
            const auto& trail = empty_class.get_trail();
            for (std::size_t ti = 1, s = trail.size(); ti < s; ++ti) {
                Lit f = lit::negate(trail[ti]);
                event_handler->feature_literal_pair_infeasible(l, f);
            }
            empty_class.pop_level();
        }
        empty_class.incorporate_or_throw();
    }

    std::size_t total_memory_usage() const noexcept {
        std::size_t total_prop_usage = 0;
        for (const SharedDBPropagator& c : color_classes) {
            total_prop_usage += c.total_memory_usage();
        }
        std::size_t total_with_literal_usage = 0;
        for (const auto& v : classes_with_literal) {
            total_with_literal_usage += v.capacity() * sizeof(Index);
        }
        return empty_class.total_memory_usage() + total_prop_usage +
               (color_classes.capacity() - color_classes.size()) *
                   sizeof(SharedDBPropagator) +
               classes_with_literal.capacity() * sizeof(std::vector<Index>) +
               total_with_literal_usage +
               class_spawners.capacity() * sizeof(Vertex) +
               m_phases.capacity() / CHAR_BIT +
               m_decisions.capacity() * sizeof(Lit) + sizeof(ColorClasses);
    }

    /* void print_memory_stats() const {
        std::size_t total_prop_usage = 0;
        std::size_t max_idx = 0, i = 0, max_usage = 0;
        for (const SharedDBPropagator& c : color_classes) {
            std::size_t usage = c.total_memory_usage();
            total_prop_usage += usage;
            if (usage > max_usage) {
                max_usage = usage;
                max_idx = i;
            }
            ++i;
        }
        std::size_t total_with_literal_usage =
            classes_with_literal.capacity() * sizeof(std::vector<Index>);
        for (const auto& v : classes_with_literal) {
            total_with_literal_usage += v.capacity() * sizeof(Index);
        }
        std::cout << "    Propagators " << mibibytes(total_prop_usage) << " ("
                  << (mibibytes(total_prop_usage) / color_classes.size())
                  << " MiB per propagator)\n";
        std::cout << "    Largest propagator: \n";
        color_classes[max_idx].print_memory_stats();
        std::cout << "    Empty class "
                  << mibibytes(empty_class.total_memory_usage()) << " MiB\n";
        std::cout << "    CWL lists   " << mibibytes(total_with_literal_usage)
                  << " MiB\n";
    } */

    void reset_coloring() {
        empty_class.incorporate_new_clauses_at_level_0();
        color_classes.clear();
        for (auto& cc : classes_with_literal) {
            cc.clear();
        }
        class_spawners.clear();
    }

    template <typename FullAssignment>
    void incorporate_colors(const std::vector<FullAssignment>& colors,
                            const std::vector<Vertex>& spawners) {
        for (std::size_t i = 0, s = colors.size(); i < s; ++i) {
            SharedDBPropagator propagator{empty_class};
            propagator.incorporate_assignment(colors[i]);
            p_add_class(std::move(propagator), spawners[i]);
        }
    }
};

struct QueueVertexEntry {
    Lit lmin, lmax;
    std::uint32_t class_count;
    Bitset available_classes;

    QueueVertexEntry(Lit lmin, Lit lmax, std::uint32_t class_count,
                     Bitset available_classes)
        : lmin(lmin), lmax(lmax), class_count(class_count),
          available_classes(std::move(available_classes)) {}

    bool class_unavailable(Index cls) noexcept {
        if (available_classes[cls]) {
            available_classes[cls] = false;
            --class_count;
            return true;
        }
        return false;
    }

    bool class_available(Index cls) noexcept {
        if (!available_classes[cls]) {
            available_classes[cls] = true;
            ++class_count;
            return true;
        }
        return false;
    }

    std::size_t total_memory_usage() const noexcept {
        return sizeof(QueueVertexEntry) + available_classes.bytes_used();
    }
};

struct EntryComesBefore {
    bool operator()(const QueueVertexEntry& e1,
                    const QueueVertexEntry& e2) const noexcept {
        return e1.class_count < e2.class_count;
    }
};

template <typename IndexMap, typename Compare = EntryComesBefore>
class IndexedVertexHeap {
  public:
    QueueVertexEntry& operator[](std::size_t s) noexcept {
        assert(s < m_entries.size());
        return m_entries[s];
    }

    const QueueVertexEntry& operator[](std::size_t s) const noexcept {
        assert(s < m_entries.size());
        return m_entries[s];
    }

    Index index_of(Lit l1, Lit l2) const noexcept {
        return ordered_index_of((std::min)(l1, l2), (std::max)(l1, l2));
    }

    Index ordered_index_of(Lit lmin, Lit lmax) const noexcept {
        return m_indices.index_of(lmin, lmax);
    }

    std::vector<QueueVertexEntry>& entries() noexcept { return m_entries; }

    const std::vector<QueueVertexEntry>& entries() const noexcept {
        return m_entries;
    }

    QueueVertexEntry& top() noexcept { return m_entries.front(); }

    const QueueVertexEntry& top() const noexcept { return m_entries.front(); }

    bool empty() const noexcept { return m_entries.empty(); }

    void pop() {
        assert(!empty());
        p_swap(0, m_entries.size() - 1);
        p_pop_back();
        lowered_priority(0); // works for empty queue as well
    }

    std::size_t size() const noexcept { return m_entries.size(); }

    template <typename... Args> void push(Args&&... args) {
        Index idx = m_entries.size();
        m_entries.emplace_back(std::forward<Args>(args)...);
        m_indices.add(m_entries.back().lmin, m_entries.back().lmax, idx);
        increased_priority(idx);
    }

    void lowered_priority(std::size_t idx) noexcept {
        QueueVertexEntry* cur = &m_entries[idx];
        std::size_t cpos = idx, hsize = m_entries.size();
        while ((cpos = (cpos << 1) + 1) < hsize) {
            QueueVertexEntry* c1 = &m_entries[cpos - 1];
            QueueVertexEntry* c2 = &m_entries[cpos];
            if (m_comp(*c1, *c2)) {
                if (m_comp(*c1, *cur)) {
                    p_swap(cpos - 1, idx);
                    idx = cpos - 1;
                    cur = c1;
                } else {
                    return;
                }
            } else {
                if (m_comp(*c2, *cur)) {
                    p_swap(cpos, idx);
                    idx = cpos;
                    cur = c2;
                } else {
                    return;
                }
            }
        }
        if (cpos == hsize) {
            QueueVertexEntry* c1 = &m_entries[cpos - 1];
            if (m_comp(*c1, *cur)) {
                p_swap(idx, cpos - 1);
            }
        }
    }

    void increased_priority(std::size_t idx) noexcept {
        QueueVertexEntry* cur = &m_entries[idx];
        for (;;) {
            std::size_t par = idx >> 1;
            QueueVertexEntry* p = &m_entries[par];
            if (m_comp(*cur, *p)) {
                p_swap(par, idx);
                idx = par;
                cur = p;
            } else {
                return;
            }
        }
    }

    void clear() {
        m_entries.clear();
        m_indices.clear();
    }

    template <typename... CompareArgs>
    IndexedVertexHeap(std::size_t n_concrete, CompareArgs&&... c)
        : m_entries(), m_indices(n_concrete),
          m_comp(std::forward<CompareArgs>(c)...) {}

    std::size_t total_memory_usage() const noexcept {
        std::size_t total_queue_usage = 0;
        for (const auto& e : m_entries)
            total_queue_usage += e.total_memory_usage();
        return m_indices.total_memory_usage() + total_queue_usage +
               (m_entries.capacity() - m_entries.size()) *
                   sizeof(QueueVertexEntry);
    }

  private:
    std::vector<QueueVertexEntry> m_entries;
    IndexMap m_indices;
    Compare m_comp;

    void p_pop_back() noexcept {
        const auto& b = m_entries.back();
        m_indices.remove(b.lmin, b.lmax);
        m_entries.pop_back();
    }

    void p_swap(std::size_t i1, std::size_t i2) noexcept {
        auto& e1 = m_entries[i1];
        auto& e2 = m_entries[i2];
        m_indices.swap(e1.lmin, e1.lmax, e2.lmin, e2.lmax);
        std::swap(e1, e2);
    }
};

template <typename T> class LiteralPairMatrix {
  public:
    explicit LiteralPairMatrix(std::size_t n_concrete, const T& init)
        : row_length(n_concrete * 2),
          buffer(std::make_unique<T[]>(row_length * (n_concrete - 1))) {
        set_all(init);
    }

    void set_all(const T& value) noexcept {
        std::size_t n_concrete = row_length / 2;
        std::fill_n(buffer.get(), row_length * (n_concrete - 1), value);
    }

    std::size_t offset(Lit lmin, Lit lmax) const noexcept {
        std::size_t pairs_below = lmin >> 1;
        std::size_t row_offset =
            2 * pairs_below * (row_length - pairs_below - 1);
        std::size_t even_odd_offset = -(lmin & 1) & (row_length - lmin);
        return row_offset + even_odd_offset + lmax - lmin - 2;
    }

    T& operator()(Lit lmin, Lit lmax) noexcept {
        return buffer[offset(lmin, lmax)];
    }

    const T& operator()(Lit lmin, Lit lmax) const noexcept {
        return buffer[offset(lmin, lmax)];
    }

    template <typename Callable /*(Lit lmin, Lit lmax, T&)*/>
    void iterate_all(Callable&& callable) {
        T* current = buffer.get();
        for (Lit vmin = 0, nc = row_length / 2; vmin < nc; ++vmin) {
            Lit lmin = 2 * vmin;
            for (Lit vmax = vmin + 1; vmax < nc; ++vmax, current += 2) {
                Lit lmax = 2 * vmax;
                callable(lmin, lmax, current[0]);
                callable(lmin, lmax + 1, current[1]);
            }
            ++lmin;
            for (Lit vmax = vmin + 1; vmax < nc; ++vmax, current += 2) {
                Lit lmax = 2 * vmax;
                callable(lmin, lmax, current[0]);
                callable(lmin, lmax + 1, current[1]);
            }
        }
    }

    /**
     * @brief Iterates the 'rows' of the matrix.
     *        Each row corresponds to the entries
     *        for some literal row_lit and has two entries
     *        for each variable strictly greater than the
     *        one of row_lit.
     *
     * @tparam RowCallbacks
     * @param row_out
     */
    template <typename RowCallback /*(Lit row_lit, Lit first_other, T* row)*/>
    void iterate_rows(RowCallback&& row_out) {
        T* current = buffer.get();
        for (Lit vmin = 0, nc = row_length / 2; vmin < nc; ++vmin) {
            Lit lmin = 2 * vmin;
            std::forward<RowCallback>(row_out)(lmin, lmin + 2, current);
            current += 2 * (nc - vmin - 1);
            std::forward<RowCallback>(row_out)(lmin + 1, lmin + 2, current);
            current += 2 * (nc - vmin - 1);
        }
    }

    template <typename Callable /*(Lit lmin, Lit lmax, const T&)*/>
    void iterate_all(Callable&& callable) const noexcept {
        const T* current = buffer.get();
        for (Lit vmin = 0, nc = row_length / 2; vmin < nc; ++vmin) {
            Lit lmin = 2 * vmin;
            for (Lit vmax = vmin + 1; vmax < nc; ++vmax, current += 2) {
                Lit lmax = 2 * vmax;
                callable(lmin, lmax, current[0]);
                callable(lmin, lmax + 1, current[1]);
            }
            ++lmin;
            for (Lit vmax = vmin + 1; vmax < nc; ++vmax, current += 2) {
                Lit lmax = 2 * vmax;
                callable(lmin, lmax, current[0]);
                callable(lmin, lmax + 1, current[1]);
            }
        }
    }

    std::size_t total_memory_usage() const noexcept {
        const std::size_t n_concrete = row_length / 2;
        return sizeof(LiteralPairMatrix) +
               (row_length * (n_concrete - 1) * sizeof(T));
    }

  private:
    std::size_t row_length;
    std::unique_ptr<T[]> buffer;
};

class MatrixIndexMap {
  public:
    explicit MatrixIndexMap(std::size_t n_concrete)
        : m_matrix(n_concrete, NIL) {}

    void swap(Lit lmin1, Lit lmax1, Lit lmin2, Lit lmax2) {
        std::swap(m_matrix(lmin1, lmax1), m_matrix(lmin2, lmax2));
    }

    void remove(Lit lmin, Lit lmax) { m_matrix(lmin, lmax) = NIL; }

    void add(Lit lmin, Lit lmax, Index index) { m_matrix(lmin, lmax) = index; }

    Index index_of(Lit lmin, Lit lmax) const noexcept {
        return m_matrix(lmin, lmax);
    }

    std::size_t total_memory_usage() const noexcept {
        return m_matrix.total_memory_usage();
    }

    void clear() noexcept { m_matrix.set_all(NIL); }

  private:
    LiteralPairMatrix<Index> m_matrix;
};

class ColoringHeuristicSolver {
  public:
    ColoringHeuristicSolver(ClauseDB* all_clauses, Lit n_concrete,
                            ThreadGroup<void>* thread_pool)
        : all_clauses(all_clauses), n_concrete(n_concrete),
          color_classes(this, all_clauses),
          internal_inf_map(std::in_place, n_concrete),
          inf_map(&*internal_inf_map),
          vertex_queue(n_concrete, EntryComesBefore{}),
          class_completer(n_concrete, all_clauses->num_vars(), this),
          explicit_partners_of(2 * n_concrete, std::vector<Lit>{}),
          m_old_true(2 * all_clauses->num_vars(), false),
          m_new_true(2 * all_clauses->num_vars(), false),
          m_buffer(2 * all_clauses->num_vars(), false),
          m_bitset_buffer(thread_pool) {}

    ColoringHeuristicSolver(ClauseDB* all_clauses, Lit n_concrete,
                            ThreadGroup<void>* thread_pool,
                            PairInfeasibilityMap* infmap)
        : all_clauses(all_clauses), n_concrete(n_concrete),
          color_classes(this, all_clauses), internal_inf_map(std::nullopt),
          inf_map(infmap), vertex_queue(n_concrete, EntryComesBefore{}),
          class_completer(n_concrete, all_clauses->num_vars(), this),
          explicit_partners_of(2 * n_concrete, std::vector<Lit>{}),
          m_old_true(2 * all_clauses->num_vars(), false),
          m_new_true(2 * all_clauses->num_vars(), false),
          m_buffer(2 * all_clauses->num_vars(), false),
          m_bitset_buffer(thread_pool) {}

    void initialize_feasibilities() {
        color_classes.initialize_feasibilities();
    }

    template <typename RNG>
    Vertex random_from_class(Index index, RNG& rng) const {
        const auto& t = color_classes[index].get_trail();
        const auto nclit = 2 * n_concrete;
        Lit l1, l2;
        std::uniform_int_distribution<std::size_t> idx_gen(0, t.size() - 1);
        do {
            l1 = t[idx_gen(rng)];
        } while (l1 >= nclit);
        do {
            l2 = t[idx_gen(rng)];
        } while (l2 >= nclit || l1 == l2);
        return {std::min(l1, l2), std::max(l1, l2)};
    }

    void extract_feasibilities() {
        if (extracted_feasibilities)
            return;

#ifdef SAMMY_CUDA_SUPPORTED
        if (n_concrete > 2048 && color_classes.all().size() > 512 &&
            should_use_cuda())
        {
            inf_map->cuda_incorporate_complete_classes(
                literals_in_class, m_bitset_buffer.thread_group());
            extracted_feasibilities = true;
            return;
        }
#endif

        inf_map->incorporate_complete_classes(literals_in_class,
                                              m_bitset_buffer.thread_group());
        extracted_feasibilities = true;
    }

    std::vector<Vertex> extract_coloring_order() {
        if (!extracted_feasibilities)
            extract_feasibilities();
        coloring_order.erase(
            std::remove_if(
                coloring_order.begin(), coloring_order.end(),
                [this](Vertex v) { return (*inf_map)(v.first, v.second); }),
            coloring_order.end());
        EdgeSet elements;
        coloring_order.erase(
            std::remove_if(
                coloring_order.begin(), coloring_order.end(),
                [&](Vertex v) { return !elements.insert(v).second; }),
            coloring_order.end());
        return coloring_order;
    }

    void feature_literal_infeasible(Lit l) {
        if (l < 2 * n_concrete) {
            inf_map->literal_infeasible(l);
        }
    }

    void feature_literal_pair_infeasible(Lit l1, Lit l2) {
        if ((std::max)(l1, l2) < 2 * n_concrete) {
            inf_map->literal_pair_infeasible(l1, l2);
        }
    }

    bool is_covered(Lit lmin, Lit lmax) const noexcept {
        if ((*inf_map)(lmin, lmax))
            return true;
        for (Index cls_index : color_classes.with_literal(lmin)) {
            if (literals_in_class[cls_index][lmax])
                return true;
        }
        return false;
    }

    Index find_covering_class(Lit lmin, Lit lmax) const {
        if ((*inf_map)(lmin, lmax))
            return NIL;
        for (Index cls_index : color_classes.with_literal(lmin)) {
            if (literals_in_class[cls_index][lmax])
                return cls_index;
        }
        throw std::out_of_range(
            "Vertex (" + std::to_string(lmin) + ", " + std::to_string(lmax) +
            ") is not covered and not marked as infeasible!");
    }

    void find_color_class_for_first_prefer_half(Lit lmin, Lit lmax) {
        // check if already colored or infeasible
        if (is_covered(lmin, lmax))
            return;
        // need to copy (concurrent modification!)
        candidates = color_classes.with_literal(lmin);
        // check classes containing lmin
        if (p_find_first_suitable_candidate_one_new(lmin, lmax, lmax))
            return;
        // check classes containing lmax
        candidates = color_classes.with_literal(lmax);
        if (p_find_first_suitable_candidate_one_new(lmin, lmax, lmin))
            return;
        // check all classes
        if (p_find_first_suitable_candidate_from_all(lmin, lmax))
            return;
        // if all else fails, create new class
        if (color_classes.spawn_class(Vertex(lmin, lmax))) {
            coloring_order.emplace_back(lmin, lmax);
        }
    }

    template <typename Callable /*(Lit lmin, Lit lmax)*/>
    void iterate_all_uncovered(Callable&& callable) const {
        Lit nclit = 2 * n_concrete;

#ifdef SAMMY_CUDA_SUPPORTED
        if (nclit > 2048 && color_classes.all().size() > 512 &&
            should_use_cuda())
        {
            try {
                cuda_iterate_all_uncovered(
                    literals_in_class, color_classes.get_classes_with_literal(),
                    inf_map, callable);
                return;
            } catch (const CUDAError& err) {
                std::cerr << "Not using CUDA because of an error: "
                          << err.what() << "\n";
                had_cuda_error(err);
            }
        }
#endif

        Bitset row_uncolored(nclit, true);
        for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
            row_uncolored.set();
            row_uncolored ^= (*inf_map)[lmin];
            const auto& with_lit = color_classes.with_literal(lmin);
            sammy::bitwise_filter(
                m_bitset_buffer, row_uncolored,
                p_make_bitset_transform_iterator(with_lit.begin()),
                p_make_bitset_transform_iterator(with_lit.end()));
            for (Lit lmax : row_uncolored.ones_from(lmin + 1)) {
                (std::forward<Callable>(callable))(lmin, lmax);
            }
        }
    }

    template <typename Callable /*(Lit lmin, Lit lmax)*/>
    void iterate_all_uniquely_covered(Callable&& callable) {
        Lit nclit = 2 * n_concrete;
        Bitset multi_covered(nclit, false);
        Bitset prev_covered(nclit, false);
        Bitset tmp(nclit, false);
        for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
            multi_covered.reset();
            prev_covered.reset();
            const auto& with_lit = color_classes.with_literal(lmin);
            for (Index cci : with_lit) {
                const auto& ccbits = literals_in_class[cci];
                tmp = prev_covered;
                tmp &= ccbits;
                multi_covered |= tmp;
                prev_covered |= ccbits;
            }
            tmp = multi_covered;
            tmp.flip();
            tmp -= (*inf_map)[lmin];
            for (Lit lmax : tmp.ones_from(lmin + 1)) {
                (std::forward<Callable>(callable))(lmin, lmax);
            }
        }
    }

    template <typename FullAssignment>
    void replace_colors(const std::vector<FullAssignment>& colors,
                        const std::vector<Vertex>& spawners) {
        reset_coloring();
        color_classes.incorporate_colors(colors, spawners);
    }

    PartialSolution get_partial_solution() const {
        return PartialSolution{all_clauses->num_vars(), inf_map,
                               all_classes().begin(), all_classes().end()};
    }

    LiteralPairMatrix<std::uint16_t> saturated_coverage_count_matrix() const {
        std::vector<std::uint16_t> coverage_buf(2 * all_clauses->num_vars(), 0);
        LiteralPairMatrix<std::uint16_t> result(n_concrete, 0);
        std::uint16_t* cv_buf_beg = &coverage_buf[2];
        std::uint16_t* cv_buf_end = coverage_buf.data() + (2 * n_concrete);
        result.iterate_rows([&](Lit row_lit, Lit /*next*/, std::uint16_t* row) {
            if (color_classes.all().size() <=
                std::numeric_limits<std::uint16_t>::max())
            {
                for (Index cci : color_classes.with_literal(row_lit)) {
                    const SharedDBPropagator& cc = color_classes[cci];
                    for (Lit l : cc.get_trail())
                        ++coverage_buf[l];
                }
            } else {
                for (Index cci : color_classes.with_literal(row_lit)) {
                    const SharedDBPropagator& cc = color_classes[cci];
                    for (Lit l : cc.get_trail())
                        p_saturated_add(coverage_buf[l], 1);
                }
            }
            for (std::uint16_t* cv_cur = cv_buf_beg; cv_cur != cv_buf_end;
                 ++cv_cur, ++row)
            {
                p_saturated_add(*row, *cv_cur);
                *cv_cur = 0;
            }
            if (row_lit & 1)
                cv_buf_beg += 2;
        });
        return result;
    }

    void find_color_class_for(Lit lmin, Lit lmax) {
        find_color_class_for_first_prefer_half(lmin, lmax);
    }

    bool complete_class(Index index) {
        SharedDBPropagator& prop = color_classes[index];
        if (prop.get_trail().size() == all_clauses->num_vars())
            return true;
        p_prop_to_buffer(prop, m_old_true);
        bool compl_result = class_completer.complete_class(prop);
        p_prop_to_buffer(prop, m_new_true);
        color_classes.class_changed(index, m_old_true, m_new_true, m_buffer);
        return compl_result;
    }

    bool complete_classes() {
        Index nc = color_classes.all().size();
        for (Index i = 0; i < nc; ++i) {
            if (!complete_class(i)) {
                if (!quiet) {
                    std::cout << "Could not complete class " << i
                              << " without rolling back a decision!\n";
                }
                return false;
            }
        }
        return true;
    }

    void add_greedy_false_class() {
        SharedDBPropagator prop = color_classes.make_empty_class();
        GreedySAT gsat(&prop, n_concrete, PreferFalse{});
        if (gsat.solve()) {
            add_color_class(std::move(prop));
        }
    }

    void color_nonlazy(bool use_initialization = true) {
        if (use_initialization) {
            add_greedy_false_class();
            auto initial = default_initial_vertices();
            p_bulk_add_vertices(initial.cbegin(), initial.cend());
            color_vertices_in_queue();
        }
        for (;;) {
            add_all_uncolored_to_queue();
            color_vertices_in_queue();
            if (everything_covered()) {
                if (complete_classes())
                    return;
            }
        }
    }

    void color_nonlazy(const std::vector<Vertex>& initial_vertices) {
        for (Vertex v : initial_vertices) {
            add_vertex_to_queue(v.first, v.second);
        }
        color_vertices_in_queue();
        for (;;) {
            add_all_uncolored_to_queue();
            color_vertices_in_queue();
            if (everything_covered()) {
                if (complete_classes())
                    return;
            }
        }
    }

    void color_lazy(const std::vector<Vertex>& initial_vertices_ = {}) {
        bool had_initial_vertices = true;
        if (initial_vertices_.empty() && coloring_order.empty()) {
            add_greedy_false_class();
            had_initial_vertices = false;
        }
        const std::vector<Vertex>& initial_vertices =
            !had_initial_vertices ? default_initial_vertices()
                                  : initial_vertices_;
        for (Vertex v : initial_vertices) {
            add_vertex_to_queue(v.first, v.second);
        }
        color_vertices_in_queue();
        std::size_t lazy_col_explicit = coloring_order.size();
        std::size_t lazy_col_next = p_lazy_color_next(lazy_col_explicit);
        if (had_initial_vertices &&
            lazy_col_explicit != color_classes.all().size())
        {
            for (Vertex v : initial_vertices) {
                if (!is_covered(v.first, v.second)) {
                    // may happen due to conflict resolution
                    find_color_class_for(v.first, v.second);
                }
            }
        }
        if (!quiet) {
            std::cout << "Colored initial " << lazy_col_explicit
                      << " vertices (" << color_classes.all().size()
                      << " classes)\n";
        }
        // print_memory_stats();
        for (;;) {
            std::size_t n_unc = num_uncovered();
            if (!quiet) {
                std::cout << "Uncovered: " << n_unc << "...\n";
            }
            if (n_unc == 0) {
                if (complete_classes()) {
                    return;
                }
                continue;
            }
            if (n_unc <= 2 * lazy_col_next) {
                add_all_uncolored_to_queue();
                lazy_col_explicit += n_unc;
            } else {
                lazy_col_explicit += lazy_sample_to_queue(lazy_col_next, n_unc);
            }
            // print_memory_stats();
            color_vertices_in_queue();
            lazy_col_next = p_lazy_color_next(lazy_col_explicit);
            if (!quiet) {
                std::cout << "Colored " << lazy_col_explicit
                          << " explicit vertices ("
                          << color_classes.all().size() << " classes) ...\n";
            }
        }
    }

    std::size_t lazy_sample_to_queue(std::size_t approx_vertices,
                                     std::size_t n_unc) {
        double sample_rate = double(approx_vertices) / n_unc;
        auto& rng = sammy::rng();
        std::geometric_distribution<std::size_t> skip_sample_dist(sample_rate);
        std::size_t skip = skip_sample_dist(rng);
        std::size_t count = 0;
        std::vector<Vertex> vbuf;
        iterate_all_uncovered([&](Lit lmin, Lit lmax) {
            if (!skip) {
                skip = skip_sample_dist(rng);
                vbuf.emplace_back(lmin, lmax);
                ++count;
            } else {
                --skip;
            }
        });
        p_bulk_add_vertices(vbuf.begin(), vbuf.end());
        return count;
    }

    std::vector<Vertex> default_initial_vertices() {
        std::vector<Vertex> result;
        result.reserve(n_concrete);
        std::vector<Lit> v2;
        v2.reserve(n_concrete);
        for (Lit l = 0; l < n_concrete; ++l) {
            v2.push_back(lit::positive_lit(l));
        }
        std::shuffle(v2.begin(), v2.end(), sammy::rng());
        SharedDBPropagator prop = color_classes.make_empty_class();
        for (Lit v1 = 0; v1 < n_concrete; ++v1) {
            Lit l1 = lit::positive_lit(v1);
            if (!prop.is_open(l1))
                continue;
            if (!prop.push_level(l1)) {
                prop.resolve_or_throw();
                continue;
            }
            auto start_from = v2.begin() + v1;
            for (; start_from != v2.end(); ++start_from) {
                Lit o = *start_from;
                if (!prop.is_open(o))
                    continue;
                bool pres = prop.push_level(o);
                prop.pop_level();
                if (!pres)
                    continue;
                break;
            }
            if (start_from == v2.end()) {
                for (start_from = v2.begin(); start_from != v2.end();
                     ++start_from)
                {
                    Lit o = *start_from;
                    if (!prop.is_open(o))
                        continue;
                    bool pres = prop.push_level(o);
                    prop.pop_level();
                    if (!pres)
                        continue;
                    break;
                }
            }
            if (start_from != v2.end()) {
                result.emplace_back(std::min(l1, *start_from),
                                    std::max(l1, *start_from));
            }
            prop.pop_level();
        }
        return result;
    }

    const std::vector<Vertex>& class_spawners() const noexcept {
        return color_classes.spawners();
    }

    void new_color_class(Index index) {
        for (auto& e : vertex_queue.entries()) {
            e.class_count += 1;
            e.available_classes.push_back(true);
        }
        literals_in_class.push_back(Bitset(2 * n_concrete, false));
        const auto& cc = color_classes[index];
        for (Lit lpos : cc.get_trail()) {
            literal_added_to_class(index, lpos);
        }
    }

    void literal_added_to_class(Index cindex, Lit l) {
        if (l >= 2 * n_concrete)
            return;
        literals_in_class[cindex][l] = true;
        l = lit::negate(l);
        auto& partners = explicit_partners_of[l];
        auto pend =
            std::remove_if(partners.begin(), partners.end(), [&](Lit o) {
                Index idx = vertex_queue.index_of(l, o);
                if (idx == NIL)
                    return true;
                if (vertex_queue[idx].class_unavailable(cindex)) {
                    vertex_queue.increased_priority(idx);
                }
                return false;
            });
        partners.erase(pend, partners.end());
    }

    void literal_removed_from_class(Index cindex, Lit l) {
        if (l >= 2 * n_concrete)
            return;
        literals_in_class[cindex][l] = false;
        SharedDBPropagator& cc = color_classes[cindex];
        l = lit::negate(l);
        auto& partners = explicit_partners_of[l];
        auto pend =
            std::remove_if(partners.begin(), partners.end(), [&](Lit o) {
                Index idx = vertex_queue.index_of(l, o);
                if (idx == NIL)
                    return true;
                if (!cc.is_false(o) &&
                    vertex_queue[idx].class_available(cindex)) {
                    vertex_queue.lowered_priority(idx);
                }
                return false;
            });
        partners.erase(pend, partners.end());
    }

    bool is_pair_known_feasible(Lit lmin, Lit lmax) const noexcept {
        return inf_map->is_definitely_feasible(lmin, lmax);
    }

    void pair_definitely_feasible(Lit lmin, Lit lmax) noexcept {
        inf_map->set_definitely_feasible(lmin, lmax);
    }

    void color_vertices_in_queue() {
        while (!vertex_queue.empty()) {
            QueueVertexEntry ve = std::move(vertex_queue.top());
            vertex_queue.pop();
            if (ve.class_count == 0) {
                if (color_classes.spawn_class(Vertex(ve.lmin, ve.lmax))) {
                    coloring_order.emplace_back(ve.lmin, ve.lmax);
                }
            } else {
                find_color_class_for(ve.lmin, ve.lmax);
            }
        }
    }

    void initialize_with_start_vertices(const std::vector<Vertex>& vertices) {
        for (Vertex v : vertices) {
            find_color_class_for(v.first, v.second);
        }
    }

    void add_vertex_to_queue(Lit lmin, Lit lmax) {
        const auto& cc_all = color_classes.all();
        auto s = cc_all.size();
        std::uint32_t count = s;
        Bitset bset(s, true);
        for (Index i = 0, s = cc_all.size(); i < s; ++i) {
            const auto& cf = cc_all[i];
            if (cf.is_false(lmin) || cf.is_false(lmax)) {
                bset[i] = false;
                --count;
            }
        }
        vertex_queue.push(lmin, lmax, count, std::move(bset));
    }

    template <typename Iterator>
    void add_vertices_to_queue(Iterator begin, Iterator end) {
        p_bulk_add_vertices(begin, end);
    }

    void add_all_uncolored_to_queue() {
        std::vector<Vertex> vbuf;
        iterate_all_uncovered(
            [&](Lit lmin, Lit lmax) { vbuf.emplace_back(lmin, lmax); });
        p_bulk_add_vertices(vbuf.begin(), vbuf.end());
    }

    void add_color_class(const SharedDBPropagator& cc) {
        SharedDBPropagator prop(cc);
        color_classes.add_class(std::move(prop), n_concrete);
    }

    void add_color_class(const SharedDBPropagator& cc, Vertex spawner) {
        SharedDBPropagator prop(cc);
        color_classes.add_class(std::move(prop), n_concrete, spawner);
    }

    void add_color_class(SharedDBPropagator&& cc, Vertex spawner) {
        color_classes.add_class(std::move(cc), n_concrete, spawner);
    }

    void add_color_class(SharedDBPropagator&& cc) {
        color_classes.add_class(std::move(cc), n_concrete);
    }

    const std::vector<SharedDBPropagator>& all_classes() const noexcept {
        return color_classes.all();
    }

    bool everything_covered() const {
        bool res = true;
        iterate_all_uncovered([&](Lit, Lit) { res = false; });
        return res;
    }

    std::size_t num_uncovered() const {
        Lit nclit = 2 * n_concrete;

#ifdef SAMMY_CUDA_SUPPORTED
        if (nclit > 2048 && color_classes.all().size() > 512 &&
            should_use_cuda())
        {
            try {
                auto result = cuda_count_uncovered(
                    literals_in_class, color_classes.get_classes_with_literal(),
                    inf_map);
                return result;
            } catch (const CUDAError& err) {
                std::cerr << "Not using CUDA because of an error: "
                          << err.what() << "\n";
                had_cuda_error(err);
            }
        }
#endif

        Bitset row_uncolored(nclit, true);
        Bitset prv_rows(nclit, true);
        std::size_t nu = 0;
        for (Lit lmin = 0; lmin < nclit; ++lmin) {
            prv_rows[lmin] = false;
            row_uncolored = prv_rows;
            row_uncolored -= (*inf_map)[lmin];
            const auto& with_lit = color_classes.with_literal(lmin);
            sammy::bitwise_filter(
                m_bitset_buffer, row_uncolored,
                p_make_bitset_transform_iterator(with_lit.begin()),
                p_make_bitset_transform_iterator(with_lit.end()));
            nu += row_uncolored.count();
        }
        return nu;
    }

    std::size_t total_memory_usage() const noexcept {
        return color_classes.total_memory_usage() +
               (internal_inf_map ? internal_inf_map->total_memory_usage() : 0) +
               vertex_queue.total_memory_usage() +
               class_completer.total_memory_usage() +
               p_mem_usage_explicit_partners_of() +
               candidates.capacity() * sizeof(Index) +
               p_mem_usage_literals_in_class(literals_in_class) +
               m_old_true.bytes_used() + m_new_true.bytes_used() +
               m_buffer.bytes_used() + sizeof(ColoringHeuristicSolver);
    }

    /**
     * @brief Reset the coloring (dropping all color classes and related
     * information). Retains quite a bit of internal memory (vector capacities
     * etc).
     */
    void reset_coloring() {
        color_classes.reset_coloring();
        vertex_queue.clear();
        for (auto& e : explicit_partners_of)
            e.clear();
        candidates.clear();
        literals_in_class.clear();
        coloring_order.clear();
    }

    /**
     * @brief Remove color classes given by the indices.
     * Deletes the coloring order information.
     * Returns the color classes that were removed.
     */
    std::vector<SharedDBPropagator>
    remove_color_classes(std::vector<std::size_t> classes) {
        if (!vertex_queue.empty()) {
            throw std::logic_error(
                "remove_color_classes called with non-empty vertex queue!");
        }
        std::for_each(explicit_partners_of.begin(), explicit_partners_of.end(),
                      [](CVec& v) { v.clear(); });
        std::sort(classes.begin(), classes.end());
        candidates.clear();
        coloring_order.clear();
        auto result = color_classes.remove_classes(classes);
        literals_in_class.erase(remove_indices(literals_in_class.begin(),
                                               literals_in_class.end(),
                                               classes.begin(), classes.end()),
                                literals_in_class.end());
        return result;
    }

    /**
     * @brief Get the resulting sample as list of lists of internal literals.
     *
     * @return std::vector<std::vector<Lit>>
     */
    std::vector<std::vector<Lit>> internal_solution(bool only_concrete) const {
        std::vector<std::vector<Lit>> result;
        for (const auto& cc : color_classes.all()) {
            result.emplace_back(cc.get_trail());
        }
        if (only_concrete) {
            const Lit nclit = 2 * n_concrete;
            for (auto& cc : result) {
                cc.erase(std::remove_if(cc.begin(), cc.end(),
                                        [&](Lit l) { return l >= nclit; }),
                         cc.end());
            }
        }
        return result;
    }

    /**
     * @brief Get the resulting sample as list of lists of external literals.
     * @return std::vector<ExternalClause>
     */
    std::vector<std::vector<ExternalLit>>
    external_solution(bool only_concrete) const {
        std::vector<std::vector<ExternalLit>> result;
        const Lit nclit = 2 * n_concrete;
        for (const auto& cc : color_classes.all()) {
            result.emplace_back();
            auto& r = result.back();
            for (const Lit l : cc.get_trail()) {
                if (!only_concrete || l < nclit) {
                    r.push_back(lit::externalize(l));
                }
            }
        }
        return result;
    }

    void set_quiet(bool q) noexcept { quiet = q; }

  private:
    struct ClassIndexToBitset {
        const Bitset& operator()(Index cls_index) const noexcept {
            return that->literals_in_class[cls_index];
        }

        const ColoringHeuristicSolver* that;
    };

    template <typename IndexIterator>
    boost::transform_iterator<ClassIndexToBitset, IndexIterator, const Bitset&>
    p_make_bitset_transform_iterator(IndexIterator iter) const {
        ClassIndexToBitset trans{this};
        return {iter, trans};
    }

    static void p_saturated_add(std::uint16_t& x, std::uint16_t v) {
        if (x > std::numeric_limits<std::uint16_t>::max() - v) {
            x = std::numeric_limits<std::uint16_t>::max();
        } else {
            x += v;
        }
    }

    std::size_t
    p_mem_usage_literals_in_class(const std::vector<Bitset>& v) const noexcept {
        return std::transform_reduce(
            v.begin(), v.end(), std::size_t(0), std::plus<>{},
            [](const Bitset& b) { return b.bytes_used(); });
    }

    std::size_t p_mem_usage_explicit_partners_of() const noexcept {
        std::size_t element_sum = 0;
        for (const auto& v : explicit_partners_of) {
            element_sum += v.capacity() * sizeof(Lit);
        }
        return explicit_partners_of.capacity() * sizeof(std::vector<Lit>) +
               element_sum;
    }

    std::size_t p_lazy_color_next(std::size_t already_colored) const {
        if (already_colored < 1'000'000) {
            if (already_colored < 256)
                return 256;
            return already_colored;
        }
        return 1'000'000;
    }

    /**
     * Much more efficient than adding all vertices individually;
     * that indeed becomes the bottleneck in color_lazy
     * if one does not use the bulk method.
     */
    template <typename Iterator>
    void p_bulk_add_vertices(Iterator begin, Iterator end) {
        std::vector<Bitset> cl_allowing_lit =
            p_compute_classes_allowing_literal();
        for (; begin != end; ++begin) {
            Vertex v = *begin;
            auto bs = make_larger_bitset(cl_allowing_lit[v.first]);
            bs &= cl_allowing_lit[v.second];
            vertex_queue.push(v.first, v.second, bs.count(), std::move(bs));
        }
    }

    std::vector<Bitset> p_compute_classes_allowing_literal() {
        const std::size_t n_cols = color_classes.all().size();
        const std::size_t nc_lit = 2 * n_concrete;
        std::vector<Bitset> classes_allowing_lit(nc_lit, Bitset(n_cols, true));
        for (std::size_t i = 0; i < n_cols; ++i) {
            const SharedDBPropagator& prop = color_classes[i];
            for (Lit l : prop.get_trail()) {
                if (l < nc_lit) {
                    l = lit::negate(l);
                    classes_allowing_lit[l][i] = false;
                }
            }
        }
        return classes_allowing_lit;
    }

    void p_prop_to_buffer(const SharedDBPropagator& prop, Bitset& buffer) {
        buffer.reset();
        for (Lit l : prop.get_trail()) {
            buffer[l] = true;
        }
    }

    bool p_find_first_suitable_candidate_one_new(Lit lmin, Lit lmax, Lit lnew) {
        for (Index cls_index : candidates) {
            auto& cls = color_classes[cls_index];
            if (!cls.is_false(lnew) &&
                color_classes.push_vertex(cls_index, lmin, lmax))
            {
                coloring_order.emplace_back(lmin, lmax);
                return true;
            }
        }
        return false;
    }

    bool p_find_first_suitable_candidate_from_all(Lit lmin, Lit lmax) {
        for (Index i = 0, nc = color_classes.num_classes(); i < nc; ++i) {
            auto& cls = color_classes[i];
            if (!cls.is_false(lmin) && !cls.is_false(lmax)) {
                if (color_classes.push_vertex(i, lmin, lmax)) {
                    coloring_order.emplace_back(lmin, lmax);
                    return true;
                }
            }
        }
        return false;
    }

    bool p_find_first_suitable_candidate(Lit lmin, Lit lmax) {
        for (Index cls_index : candidates) {
            auto& cls = color_classes[cls_index];
            if (!cls.is_false(lmin) && !cls.is_false(lmax)) {
                if (color_classes.push_vertex(cls_index, lmin, lmax)) {
                    return true;
                }
            }
        }
        return false;
    }

    ClauseDB* all_clauses;
    Lit n_concrete;
    ColorClasses<ColoringHeuristicSolver> color_classes;
    std::optional<PairInfeasibilityMap> internal_inf_map;
    PairInfeasibilityMap* inf_map;
    IndexedVertexHeap<MatrixIndexMap> vertex_queue;
    ClassCompleter<ColoringHeuristicSolver> class_completer;
    std::vector<std::vector<Lit>> explicit_partners_of;
    std::vector<Index> candidates;
    std::vector<Bitset> literals_in_class;
    Bitset m_old_true, m_new_true, m_buffer;
    mutable BitsetOperationsBuffer m_bitset_buffer;
    std::vector<Vertex> coloring_order;
    bool extracted_feasibilities = false;
    bool quiet = false;
};

} // namespace sammy

#endif
