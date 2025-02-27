#ifndef SAMMY_PARTIAL_SOLUTION_H_INCLUDED_
#define SAMMY_PARTIAL_SOLUTION_H_INCLUDED_

#include "dynamic_bitset.h"
#include "literals.h"
#include "pair_infeasibility_map.h"
#include "thread_interrupt.h"

namespace sammy {

class PartialSolution {
  public:
    explicit PartialSolution(Var num_vars, const PairInfeasibilityMap* inf_map)
        : PartialSolution(num_vars, inf_map,
                          static_cast<const DynamicBitset*>(nullptr),
                          static_cast<const DynamicBitset*>(nullptr)) {}

    template <typename ClassIterator>
    explicit PartialSolution(Var num_vars, const PairInfeasibilityMap* inf_map,
                             ClassIterator classes_begin,
                             ClassIterator classes_end)
        : m_assignments_with_literal(2 * inf_map->get_n_concrete()),
          m_inf_map(inf_map), m_num_vars(num_vars),
          m_row_uncolored(2 * inf_map->get_n_concrete(), true) {
        using ClassType = std::remove_reference_t<std::remove_cv_t<
            typename std::iterator_traits<ClassIterator>::value_type>>;

        m_assignments.reserve(std::distance(classes_begin, classes_end));
        m_true_concrete_literals.reserve(
            std::distance(classes_begin, classes_end));
        if constexpr (std::is_same_v<ClassType, SharedDBPropagator>) {
            p_from_props(classes_begin, classes_end);
        } else {
            p_from_bitsets(classes_begin, classes_end);
        }
    }

    template <typename Callable /*(Lit lmin, Lit lmax)*/>
    void iterate_all_uncovered(Callable&& callable,
                               bool interruptible = false) {
        Lit nclit = 2 * m_inf_map->get_n_concrete();
        std::size_t count = 0;
        for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
            m_row_uncolored.set();
            m_row_uncolored ^= (*m_inf_map)[lmin];
            const auto& with_lit = m_assignments_with_literal[lmin];
            for (Index i : with_lit) {
                m_row_uncolored.binary_subtract(m_true_concrete_literals[i],
                                                lmin + 1);
            }
            for (Lit lmax : m_row_uncolored.ones_from(lmin + 1)) {
                (std::forward<Callable>(callable))(lmin, lmax);
            }
            if (interruptible && ++count == 512) {
                count = 0;
                throw_if_interrupted();
            }
        }
    }

    std::size_t get_n_concrete() const noexcept {
        return m_inf_map->get_n_concrete();
    }

    std::size_t count_uncovered() noexcept {
        std::size_t result = 0;
        Lit nclit = 2 * m_inf_map->get_n_concrete();
        for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
            m_row_uncolored.set();
            m_row_uncolored ^= (*m_inf_map)[lmin];
            const auto& with_lit = m_assignments_with_literal[lmin];
            for (Index i : with_lit) {
                m_row_uncolored.binary_subtract(m_true_concrete_literals[i],
                                                lmin + 1);
            }
            result += m_row_uncolored.count_from(lmin + 1);
        }
        return result;
    }

    template <typename Callable /*(Lit lmin, Lit lmax)*/>
    void iterate_all_uniquely_covered(Callable&& callable,
                                      bool cover_is_incomplete = false) {
        if (cover_is_incomplete) {
            p_iterate_uniquely_covered_incomplete(
                std::forward<Callable>(callable));
        } else {
            p_iterate_uniquely_covered_complete(
                std::forward<Callable>(callable));
        }
    }

    template <typename Callable /*(Index class_index, Lit lmin, Lit lmax)*/>
    void iterate_all_uniquely_covered_with_class(Callable&& callable) {
        Lit nclit = 2 * m_inf_map->get_n_concrete();
        // track everything that is multiply-covered
        DynamicBitset multi_covered = make_bitset(nclit, false);
        // track everything that is covered or infeasible
        DynamicBitset covered = make_bitset(nclit, false);
        // temporary buffer
        DynamicBitset tmp = make_bitset(nclit, false);
        for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
            const auto& with_lit = m_assignments_with_literal[lmin];
            covered = (*m_inf_map)[lmin];
            multi_covered = covered;
            for (Index cci : with_lit) {
                const auto& ccbits = m_true_concrete_literals[cci];
                tmp = covered;
                tmp &= ccbits; // tmp = set of previously covered that are by
                               // config cci
                multi_covered |= tmp;
                covered |= ccbits;
            }
            multi_covered.flip(); // multi_covered = everything that's not
                                  // covered multiple times
            multi_covered &=
                covered; // multi_covered = everything that's not covered
                         // multiple times and not uncovered
            for (Index cci : with_lit) {
                tmp = m_true_concrete_literals[cci];
                tmp &= multi_covered;
                for (Lit lmax : tmp.ones_from(lmin + 1)) {
                    (std::forward<Callable>(callable))(cci, lmin, lmax);
                }
            }
        }
    }

    Index add_class(const SharedDBPropagator& prop) {
        Index i = m_assignments.size();
        p_add_class_from_prop(i, prop);
        return i;
    }

    template <typename BitsetType>
    Index add_assignment(BitsetType&& assignment) {
        Index i = m_assignments.size();
        p_add_class_from_bitset(i, std::forward<BitsetType>(assignment));
        return i;
    }

    std::vector<DynamicBitset>
    remove_assignments(const std::vector<Index>& indices) {
        std::vector<DynamicBitset> result;
        result.reserve(indices.size());
        DynamicBitset is_deleted(m_assignments.size(), false);
        for (Index i : indices)
            is_deleted[i].set();
        std::vector<Index> old_to_new(m_assignments.size(), NIL);
        for (std::size_t i = 0, o = 0, s = m_assignments.size(); i < s; ++i) {
            if (!is_deleted[i]) {
                old_to_new[i] = o++;
            } else {
                result.emplace_back(std::move(m_assignments[i]));
                m_true_concrete_literals[i].clear();
            }
        }
        m_assignments.erase(std::remove_if(m_assignments.begin(),
                                           m_assignments.end(),
                                           [](const DynamicBitset& b) {
                                               return b.size() == 0;
                                           }),
                            m_assignments.end());
        m_true_concrete_literals.erase(
            std::remove_if(
                m_true_concrete_literals.begin(),
                m_true_concrete_literals.end(),
                [](const DynamicBitset& b) { return b.size() == 0; }),
            m_true_concrete_literals.end());
        std::for_each(
            m_assignments_with_literal.begin(),
            m_assignments_with_literal.end(), [&](std::vector<Index>& sinds) {
                std::transform(sinds.begin(), sinds.end(), sinds.begin(),
                               [&](Index x) { return old_to_new[x]; });
                sinds.erase(std::remove(sinds.begin(), sinds.end(), NIL),
                            sinds.end());
            });
        return result;
    }

    std::size_t size() const noexcept { return m_assignments.size(); }

    bool empty() const noexcept { return m_assignments.empty(); }

    const DynamicBitset&
    get_assignment(std::size_t class_index) const noexcept {
        return m_assignments[class_index];
    }

    const DynamicBitset&
    get_true_concrete_literals(std::size_t class_index) const noexcept {
        return m_true_concrete_literals[class_index];
    }

    const std::vector<DynamicBitset>& assignments() const noexcept {
        return m_assignments;
    }

    template <typename BitsetType>
    std::vector<BitsetType> assignments_as() const {
        std::vector<BitsetType> result;
        for (const DynamicBitset& assignment : m_assignments) {
            BitsetType a(m_num_vars, false);
            for (Var v : range(m_num_vars)) {
                if (assignment[v])
                    a[v] = true;
            }
            result.emplace_back(std::move(a));
        }
        return result;
    }

    Index find_covering_class(Lit lmin, Lit lmax) const {
        if ((*m_inf_map)(lmin, lmax))
            return NIL;
        const auto& l1 = m_assignments_with_literal[lmin];
        const auto& l2 = m_assignments_with_literal[lmax];
        const auto& list = l1.size() < l2.size() ? l1 : l2;
        const Lit find_lit = l1.size() < l2.size() ? lmax : lmin;
        for (Index cls_index : list) {
            if (m_true_concrete_literals[cls_index][find_lit])
                return cls_index;
        }
        throw std::out_of_range(
            "Vertex (" + std::to_string(lmin) + ", " + std::to_string(lmax) +
            ") is not covered and not marked as infeasible!");
    }

    template <typename Callable /*(Index)*/>
    void find_covering_classes(Lit lmin, Lit lmax, Callable&& callable) const {
        if ((*m_inf_map)(lmin, lmax))
            return;
        const auto& l1 = m_assignments_with_literal[lmin];
        const auto& l2 = m_assignments_with_literal[lmax];
        const auto& list = l1.size() < l2.size() ? l1 : l2;
        const Lit find_lit = l1.size() < l2.size() ? lmax : lmin;
        for (Index cls_index : list) {
            if (m_true_concrete_literals[cls_index][find_lit]) {
                std::forward<Callable>(callable)(cls_index);
            }
        }
    }

    std::size_t bytes_used() const noexcept {
        std::size_t result = sizeof(PartialSolution);
        result += m_assignments.capacity() * sizeof(DynamicBitset);
        if (!m_assignments.empty()) {
            result += m_assignments[0].bytes_used() * m_assignments.size();
        }
        result += m_true_concrete_literals.capacity() * sizeof(DynamicBitset);
        if (!m_true_concrete_literals.empty()) {
            result += m_true_concrete_literals[0].bytes_used() *
                      m_true_concrete_literals.size();
        }
        result +=
            sizeof(std::vector<Index>) * m_assignments_with_literal.capacity();
        for (const auto& v : m_assignments_with_literal) {
            result += sizeof(Index) * v.capacity();
        }
        return result * (CHAR_BIT / 8);
    }

  private:
    template <typename Callable>
    void p_iterate_uniquely_covered_complete(Callable&& callable) {
        Lit nclit = 2 * m_inf_map->get_n_concrete();
        DynamicBitset multi_covered = make_bitset(nclit, false);
        DynamicBitset prev_covered = make_bitset(nclit, false);
        DynamicBitset tmp = make_bitset(nclit, false);
        for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
            multi_covered.reset();
            prev_covered.reset();
            const auto& with_lit = m_assignments_with_literal[lmin];
            for (Index cci : with_lit) {
                const auto& ccbits = m_true_concrete_literals[cci];
                tmp = prev_covered;
                tmp &= ccbits;
                multi_covered |= tmp;
                prev_covered |= ccbits;
            }
            tmp = multi_covered;
            tmp.flip();
            tmp -= (*m_inf_map)[lmin];
            for (Lit lmax : tmp.ones_from(lmin + 1)) {
                (std::forward<Callable>(callable))(lmin, lmax);
            }
        }
    }

    template <typename Callable>
    void p_iterate_uniquely_covered_incomplete(Callable&& callable) {
        Lit nclit = 2 * m_inf_map->get_n_concrete();
        // track everything that is multiply-covered
        DynamicBitset multi_covered = make_bitset(nclit, false);
        // track everything that is covered or infeasible
        DynamicBitset covered = make_bitset(nclit, false);
        // temporary buffer
        DynamicBitset tmp = make_bitset(nclit, false);
        for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
            const auto& with_lit = m_assignments_with_literal[lmin];
            covered = (*m_inf_map)[lmin];
            multi_covered = covered;
            for (Index cci : with_lit) {
                const auto& ccbits = m_true_concrete_literals[cci];
                tmp = covered;
                tmp &= ccbits; // tmp = set of previously covered that are by
                               // config cci
                multi_covered |= tmp;
                covered |= ccbits;
            }
            multi_covered.flip(); // multi_covered = everything that's not
                                  // covered multiple times
            multi_covered &=
                covered; // multi_covered = everything that's not covered
                         // multiple times and not uncovered
            for (Lit lmax : multi_covered.ones_from(lmin + 1)) {
                (std::forward<Callable>(callable))(lmin, lmax);
            }
        }
    }

    template <typename Iterator>
    void p_from_props(Iterator begin, Iterator end) {
        for (Index iclass = 0; begin != end; ++begin, ++iclass) {
            const auto& propagator = *begin;
            p_add_class_from_prop(iclass, propagator);
        }
    }

    void p_add_class_from_prop(Index class_index,
                               const SharedDBPropagator& prop) {
        Var n_conc = m_inf_map->get_n_concrete();
        DynamicBitset true_concrete(2 * n_conc, false);
        DynamicBitset assignment(m_num_vars, false);
        for (Var v : range(n_conc)) {
            bool vtrue = prop.is_true(lit::positive_lit(v));
            Lit t = vtrue ? lit::positive_lit(v) : lit::negative_lit(v);
            true_concrete[t].set();
            m_assignments_with_literal[t].push_back(class_index);
            assignment[v] = vtrue;
        }
        for (Var v : range(n_conc, m_num_vars)) {
            bool vtrue = prop.is_true(lit::positive_lit(v));
            assignment[v] = vtrue;
        }
        m_assignments.emplace_back(std::move(assignment));
        m_true_concrete_literals.emplace_back(std::move(true_concrete));
    }

    template <typename BitsetType>
    void p_add_class_from_bitset(Index class_index, BitsetType&& assignment) {
        Var n_conc = m_inf_map->get_n_concrete();
        DynamicBitset true_concrete(2 * n_conc, false);
        for (Var v : range(n_conc)) {
            Lit t = assignment[v] ? lit::positive_lit(v) : lit::negative_lit(v);
            true_concrete[t].set();
            m_assignments_with_literal[t].push_back(class_index);
        }
        m_assignments.emplace_back(std::forward<BitsetType>(assignment));
        m_true_concrete_literals.emplace_back(std::move(true_concrete));
    }

    template <typename Iterator>
    void p_from_bitsets(Iterator begin, Iterator end) {
        std::size_t index = 0;
        for (const auto& c : IteratorRange{begin, end}) {
            p_add_class_from_bitset(index++, c);
        }
    }

    std::vector<DynamicBitset> m_assignments;
    std::vector<DynamicBitset> m_true_concrete_literals;
    std::vector<std::vector<Index>> m_assignments_with_literal;
    const PairInfeasibilityMap* m_inf_map;
    Var m_num_vars;
    DynamicBitset m_row_uncolored;
};

} // namespace sammy

#endif
