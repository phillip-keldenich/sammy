#ifndef SAMMY_CLASS_COMPLETER_H_INCLUDED_
#define SAMMY_CLASS_COMPLETER_H_INCLUDED_

#include "literals.h"
#include "pair_infeasibility_map.h"

namespace sammy {

/**
 * A class that is used to complete 'color classes', i.e.,
 * turn partial configurations into complete valid configurations,
 * to finalize our partial assignments.
 * May fail to complete classes if the partial assignment turns out
 * to be invalid; in that case, we may have to cover vertices
 * that become uncovered during conflict resolution.
 */
template <typename EventHandler> class ClassCompleter {
  public:
    ClassCompleter(Lit n_concrete, Lit n_all, EventHandler* handler)
        : m_handler(handler), n_concrete(n_concrete), n_all(n_all),
          m_phases(n_all, false), m_given_literals_bits(2 * n_concrete, false),
          m_failed_literals_bits(2 * n_concrete, false), m_failed_qpos(0) {}

    /**
     * Try to turn the partial assignment represented by
     * the given propatator into a complete, valid configuration.
     * @return true if the assignment could be completed without
     *              changing any of the already-assigned literals.
     */
    bool complete_class(SharedDBPropagator& prop) {
        if (prop.get_trail().size() == n_all)
            return true;
        p_complete_init(prop);
        bool any_open = true;
        while (any_open) {
            if (p_push_from_failed_queue(prop)) {
                continue;
            }
            any_open = p_scan_all(prop);
        }
        for (Lit l : m_failed_concretes) {
            if (prop.is_false(l)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Compute the complete memory usage, in bytes, that this class completer
     * uses.
     */
    std::size_t total_memory_usage() const noexcept {
        return m_phases.bytes_used() + m_given_literals_bits.bytes_used() +
               m_failed_literals_bits.bytes_used() +
               m_given_concretes.capacity() * sizeof(Lit) +
               m_failed_concretes.capacity() * sizeof(Lit);
    }

  private:
    bool p_push_from_failed_queue(SharedDBPropagator& prop) {
        while (m_failed_qpos < m_failed_concretes.size()) {
            Lit l = m_failed_concretes[m_failed_qpos++];
            if (prop.is_open(l)) {
                p_push(prop, l);
                return true;
            }
        }
        return false;
    }

    bool p_scan_all(SharedDBPropagator& prop) {
        bool any_open = false;
        for (Lit v = 0; v < n_all; ++v) {
            Lit l = m_phases[v] ? lit::positive_lit(v) : lit::negative_lit(v);
            if (prop.is_open(l)) {
                any_open = true;
                if (!p_push(prop, l))
                    break;
            }
        }
        return any_open;
    }

    bool p_push(SharedDBPropagator& prop, Lit l) {
        struct AssignmentHandler {
            void assignment_undone(Lit l) const {
                if (lit::var(l) < that->n_concrete &&
                    that->m_given_literals_bits[l])
                {
                    if (!that->m_failed_literals_bits[l]) {
                        that->m_failed_literals_bits[l] = true;
                        that->m_failed_concretes.push_back(l);
                    }
                }
            }
            void assignment_forced(Lit) const {}

            ClassCompleter* that;
        };
        if (prop.push_level(l))
            return true;
        m_phases[lit::var(l)] = lit::negative(l);
        AssignmentHandler handler{this};
        if (!prop.resolve_conflicts(handler)) {
            throw UNSATError();
        }
        m_failed_qpos = 0;
        return false;
    }

    void p_complete_init(SharedDBPropagator& prop) {
        m_given_concretes = prop.get_trail();
        auto new_end =
            std::remove_if(m_given_concretes.begin(), m_given_concretes.end(),
                           [&](Lit l) { return lit::var(l) >= n_concrete; });
        m_given_concretes.erase(new_end, m_given_concretes.end());
        p_prop_to_buffer(prop, m_given_literals_bits);
        m_phases.reset();
        m_failed_concretes.clear();
        m_failed_qpos = 0;
        m_failed_literals_bits.reset();
    }

    void p_prop_to_buffer(const SharedDBPropagator& prop, Bitset& buffer) {
        buffer.reset();
        for (Lit l : prop.get_trail()) {
            if (lit::var(l) < n_concrete) {
                buffer[l] = true;
            }
        }
    }

    EventHandler* m_handler;
    Lit n_concrete, n_all;
    Bitset m_phases;
    Bitset m_given_literals_bits;
    Bitset m_failed_literals_bits;
    std::vector<Lit> m_given_concretes;
    std::vector<Lit> m_failed_concretes;
    std::size_t m_failed_qpos;
};

} // namespace sammy

#endif
