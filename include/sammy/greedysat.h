#ifndef SAMMY_GREEDYSAT_H_INCLUDED_
#define SAMMY_GREEDYSAT_H_INCLUDED_

#include "literals.h"
#include "shared_db_propagator.h"

namespace sammy {

/**
 * A very simple greedy SAT solver based on SharedDBPropagator.
 * Whenever a concrete feature is unassigned, it picks an
 * arbitrary unassigned literal l and assigns it to the preferred
 * polarity.
 * The preferred polarity is given by PreferredAssignmentFn::operator[](l).
 * When all concrete features are assigned, any open non-concrete features
 * (of which there usually should be none) are greedily set to false
 * one after another, with propagation in between.
 */
template <typename PreferredAssignmentFn> class GreedySAT {
  public:
    GreedySAT(SharedDBPropagator* prop, Lit n_concrete,
              PreferredAssignmentFn&& assignment)
        : prop(prop),
          preferred_assignment(std::forward<PreferredAssignmentFn>(assignment)),
          n_concrete(n_concrete) {}

    bool solve() {
        while (prop->get_trail().size() < prop->db().num_vars()) {
            bool changed = false;
            for (Lit v = 0; v < n_concrete; ++v) {
                Lit l = preferred_assignment[v] ? lit::positive_lit(v)
                                                : lit::negative_lit(v);
                if (prop->is_open(l)) {
                    changed = true;
                    if (!prop->push_level(l) && !prop->resolve_conflicts()) {
                        throw UNSATError();
                    }
                }
            }
            if (changed)
                continue;
            Lit n_all = prop->db().num_vars();
            for (Lit v = n_concrete; v < n_all; ++v) {
                Lit l = lit::negative_lit(v);
                if (prop->is_open(l)) {
                    if (!prop->push_level(l)) {
                        prop->resolve_or_throw();
                    }
                }
            }
        }
        return true;
    }

  private:
    SharedDBPropagator* prop;
    PreferredAssignmentFn preferred_assignment;
    Lit n_concrete;
};

template <bool B> struct PreferValue {
    bool operator[](Lit /*var*/) const noexcept { return B; }
};
using PreferFalse = PreferValue<false>;
using PreferTrue = PreferValue<true>;

} // namespace sammy

#endif
