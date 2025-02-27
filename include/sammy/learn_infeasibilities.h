#ifndef SAMMY_LEARN_INFEASIBILITIES_H_INCLUDED_
#define SAMMY_LEARN_INFEASIBILITIES_H_INCLUDED_

#include "literals.h"
#include "pair_infeasibility_map.h"
#include "shared_db_propagator.h"

namespace sammy {

namespace detail {

static inline void
learn_unary_infeasibilities(SharedDBPropagator& propagator,
                            const PairInfeasibilityMap* infeasibilities) {
    const auto& inf = *infeasibilities;
    const Lit nclit = 2 * inf.get_n_concrete();
    for (Lit l = 0; l < nclit; ++l) {
        if (inf[l][l]) {
            if (!propagator.is_false(l)) {
                Lit new_unary[1] = {lit::negate(l)};
                propagator.db().add_clause(+new_unary, new_unary + 1);
                propagator.incorporate_or_throw();
            }
        }
    }
}

} // namespace detail

/**
 * @brief Learn clauses so that all infeasible pairs
 *        are recognized easily by propagation, no matter
 *        which literal of the pair is pushed.
 *        Only learns clauses as they're needed, does
 *        not blindly create a new binary clause for each
 *        literal pair. Similarly handles missing unary clauses.
 *
 * @param clauses
 * @param infeasibilities
 */
inline void learn_infeasibilities(ClauseDB& clauses,
                                  const PairInfeasibilityMap* infeasibilities) {
    const Lit nclit = 2 * infeasibilities->get_n_concrete();
    const auto& inf = *infeasibilities;
    SharedDBPropagator propagator(&clauses);
    detail::learn_unary_infeasibilities(propagator, infeasibilities);
    for (Lit l = 0; l < nclit; ++l) {
        // the literal is false at level 0
        if (propagator.is_false(l)) {
            if (!inf[l][l]) {
                throw std::logic_error(
                    "Non-infeasible literal false at level 0!");
            }
            continue;
        }
        // the literal is true or open at level 0
        bool pushed_outer = false;
        if (propagator.is_open(l)) {
            if (!propagator.push_level(l)) {
                throw std::logic_error(
                    "Non-infeasible literal conflicts at level 1!");
            }
            pushed_outer = true;
        }
        // we are either at level 0 or 1, and l is now set to true
        for (Lit l2 : inf[l].ones()) {
            if (!propagator.is_false(l2)) {
                Lit new_binary[2] = {lit::negate(l), lit::negate(l2)};
                clauses.add_clause(+new_binary, new_binary + 2);
                if (pushed_outer)
                    propagator.pop_level();
                propagator.incorporate_or_throw();
                pushed_outer = false;
                if (propagator.is_open(l)) {
                    if (!propagator.push_level(l)) {
                        throw std::logic_error(
                            "Non-infeasible literal conflicts at level 1!");
                    }
                    pushed_outer = true;
                } else if (propagator.is_false(l)) {
                    throw std::logic_error(
                        "Non-infeasible literal false at level 0!");
                }
                assert(propagator.is_false(l2));
            }
        }
        if (pushed_outer)
            propagator.pop_level();
        assert(propagator.get_current_level() == 0);
    }
}

} // namespace sammy

#endif
