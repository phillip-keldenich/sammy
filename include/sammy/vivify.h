#ifndef SAMMY_VIVIFY_H_INCLUDED_
#define SAMMY_VIVIFY_H_INCLUDED_

#include "shared_db_propagator.h"
#include "simplification_stats.h"
#include "simplify_datastructure.h"
#include <cassert>

namespace sammy {

class ClauseVivifier {
  public:
    explicit ClauseVivifier(SimplifyDatastructure* simplify_ds)
        : simplify_ds(simplify_ds), clauses((simplify_ds->sort_clauses(),
                                             simplify_ds->extract_clause_db())),
          propagator(&clauses), first_long_index(p_find_first_long()),
          stats(&stats_buffer) {}

    void set_stats(SimplificationStats* stats) { this->stats = stats; }

    bool run() {
        if (!propagator.get_trail().empty())
            return true;
        CRef c = 1, n = clauses.literal_db_size();
        CRef coffset = 0;
        bool changed = false;
        abort_run = false;
        for (; c < n; ++coffset, c = clauses.next_clause(c)) {
            auto size = clauses.clause_length(c);
            if (size > 3) {
                break;
            }
            changed |= p_vivify_ternary(clauses.lits_of(c), coffset);
            if (abort_run)
                return changed;
        }
        for (; c < n; ++coffset, c = clauses.next_clause(c)) {
            changed |= p_vivify_longer(clauses.lits_of(c), coffset);
            if (abort_run)
                return changed;
        }
        return changed;
    }

  private:
    CRef p_find_first_long() {
        auto& sc = simplify_ds->clauses();
        SCVec dummy(3, NIL);
        auto first_long_iter = std::lower_bound(
            sc.begin(), sc.end(), dummy, [](const SCVec& v1, const SCVec& v2) {
                return v1.size() < v2.size();
            });
        return CRef(first_long_iter - sc.begin());
    }

    bool p_vivify_longer(ClauseDB::Lits lits_, CRef coffset) {
        const Lit* lits = lits_.begin();
        auto [first_diff_pushed, first_diff_lits] =
            std::mismatch(pushed.begin(), pushed.end(), lits);
        if (first_diff_pushed == pushed.end()) {
            // the latter clause is subsumed by the first;
            // let subsumption detection handle this
            return false;
        }
        p_vivify_longer_pop_mismatched(first_diff_pushed);
        if (abort_run) {
            return true;
        }
        auto& scclause = simplify_ds->clauses()[first_long_index + coffset];
        for (auto pos = first_diff_lits; pos != lits_.end() - 1; ++pos) {
            Lit lcur = *pos;
            if (!propagator.is_open(lcur)) {
                if (propagator.is_true(lcur)) {
                    // (-l1 & ... & -lk -> lcur) => shorten clause
                    scclause.assign(lits_.begin(), pos + 1);
                } else {
                    // (-l1 & ... & -lk -> -lcur) === (l1 v ... v lk v -lcur)
                    // => remove lcur by resolution
                    scclause.assign(lits_.begin(), pos);
                    scclause.insert(scclause.end(), pos + 1, lits_.end());
                    clauses.add_clause(scclause.begin(), scclause.end());
                }
                stats->clauses_strengthened_by_vivification += 1;
                return true;
            }
            if (!propagator.push_level(lit::negate(lcur))) {
                // (-l1 & ... & -lk & -lcur) -> conflict => shorten clause
                scclause.assign(lits_.begin(), pos + 1);
                clauses.add_clause(scclause.begin(), scclause.end());
                stats->clauses_strengthened_by_vivification += 1;
                propagator.pop_level();
                return true;
            }
            pushed.push_back(lcur);
        }
        return false;
    }

    void p_vivify_longer_pop_mismatched(
        std::vector<Lit>::const_iterator first_diff_pushed) {
        // pop the additional levels from the propagator
        for (auto i = first_diff_pushed; i != pushed.end(); ++i) {
            propagator.pop_level();
        }
        pushed.erase(first_diff_pushed, pushed.end());
        if (pushed.empty()) {
            propagator.incorporate_or_throw();
            if (!propagator.get_trail().empty()) {
                abort_run = true;
            }
        }
    }

    bool p_vivify_ternary(ClauseDB::Lits lits_, CRef coffset) {
        assert(propagator.get_current_level() == 0);
        const Lit* lits = lits_.begin();
        return p_vivify_ternary_pair_conflicts(lits[0], lits[1], lits[2],
                                               coffset) ||
               p_vivify_ternary_pair_conflicts(lits[0], lits[2], lits[1],
                                               coffset) ||
               p_vivify_ternary_pair_conflicts(lits[1], lits[2], lits[0],
                                               coffset);
    }

    bool p_vivify_ternary_pair_conflicts(Lit l1, Lit l2, Lit l3, CRef coffset) {
        if (!propagator.push_level(lit::negate(l1))) {
            stats->failed_literals += 1;
            propagator.pop_level();
            simplify_ds->clauses().emplace_back(1, l1);
            abort_run = true;
            return true;
        }
        Lit replace_clause[2] = {l1, l2};
        bool replace = false;
        if (!propagator.is_open(l2)) {
            // -l1 -> l2 or -l1 -> -l2
            // -l1 -> l2 === l1 v l2 => subsumes clause.
            // -l1 -> -l2 === l1 v -l2 => resolution with clause on l2 gives l1
            // v l3
            replace = true;
            if (propagator.is_false(l2)) {
                replace_clause[1] = l3;
            }
        }
        if (!replace) {
            if (!propagator.push_level(lit::negate(l2))) {
                // -l1 & -l2 -> conflict
                // learn l1 v l2 => subsumes clause
                replace = true;
            }
            propagator.pop_level();
        }
        if (replace) {
            clauses.add_clause(+replace_clause, replace_clause + 2);
            simplify_ds->clauses()[first_long_index + coffset].assign(
                +replace_clause, replace_clause + 2);
            stats->clauses_strengthened_by_vivification += 1;
        }
        propagator.pop_level();
        if (replace) {
            propagator.incorporate_or_throw();
            if (!propagator.get_trail().empty()) {
                abort_run = true;
            }
        }
        return replace;
    }

    SimplifyDatastructure* simplify_ds;
    ClauseDB clauses;
    SharedDBPropagator propagator;
    std::vector<Lit> pushed;
    CRef first_long_index;
    SimplificationStats stats_buffer;
    SimplificationStats* stats;
    bool abort_run = false;
};

inline bool vivify(SimplifyDatastructure& simplifier,
                   SimplificationStats* stats = nullptr) {
    ClauseVivifier vivifier(&simplifier);
    if (stats)
        vivifier.set_stats(stats);
    return vivifier.run();
}

} // namespace sammy

#endif
