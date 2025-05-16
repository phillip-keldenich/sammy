#ifndef SAMMY_BOUNDED_VARIABLE_ELIMINATION_H_INCLUDED_
#define SAMMY_BOUNDED_VARIABLE_ELIMINATION_H_INCLUDED_

#include "literals.h"
#include "simplify_datastructure.h"

namespace sammy {

/**
 * @brief A simple class to track eliminated clauses.
 * 
 */
class EliminatedClausesTracker {
    public:
        EliminatedClausesTracker(SimplifyDatastructure* simplify_ds)
            : simplify_ds(simplify_ds),
             clauses_with_literal(2 * simplify_ds->original_num_vars()) // num literals = 2 * num vars
            {
                // initialize the clauses_with_literal vector
                const auto& clauses = simplify_ds->clauses();
                for (CRef ci = 0, cn = clauses.size(); ci < cn; ++ci) {
                    const auto& clause = clauses[ci];
                    for (Lit l : clause) {
                        clauses_with_literal[l].push_back(ci);
                    }
                }

            }

        void mark_eliminated(CRef c) {
            p_ensure_clause_eliminated_size();
            is_clause_eliminated[c] = true;
        }

        template <typename Range>
        void mark_eliminated_in_batch(Range&& range) {
            static_assert(std::is_same_v<decltype(std::begin(range)), decltype(std::end(range))>, 
                            "Range must be iterable with begin() and end() methods.");
            for (const CRef& c : range) {
                mark_eliminated(c);
            };
        }

        bool is_eliminated(CRef c) const {
            return c < is_clause_eliminated.size() && is_clause_eliminated[c];
        }

        void add_clause_with_literal(Lit l, CRef c) {
            clauses_with_literal[l].push_back(c);
        }

        std::vector<CRef>& get_clauses_with_literal(Lit l) {
            return clauses_with_literal[l];
        }

    /**
     * @brief Will calculate an upper bound on the number of resolvents.
     * The upper bound will be the product of the number of clauses with
     * the positive literal and the number of clauses with the negative
     * literal. This is a very rough estimate as usually we can eliminate
     * clauses, but it is a good starting point. If the result is zero,
     * the variable is guaranteed to be pure.
     * 
     * @param v The variable which we think about eliminating.
     * @param remove_eliminated As we are removing clauses, we may want to
     * remove eliminated clauses before calculating the number of resolvents.
     * @return std::size_t An upper bound on the number of resolvents.
     */
    std::size_t estimate_max_resolvents(Var v, bool remove_eliminated = true) {
        Lit p = lit::positive_lit(v);
        Lit n = lit::negative_lit(v);
        if (remove_eliminated) {
            p_remove_eliminated_clauses_for(p);
            p_remove_eliminated_clauses_for(n);
        }
        return clauses_with_literal[p].size() * clauses_with_literal[n].size();
    }

    private:
        void p_ensure_clause_eliminated_size() {
            std::size_t nc = simplify_ds->clauses().size();
            if (is_clause_eliminated.size() < nc) {
                is_clause_eliminated.reserve(std::max(2 * is_clause_eliminated.size(), nc + 32));
                is_clause_eliminated.resize(nc, false);
            }
        }

        void p_remove_eliminated_clauses_for(Lit l) {
            // Will remove eliminated clauses from the clause list for literal l.
            auto is_eliminated_ = [&](CRef c) { return is_eliminated(c); };
            auto& list = clauses_with_literal[l];
            list.erase(std::remove_if(list.begin(), list.end(), is_eliminated_),
                    list.end());
        }

        std::vector<bool> is_clause_eliminated;
        SimplifyDatastructure* simplify_ds; // Ensure this matches the actual class name
        std::vector<std::vector<CRef>> clauses_with_literal;
};


/**
 * Can we use BVE (bounded variable elimination) on non-concrete features?
 *  - Let (l1,l2) be a feasible interaction.
 *  - Assume we eliminate some non-concrete y by VE.
 *  - This adds the resolvents of all clauses with y and all clauses with -y.
 *  - As such, these clauses preserve logical equivalence: any solution to
 *    the old formula is a solution (including y) to the new formula.
 *  - Dropping y transforms this into a solution of the new formula, since
 *    removal of clauses maintains the solution in that direction.
 *  - Let (p1, p2) be feasible after the transformation.
 *  - Let Q be the model including (p1, p2) on the new formula.
 *  - Let Q+ be Q extended by y = True.
 *  - Analogously, let Q- be the model Q extended by y = False.
 *  - If Q+ is not a model of the original formula, there is a clause
 *    that is not satisfied; this must be one of the dropped clauses containing
 * -y.
 *  - Thus, all resolvents of this clause with clauses containing y must be
 *    satisfied in Q by satisfying the parts coming from clauses containing y.
 *  - Therefore, Q- must be a model of the original formula.
 * We should be allowed to use BVE on non-concrete features!
 * Do we introduce another dummy? Or some other way to represent elimination?
 */

using ScoredVar = std::pair<Var, std::size_t>;

class CappedScoreList: public std::vector<ScoredVar> {
  public:
    CappedScoreList(std::size_t max_size) : max_size(max_size) {}

    void add(ScoredVar e) {
        // Add the element into the list, expecting it to be sorted
        // in increasing order by the second element.
        // Insert the new element in the right place
        auto it = std::lower_bound(this->begin(), this->end(), e,
                                   [](const ScoredVar& a, const ScoredVar& b) {
                                       return a.second < b.second;
                                   });
        insert(it, e);
        // If the list is too long, remove the last element
        if (this->size() > max_size) {
            pop_back();
        }
    }

  private:
    std::size_t max_size;
};

class BoundedVariableEliminator {
  public:
    explicit BoundedVariableEliminator(SimplifyDatastructure* simplify_ds)
        : simplify_ds(simplify_ds),
            clause_tracker(simplify_ds),
          stats(&stats_buffer) {
    }

    void set_stats(SimplificationStats* stats) noexcept { this->stats = stats; }

    bool run_elimination(std::size_t max_gap) {
        bool changed = false, once_changed = false;
        do {
            changed = false;
            auto [best_scores, eliminated_pure] = p_score_and_eliminate_pure();
            once_changed |= eliminated_pure;
            for (const auto& e : best_scores) {
                assert(e.second > 0); // pure variables should have been eliminated
                Var var = e.first;
                if (p_elimination_cost(var) <= static_cast<int>(max_gap)) {
                    p_eliminate_variable(var);
                    changed = true;
                    break;
                }
            }
            once_changed |= changed;
        } while (changed);
        return once_changed;
    }

  private:
    // This function will provide a scored list of variables to eliminate.
    // If it finds pure variables during this, it will eliminate them directly.
    std::pair<CappedScoreList, bool> p_score_and_eliminate_pure() {
        // whether we eliminated a pure variable in any round
        bool eliminated_pure = false;
        // whether we eliminated a pure variable in the latest round
        bool eliminated_pure_in_round = false;
        CappedScoreList best_scores(20);
        do {
            // reset for the next round
            eliminated_pure_in_round = false;
            best_scores.clear();
            const auto num_vars = simplify_ds->original_num_vars();
            for (Var v = 0; v < num_vars; ++v) {
                if (simplify_ds->is_eliminated(v) || simplify_ds->is_concrete(v))
                    continue;
                auto num_resolvents = clause_tracker.estimate_max_resolvents(v);
                if (num_resolvents == 0) {
                    // if `num_resolvents` is zero, the variable is pure
                    p_eliminate_pure(v);
                    eliminated_pure_in_round = true;
                    eliminated_pure = true;
                    continue; // there may be some new pure variables in the already checked variables
                    // now, but it is faster to first continue and then do a second pass later.
                }
                // If we eliminated pure variables in this round, we will do a second pass
                if (!eliminated_pure_in_round){ best_scores.add(ScoredVar{v, num_resolvents}); }
            }
        } while(eliminated_pure_in_round);
        return {best_scores, eliminated_pure};
    }

    void p_eliminate_variable(Var v) {
        Lit p = lit::positive_lit(v), n = lit::negative_lit(v);
        const auto& plist = clause_tracker.get_clauses_with_literal(p);
        const auto& nlist = clause_tracker.get_clauses_with_literal(n);
        auto& clauses = simplify_ds->clauses();
        // Eliminate the clauses with p and n
        // and add them to the reconstruction stack
        for (CRef c1 : plist) {
            p_push_reconstruct(clauses[c1], p);
            clause_tracker.mark_eliminated(c1);
        }
        for (CRef c2 : nlist) {
            p_push_reconstruct(clauses[c2], n);
            clause_tracker.mark_eliminated(c2);
        }
        // Now we need to add the resolvents of all clauses with p and n
        auto not_p = [=](Lit l) { return l != p; };
        auto not_n = [=](Lit l) { return l != n; };
        for (CRef c1 : plist) {
            for (CRef c2 : nlist) {
                SCVec buffer;  // resolvent
                std::copy_if(clauses[c1].begin(), clauses[c1].end(),
                             std::back_inserter(buffer), not_p);
                std::copy_if(clauses[c2].begin(), clauses[c2].end(),
                             std::back_inserter(buffer), not_n);
                if (buffer.empty())
                    throw UNSATError();
                // the resolved is likely to have duplicates which we need to
                // remove
                std::sort(buffer.begin(), buffer.end());
                buffer.erase(std::unique(buffer.begin(), buffer.end()),
                             buffer.end());
                // frequently, the resolvent is a tautology, i.e., contain
                // a OR -a, and can be removed
                if (find_pair_if(buffer.begin(), buffer.end(),
                                 [](Lit pr, Lit c) {
                                     return pr == lit::negate(c);
                                 }) != buffer.end())
                {
                    continue;
                }
                CRef nclause = clauses.size();
                for (Lit l : buffer) {
                    clause_tracker.add_clause_with_literal(l, nclause);
                }
                clauses.emplace_back(std::move(buffer));
            }
        }
        // Finally, we can eliminate the variable
        simplify_ds->mark_eliminated(v);
        stats->variables_eliminated_by_resolution += 1;
    }

    // Calculate the cost of eliminating a variable. This is the number of
    // resolvents that will be added to the clause database minus the number
    // of clauses that will be removed. On badly chosen variables, the number
    // of resolvents can be very large, so we need to be careful as this optimization
    // will only have a positive impact if the number of resolvents is small.
    int p_elimination_cost(Var v) {
        Lit p = lit::positive_lit(v), n = lit::negative_lit(v);
        int actual_score = 0;
        const auto& plist = clause_tracker.get_clauses_with_literal(p);
        const auto& nlist = clause_tracker.get_clauses_with_literal(n);
        const auto& clauses = simplify_ds->clauses();
        for (CRef c1 : plist) {
            for (CRef c2 : nlist) {
                if(!p_resolvent_is_tautology(clauses[c1], clauses[c2], v)) {
                    actual_score += 1;
                }
            }
        }
        return actual_score - (plist.size() + nlist.size());
    }

    // A quick check to see if the resolvent will be a tautology
    bool p_resolvent_is_tautology(const SCVec& cpos, const SCVec& cneg,
                                  Var v) const {
        Lit p = lit::positive_lit(v), n = lit::negative_lit(v);
        for (Lit lp : cpos) {
            if (lp == p)
                continue;
            for (Lit ln : cneg) {
                if (ln == n)
                    continue;
                if (lp == lit::negate(ln)) {
                    return true;
                }
            }
        }
        return false;
    }

    void p_push_reconstruct(const SCVec& clause, Lit nelit) {
        SCVec buffer(clause);
        std::swap(*std::find(buffer.begin(), buffer.end(), nelit),
                  buffer.front());
        simplify_ds->push_reconstruction_clause(std::move(buffer));
    }

    void p_eliminate_pure(Var v) {
        Lit p = lit::positive_lit(v), n = lit::negative_lit(v);
        const auto& plist = clause_tracker.get_clauses_with_literal(p);
        const auto& nlist = clause_tracker.get_clauses_with_literal(n);
        const auto& nelist = (plist.empty() ? nlist : plist);
        Lit nelit = (plist.empty() ? n : p);
        simplify_ds->push_reconstruction_clause(SCVec(1, nelit));
        simplify_ds->mark_eliminated(v);
        stats->pure_literals_eliminated += 1;
        clause_tracker.mark_eliminated_in_batch(nelist);
    }

    SimplifyDatastructure* simplify_ds;
    EliminatedClausesTracker clause_tracker;
    SimplificationStats stats_buffer;
    SimplificationStats* stats;
};

inline bool bounded_variable_elimination(SimplifyDatastructure& simplifier,
                                         std::size_t max_gain_clauses,
                                         SimplificationStats* stats = nullptr) {
    BoundedVariableEliminator eliminator{&simplifier};
    if (stats)
        eliminator.set_stats(stats);
    bool result = eliminator.run_elimination(max_gain_clauses);
    if (result)
        simplifier.remove_eliminated_clauses();
    return result;
}

} // namespace sammy

#endif
