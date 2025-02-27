#ifndef SAMMY_BOUNDED_VARIABLE_ELIMINATION_H_INCLUDED_
#define SAMMY_BOUNDED_VARIABLE_ELIMINATION_H_INCLUDED_

#include "literals.h"
#include "simplify_datastructure.h"

namespace sammy {

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

class BoundedVariableEliminator {
  public:
    explicit BoundedVariableEliminator(SimplifyDatastructure* simplify_ds)
        : simplify_ds(simplify_ds),
          clauses_with_literal(2 * simplify_ds->original_num_vars()),
          stats(&stats_buffer) {
        p_init_cl_with_lit();
    }

    void set_stats(SimplificationStats* stats) noexcept { this->stats = stats; }

    bool is_clause_eliminated(CRef c) {
        std::size_t nc = simplify_ds->clauses().size();
        std::size_t nec = clause_eliminated.size();
        if (nec < nc) {
            clause_eliminated.reserve(std::max(2 * nec, nc + 32));
            clause_eliminated.resize(nc, false);
        }
        return clause_eliminated[c];
    }

    void set_clause_eliminated(CRef c) {
        std::size_t nc = simplify_ds->clauses().size();
        std::size_t nec = clause_eliminated.size();
        if (nec < nc) {
            clause_eliminated.reserve(std::max(2 * nec, nc + 32));
            clause_eliminated.resize(nc, false);
        }
        clause_eliminated[c] = true;
    }

    bool run_elimination(std::size_t max_gap) {
        bool changed = false, once_changed = false;
        do {
            changed = false;
            std::vector<std::pair<Var, std::size_t>> best_scores;
            changed = p_score_and_eliminate_pure(best_scores);
            if (!changed) {
                p_compress_scored(best_scores);
                bool eliminated_pure = false;
                for (const auto& e : best_scores) {
                    if (e.second == 0) {
                        changed = eliminated_pure = true;
                        p_eliminate_pure(e.first);
                        continue;
                    }
                    if (eliminated_pure)
                        break;
                    auto [score, baseline] = p_actual_score(e.first);
                    if (score <= baseline + max_gap) {
                        p_eliminate_variable(e.first);
                        changed = true;
                        break;
                    }
                }
            }
            once_changed |= changed;
        } while (changed);
        return once_changed;
    }

  private:
    bool p_score_and_eliminate_pure(
        std::vector<std::pair<Var, std::size_t>>& best_scores) {
        bool changed = false;
        const Var nv = simplify_ds->original_num_vars();
        for (Var v = 0; v < nv; ++v) {
            if (simplify_ds->is_eliminated(v) || simplify_ds->is_concrete(v))
                continue;
            Lit p = lit::positive_lit(v), n = lit::negative_lit(v);
            const auto& plist = clauses_with_literal[p];
            const auto& nlist = clauses_with_literal[n];
            if (plist.empty() || nlist.empty()) {
                p_eliminate_pure(v);
                changed = true;
                continue;
            }
            if (!changed)
                p_add_score(best_scores, v, plist.size() * nlist.size(), 20);
        }
        return changed;
    }

    void p_eliminate_variable(Var v) {
        Lit p = lit::positive_lit(v), n = lit::negative_lit(v);
        const auto& plist = clauses_with_literal[p];
        const auto& nlist = clauses_with_literal[n];
        auto& clauses = simplify_ds->clauses();
        for (CRef c1 : plist) {
            p_push_reconstruct(clauses[c1], p);
            set_clause_eliminated(c1);
        }
        for (CRef c2 : nlist) {
            p_push_reconstruct(clauses[c2], n);
            set_clause_eliminated(c2);
        }
        auto not_p = [=](Lit l) { return l != p; };
        auto not_n = [=](Lit l) { return l != n; };
        for (CRef c1 : plist) {
            for (CRef c2 : nlist) {
                SCVec buffer;
                std::copy_if(clauses[c1].begin(), clauses[c1].end(),
                             std::back_inserter(buffer), not_p);
                std::copy_if(clauses[c2].begin(), clauses[c2].end(),
                             std::back_inserter(buffer), not_n);
                if (buffer.empty())
                    throw UNSATError();
                std::sort(buffer.begin(), buffer.end());
                buffer.erase(std::unique(buffer.begin(), buffer.end()),
                             buffer.end());
                if (find_pair_if(buffer.begin(), buffer.end(),
                                 [](Lit pr, Lit c) {
                                     return pr == lit::negate(c);
                                 }) != buffer.end())
                {
                    continue;
                }
                CRef nclause = clauses.size();
                for (Lit l : buffer) {
                    clauses_with_literal[l].push_back(nclause);
                }
                clauses.emplace_back(std::move(buffer));
            }
        }
        simplify_ds->mark_eliminated(v);
        stats->variables_eliminated_by_resolution += 1;
    }

    std::pair<std::size_t, std::size_t> p_actual_score(Var v) {
        Lit p = lit::positive_lit(v), n = lit::negative_lit(v);
        std::size_t actual_score = 0;
        const auto& plist = clauses_with_literal[p];
        const auto& nlist = clauses_with_literal[n];
        const auto& clauses = simplify_ds->clauses();
        for (CRef c1 : plist) {
            for (CRef c2 : nlist) {
                actual_score +=
                    (p_resolvent_is_tautology(clauses[c1], clauses[c2], v) ? 0
                                                                           : 1);
            }
        }
        return {actual_score, plist.size() + nlist.size()};
    }

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

    void p_compress_scored(std::vector<std::pair<Var, std::size_t>>& scores) {
        for (auto& e : scores) {
            e.second = p_compress_lists(e.first);
        }
        auto compare = [](const std::pair<Var, std::size_t>& v1,
                          const std::pair<Var, std::size_t>& v2) {
            return v1.second < v2.second;
        };
        std::sort(scores.begin(), scores.end(), compare);
    }

    void p_add_score(std::vector<std::pair<Var, std::size_t>>& scores, Var v,
                     std::size_t score, std::size_t size_cap) {
        if (scores.size() < size_cap) {
            scores.emplace_back(v, score);
            return;
        }
        if (scores.back().second <= score)
            return;
        std::pair<Var, std::size_t> val{v, score};
        auto compare = [](const std::pair<Var, std::size_t>& v1,
                          const std::pair<Var, std::size_t>& v2) {
            return v1.second < v2.second;
        };
        auto insert_pos =
            std::lower_bound(scores.begin(), scores.end(), val, compare);
        std::move_backward(insert_pos, scores.end() - 1, scores.end());
        *insert_pos = std::pair<Var, std::size_t>{v, score};
    }

    void p_push_reconstruct(const SCVec& clause, Lit nelit) {
        SCVec buffer(clause);
        std::swap(*std::find(buffer.begin(), buffer.end(), nelit),
                  buffer.front());
        simplify_ds->push_reconstruction_clause(std::move(buffer));
    }

    void p_eliminate_pure(Var v) {
        Lit p = lit::positive_lit(v), n = lit::negative_lit(v);
        const auto& plist = clauses_with_literal[p];
        const auto& nlist = clauses_with_literal[n];
        const auto& nelist = (plist.empty() ? nlist : plist);
        Lit nelit = (plist.empty() ? n : p);
        simplify_ds->push_reconstruction_clause(SCVec(1, nelit));
        simplify_ds->mark_eliminated(v);
        stats->pure_literals_eliminated += 1;
        for (CRef c : nelist) {
            set_clause_eliminated(c);
        }
    }

    void p_init_cl_with_lit() {
        const auto& clauses = simplify_ds->clauses();
        for (CRef ci = 0, cn = clauses.size(); ci < cn; ++ci) {
            const auto& clause = clauses[ci];
            for (Lit l : clause) {
                clauses_with_literal[l].push_back(ci);
            }
        }
    }

    std::size_t p_compress_lists(Var v) {
        Lit p = lit::positive_lit(v);
        Lit n = lit::negative_lit(v);
        p_compress_list(p);
        p_compress_list(n);
        return clauses_with_literal[p].size() * clauses_with_literal[n].size();
    }

    void p_compress_list(Lit l) {
        auto is_eliminated = [&](CRef c) { return is_clause_eliminated(c); };
        auto& list = clauses_with_literal[l];
        list.erase(std::remove_if(list.begin(), list.end(), is_eliminated),
                   list.end());
    }

    SimplifyDatastructure* simplify_ds;
    std::vector<bool> clause_eliminated;
    std::vector<std::vector<CRef>> clauses_with_literal;
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
