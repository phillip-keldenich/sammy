#ifndef SAMMY_SIMPLIFY_DATASTRUCTURE_H_INCLUDED_
#define SAMMY_SIMPLIFY_DATASTRUCTURE_H_INCLUDED_

#include "algorithm_ex.h"
#include "clause_db.h"
#include "eliminate_subsumed.h"
#include "literals.h"
#include <cassert>

namespace sammy {

struct SimplifiedInstance {
    std::vector<Var> new_to_old;
    ClauseDB formula;
    Var num_concrete;
};

class SimplifyDatastructure {
  public:
    explicit SimplifyDatastructure(const ClauseDB& formula, Var num_concrete)
        : m_variable_eliminated(formula.num_vars(), false),
          m_variable_concrete(formula.num_vars(), false),
          m_var_presence(formula.num_vars()),
          m_original_concrete(num_concrete) {
        std::fill_n(m_variable_concrete.begin(), num_concrete, true);
        for (Lit u : formula.unary_literals()) {
            m_clauses.emplace_back(1, u);
        }
        for (auto [l1, l2] : formula.binary_clauses()) {
            m_clauses.emplace_back(std::initializer_list<Lit>{l1, l2});
        }
        for (CRef c = 1, n = formula.literal_db_size(); c < n;
             c = formula.next_clause(c))
        {
            auto lits = formula.lits_of(c);
            m_clauses.emplace_back(lits.begin(), lits.end());
        }
        assert(!has_tautology_or_double_literal());
    }

    /**
     * A slowish (for debugging only) routine that checks if the given set
     * of literals is a (not necessarily strict) superset of any clause in the
     * formula.
     * @param set A set (must have .count(element) method) of literals.
     * @return true iff a clause of the formula subsumes the given set of
     * literals.
     */
    template <typename SetType> bool has_clause_subsuming(const SetType& set) {
        auto is_contained = [&](Lit l) -> bool { return set.count(l); };
        auto subsumes = [&](const SCVec& v) -> bool {
            return std::all_of(v.begin(), v.end(), is_contained);
        };
        return std::any_of(m_clauses.begin(), m_clauses.end(), subsumes);
    }

    /**
     * A debug method for checking that all clauses are non-tautological
     * and no clause contains several copies of the same literal.
     */
    bool has_tautology_or_double_literal() {
        for (const SCVec& cl : m_clauses) {
            m_var_presence.clear();
            for (Lit l : cl) {
                Var v = lit::var(l);
                if (!m_var_presence.check_insert(v)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * A debug method for finding duplicate binary clauses.
     */
    bool has_duplicate_binary_clause() {
        EdgeSet eset;
        for (const SCVec& cl : m_clauses) {
            if (cl.size() != 2)
                continue;
            std::pair<Lit, Lit> edge{std::min(cl[0], cl[1]),
                                     std::max(cl[0], cl[1])};
            if (!eset.insert(edge).second) {
                return true;
            }
        }
        return false;
    }

    bool is_eliminated(Var old_var) const noexcept {
        return m_variable_eliminated[old_var];
    }

    bool is_concrete(Var old_var) const noexcept {
        return m_variable_concrete[old_var];
    }

    void mark_eliminated(Var v) noexcept { m_variable_eliminated[v] = true; }

    void push_reconstruction_clause(const SCVec& c) {
        m_reconstruction_stack.emplace_back(c);
    }

    void push_reconstruction_clause(SCVec&& c) {
        m_reconstruction_stack.emplace_back(std::move(c));
    }

    const std::vector<SCVec>& get_reconstruction_stack() const noexcept {
        return m_reconstruction_stack;
    }

    void remove_eliminated_clauses() {
        auto clause_eliminated = [&](const SCVec& c) {
            return std::find_if(c.begin(), c.end(), [&](Lit l) {
                       return m_variable_eliminated[lit::var(l)];
                   }) != c.end();
        };
        m_clauses.erase(std::remove_if(m_clauses.begin(), m_clauses.end(),
                                       clause_eliminated),
                        m_clauses.end());
        eliminate_subsumed(m_clauses, original_num_vars());
    }

    bool apply_fixes_and_equalities(const std::vector<Lit>& old_to_new) {
        bool result = false;
        for (Var v = 0, nv = old_to_new.size() / 2; v < nv; ++v) {
            if (m_variable_eliminated[v])
                continue;
            Lit p = lit::positive_lit(v);
            Lit m = old_to_new[p];
            if (m == simplify::fixed_positive() ||
                m == simplify::fixed_negative())
            {
                p_apply_fix(v, m);
                result = true;
            } else if (m != p) {
                p_apply_remap(v, m);
                result = true;
            }
        }
        if (result) {
            p_fix_and_eq_rewrite_clauses(old_to_new);
        }
        return result;
    }

    ClauseDB extract_clause_db() {
        ClauseDB result{Var(m_variable_eliminated.size())};
        for (const SCVec& cl : m_clauses) {
            result.add_clause(cl.data(), cl.data() + cl.size());
        }
        return result;
    }

    void add_clauses(const ClauseDB& clause_db,
                     ClauseCounts old_clause_counts) {
        auto new_u = clause_db.unary_literals(
            old_clause_counts.unary_clause_end, clause_db.num_unaries());
        auto new_b = clause_db.binary_clauses(
            old_clause_counts.binary_clause_end, clause_db.num_binaries());
        for (Lit u : new_u) {
            m_clauses.emplace_back(1, u);
        }
        for (auto [l1, l2] : new_b) {
            m_clauses.emplace_back(std::initializer_list<Lit>{l1, l2});
        }
        for (CRef c = old_clause_counts.long_clause_end,
                  n = clause_db.literal_db_size();
             c < n; c = clause_db.next_clause(c))
        {
            auto lits = clause_db.lits_of(c);
            m_clauses.emplace_back(lits.begin(), lits.end());
        }
    }

    Var original_num_vars() const noexcept {
        return Var(m_variable_eliminated.size());
    }

    Var original_num_concrete() const noexcept {
        return m_original_concrete;
    }

    void sort_clauses() {
        auto sort_clause = [](SCVec& clause) {
            std::sort(clause.begin(), clause.end());
        };
        auto cmp_clauses = [](const SCVec& c1, const SCVec& c2) -> bool {
            return c1.size() < c2.size() ||
                   (c1.size() == c2.size() &&
                    std::lexicographical_compare(c1.begin(), c1.end(),
                                                 c2.begin(), c2.end()));
        };
        std::for_each(m_clauses.begin(), m_clauses.end(), sort_clause);
        std::sort(m_clauses.begin(), m_clauses.end(), cmp_clauses);
    }

    SimplifiedInstance compress() {
        Var old_nv = m_variable_eliminated.size();
        std::vector<Var> new_to_old;
        Var new_count = 0;
        std::vector<Lit> old_to_new(2 * old_nv, simplify::eliminated());
        for (Var v = 0; v < old_nv; ++v) {
            if (m_variable_eliminated[v])
                continue;
            if (m_variable_concrete[v]) {
                new_to_old.push_back(v);
                old_to_new[lit::positive_lit(v)] = lit::positive_lit(new_count);
                old_to_new[lit::negative_lit(v)] = lit::negative_lit(new_count);
                ++new_count;
            }
        }
        Var new_num_concrete = new_count;
        for (Var v = 0; v < old_nv; ++v) {
            if (m_variable_eliminated[v])
                continue;
            if (!m_variable_concrete[v]) {
                new_to_old.push_back(v);
                old_to_new[lit::positive_lit(v)] = lit::positive_lit(new_count);
                old_to_new[lit::negative_lit(v)] = lit::negative_lit(new_count);
                ++new_count;
            }
        }
        return SimplifiedInstance{std::move(new_to_old),
                                  p_compress_clauses(new_count, old_to_new),
                                  new_num_concrete};
    }

    std::vector<bool>
    reconstruct_solution(const SimplifiedInstance& simp,
                         const std::vector<bool>& simp_sol) const {
        Var new_num_var = simp_sol.size();
        Var old_num_var = m_variable_eliminated.size();
        std::vector<bool> result(old_num_var, false);
        for (Var v = 0; v < new_num_var; ++v) {
            result[simp.new_to_old[v]] = simp_sol[v];
        }
        auto is_satisfied = [&](Lit l) {
            return lit::negative(l) ^ result[lit::var(l)];
        };
        IteratorRange stack_order{m_reconstruction_stack.rbegin(),
                                  m_reconstruction_stack.rend()};
        for (const SCVec& rc : stack_order) {
            if (std::find_if(rc.begin(), rc.end(), is_satisfied) == rc.end()) {
                result[lit::var(rc[0])].flip();
            }
        }
        return result;
    }

    std::vector<Vertex>
    reconstruct_lb(const SimplifiedInstance& simp,
                   const std::vector<Vertex>& simp_lb) const {
        std::vector<Vertex> result;
        result.reserve(simp_lb.size());
        auto transform_literal = [&](Lit l) {
            Var ov = simp.new_to_old[lit::var(l)];
            return lit::negative(l) ? lit::negative_lit(ov)
                                    : lit::positive_lit(ov);
        };
        auto transform_vertex = [&](Vertex simp_v) {
            return Vertex{transform_literal(simp_v.first),
                          transform_literal(simp_v.second)};
        };
        std::transform(simp_lb.begin(), simp_lb.end(),
                       std::back_inserter(result), transform_vertex);
        const auto old_nclit = 2 * original_num_concrete();
        if(std::any_of(result.begin(), result.end(),
            [&](const Vertex& v) { return v.first >= old_nclit ||
                                          v.second >= old_nclit;}))
        {
            p_fix_nonconcrete_in_lb(result);
        }
        return result;
    }

    std::vector<std::vector<bool>>
    reconstruct_sample(const SimplifiedInstance& simplified_inst,
                       const std::vector<std::vector<bool>>& sample) const {
        std::vector<std::vector<bool>> result;
        result.reserve(sample.size());
        auto reconstruct = [&](const std::vector<bool>& cfg) {
            return reconstruct_solution(simplified_inst, cfg);
        };
        std::transform(sample.begin(), sample.end(), std::back_inserter(result),
                       reconstruct);
        return result;
    }

    std::vector<bool> reduce_solution(const SimplifiedInstance& simp,
                                      const std::vector<bool>& full_sol) const {
        Var nnew = simp.formula.num_vars();
        std::vector<bool> result(nnew, false);
        for (Var vnew = 0; vnew < nnew; ++vnew) {
            result[vnew] = full_sol[simp.new_to_old[vnew]];
        }
        return result;
    }

    std::vector<SCVec>& clauses() noexcept { return m_clauses; }

    const std::vector<SCVec>& clauses() const noexcept { return m_clauses; }

    void
    replace_all_binaries(const std::vector<std::pair<Lit, Lit>>& binaries) {
        auto is_long = [](const SCVec& cl) { return cl.size() > 2; };
        auto pair_to_scvec = [](const std::pair<Lit, Lit>& b) {
            return SCVec{std::initializer_list<Lit>{b.first, b.second}};
        };

        std::vector<SCVec> result;
        auto cnt_long =
            std::count_if(m_clauses.begin(), m_clauses.end(), is_long);
        result.reserve(binaries.size() + cnt_long);
        std::transform(binaries.begin(), binaries.end(),
                       std::back_inserter(result), pair_to_scvec);
        std::copy_if(m_clauses.begin(), m_clauses.end(),
                     std::back_inserter(result), is_long);
        m_clauses.swap(result);
    }

  private:
    void p_fix_nonconcrete_in_lb(std::vector<Vertex>& lb) const {
        Lit bound = 2 * original_num_concrete();
        auto process = [&] (Lit& l) {
            if(l < bound) {
                return;
            }
            Var v = lit::var(l);
            Lit lender = m_concreteness_lender.at(v);
            l = lit::negative(l) ? lit::negate(lender) : lender;
        };
        for(Vertex& v : lb) {
            process(v.first);
            process(v.second);
        }
    }

    ClauseDB p_compress_clauses(Var new_count,
                                const std::vector<Lit>& old_to_new) {
        ClauseDB result{new_count};
        CVec buffer;
        for (const SCVec& c : m_clauses) {
            buffer.clear();
            std::transform(c.begin(), c.end(), std::back_inserter(buffer),
                           [&](Lit l) { return old_to_new[l]; });
            assert(std::find(buffer.begin(), buffer.end(),
                             simplify::eliminated()) == buffer.end());
            result.add_clause(buffer.data(), buffer.data() + buffer.size());
        }
        return result;
    }

    void p_fix_and_eq_rewrite_clause(SCVec& c,
                                     const std::vector<Lit>& old_to_new) {
        auto transform_literal = [&](Lit l) { return old_to_new[l]; };
        auto not_fixed_negative = [&](Lit l) {
            return old_to_new[l] != simplify::fixed_negative();
        };
        auto fixed_positive = [&](Lit l) {
            return old_to_new[l] == simplify::fixed_positive();
        };
        if (std::any_of(c.begin(), c.end(), fixed_positive)) {
            c.clear();
            return;
        }
        c.erase(copy_transformed_if(c.begin(), c.end(), c.begin(),
                                    transform_literal, not_fixed_negative),
                c.end());
        if (c.empty())
            throw UNSATError();
        std::sort(c.begin(), c.end());
        c.erase(std::unique(c.begin(), c.end()), c.end());
        auto taut_pair = find_pair_if(c.begin(), c.end(), [](Lit l1, Lit l2) {
            return l1 == lit::negate(l2);
        });
        if (taut_pair != c.end()) {
            c.clear();
        }
    }

    void p_fix_and_eq_rewrite_clauses(const std::vector<Lit>& old_to_new) {
        std::for_each(m_clauses.begin(), m_clauses.end(), [&](SCVec& c) {
            p_fix_and_eq_rewrite_clause(c, old_to_new);
        });
        m_clauses.erase(
            std::remove_if(m_clauses.begin(), m_clauses.end(),
                           [](const SCVec& v) { return v.empty(); }),
            m_clauses.end());
    }

    void p_apply_fix(Var v, Lit l) {
        m_variable_eliminated[v] = true;
        Lit slit =
            (lit::negative(l) ? lit::negative_lit(v) : lit::positive_lit(v));
        m_reconstruction_stack.emplace_back(1, slit);
    }

    void p_apply_remap(Var v, Lit m) {
        if (m_variable_concrete[v]) {
            Var mvar = lit::var(m);
            if(!m_variable_concrete[mvar]) {
                m_variable_concrete[mvar] = true;
                bool flip_vm = lit::negative(m);
                Var lender_var = v;
                if (v >= m_original_concrete) {
                    Lit cl = m_concreteness_lender.at(v);
                    Var c = lit::var(cl);
                    if(c >= m_original_concrete) {
                        throw std::logic_error(
                            "Lender variable is not concrete");
                    }
                    bool flip_vc = lit::negative(cl);
                    // if v = c, then c = m and c is lender of m
                    // if v = -c, then c = -m and -c is lender of m
                    flip_vm ^= flip_vc;
                    lender_var = c;
                }
                m_concreteness_lender[mvar] = flip_vm ? 
                    lit::negative_lit(lender_var) :
                    lit::positive_lit(lender_var);
            }
        }
        Lit p = lit::positive_lit(v);
        Lit n = lit::negative_lit(v);
        Lit mneg = lit::negate(m);
        m_variable_eliminated[v] = true;
        m_reconstruction_stack.emplace_back(std::initializer_list<Lit>{n, m});
        m_reconstruction_stack.emplace_back(
            std::initializer_list<Lit>{p, mneg});
    }

    std::vector<SCVec> m_clauses;
    std::vector<SCVec> m_reconstruction_stack;
    std::vector<bool> m_variable_eliminated;
    std::vector<bool> m_variable_concrete;
    HashMap<Var, Lit> m_concreteness_lender;
    StampSet<Var, std::uint16_t> m_var_presence;
    Var m_original_concrete = 0;
};

} // namespace sammy

#endif
