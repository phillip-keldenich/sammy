#ifndef SAMMY_EQUALITY_GRAPH_H_INCLUDED_
#define SAMMY_EQUALITY_GRAPH_H_INCLUDED_

#include "literals.h"
#include "error.h"
#include <utility>
#include <algorithm>


namespace sammy {

/**
 * Basically a union-find datastructure
 * that allows making literals equal to other literals,
 * or to a true 'dummy' variable's literals
 * (i.e., fixing it to true or false).
 */
class EqualityGraph {
    std::vector<Var> m_component_root;
    std::vector<bool> m_component_root_negated;
    std::vector<std::uint8_t> m_rank;
    std::vector<std::pair<Var,bool>> m_path_buffer;
    Var m_dummy;  // we have a dummy variable 'True' at num_vars + 1

  public:
    explicit EqualityGraph(Var num_vars) :
        m_component_root(num_vars + 1, NIL),
        m_component_root_negated(num_vars + 1, false),
        m_rank(num_vars + 1, std::uint8_t(0)),
        m_dummy(num_vars)
    {
        m_path_buffer.reserve(64);
        // make sure the dummy will always be picked as representative
        m_rank[m_dummy] = std::numeric_limits<std::uint8_t>::max();
        for(Var v = 0; v <= num_vars; ++v) {
            m_component_root[v] = v;
        }
    }

    Var dummy() const noexcept {
        return m_dummy;
    }

    Lit find(Lit l) noexcept {
        Var v = lit::var(l);
        bool negated = lit::negative(l);
        Var r; bool rn;
        std::tie(r, rn) = p_find_pcompress(v);
        return (negated ^ rn) ? lit::negative_lit(r) : lit::positive_lit(r);
    }

    bool make_equal(Lit x, Lit y) {
        Var v1 = lit::var(x);
        Var v2 = lit::var(y);
        auto [r1, r1n] = p_find_pcompress(v1);
        auto [r2, r2n] = p_find_pcompress(v2);
        r1n ^= lit::negative(x);
        r2n ^= lit::negative(y);
        bool negated = r1n ^ r2n;
        if(r1 == r2) {
            if(negated) throw UNSATError();
            return false;
        }
        auto k1 = m_rank[v1];
        auto k2 = m_rank[v2];
        if(k1 < k2) {
            m_component_root[r1] = r2;
            m_component_root_negated[r1] = negated;
        } else {
            if(k1 == k2) {
                m_rank[r1] += 1;
            }
            m_component_root[r2] = r1;
            m_component_root_negated[r2] = negated;
        }
        return true;
    }

    bool make_true(Lit l) {
        return make_equal(l, lit::positive_lit(dummy()));
    }

    bool make_false(Lit l) {
        return make_equal(l, lit::negative_lit(dummy()));
    }

    std::vector<Lit> compute_old_to_new_map() {
        Var nv = m_component_root.size() - 1;
        std::vector<Lit> old_to_new(2 * nv, NIL);
        for(Var v = 0; v < nv; ++v) {
            Lit p = lit::positive_lit(v);
            Lit m = find(p);
            if(lit::var(m) == dummy()) {
                if(lit::negative(m)) {
                    old_to_new[p] = simplify::fixed_negative();
                    old_to_new[lit::negate(p)] = simplify::fixed_positive();
                } else {
                    old_to_new[p] = simplify::fixed_positive();
                    old_to_new[lit::negate(p)] = simplify::fixed_negative();
                }
            } else {
                old_to_new[p] = m;
                old_to_new[lit::negate(p)] = lit::negate(m);
            }
        }
        return old_to_new;
    }

  private:
    std::pair<Var,bool> p_find_pcompress(Var v) {
        Var p;
        bool pneg, rneg = false;
        for(;;) {
            p = m_component_root[v];
            pneg = m_component_root_negated[v];
            rneg ^= pneg;
            if(p == v) {
                bool cneg = false;
                for(auto i = m_path_buffer.rbegin(), e = m_path_buffer.rend(); i != e; ++i) {
                    auto [v_i, n_i] = *i;
                    cneg ^= n_i;
                    m_component_root[v_i] = p;
                    m_component_root_negated[v_i] = cneg;
                }
                m_path_buffer.clear();
                return {p, rneg};
            }
            m_path_buffer.emplace_back(v, pneg);
            v = p;
        }
    }
};

}

#endif
