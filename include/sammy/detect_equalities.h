#ifndef SAMMY_DETECT_EQUALITIES_H_INCLUDED_
#define SAMMY_DETECT_EQUALITIES_H_INCLUDED_

#include "equality_graph.h"
#include "range.h"
#include "shared_db_propagator.h"
#include "simplification_stats.h"
#include "simplify_datastructure.h"
#include "stamp_set.h"
#include <cassert>

namespace sammy {

namespace detail {

static inline auto implied_literals(const SharedDBPropagator& prop) noexcept {
    return IteratorRange{prop.current_level_begin() + 1,
                         prop.get_trail().end()};
}

static inline bool fix_unaries(const SharedDBPropagator& propagator,
                               EqualityGraph& equality_graph,
                               SimplificationStats* stats) {
    bool last_changed = false;
    for (Lit l : propagator.get_trail()) {
        if (equality_graph.make_true(l)) {
            stats->variables_fixed += 1;
            last_changed = true;
        }
    }
    return last_changed;
}

} // namespace detail

class ImplicationGraphBuilder {
  public:
    explicit ImplicationGraphBuilder(SimplifyDatastructure* simplify_ds)
        : simplify_ds(simplify_ds), clauses(simplify_ds->extract_clause_db()),
          old_counts(clauses.get_clause_counts()), propagator(&clauses),
          equality_graph(simplify_ds->original_num_vars()),
          implications_of(2 * simplify_ds->original_num_vars()),
          additional_implications_of(2 * simplify_ds->original_num_vars()),
          stats(&stats_buffer),
          reachable(2 * simplify_ds->original_num_vars()) {}

    void set_stats(SimplificationStats* stats) noexcept { this->stats = stats; }

    bool run_fixing_and_equality_iteration() {
        bool changed = false;
        bool last_changed;
        do {
            last_changed = p_init_iteration();
            for (Var v = 0, nv = simplify_ds->original_num_vars(); v < nv; ++v)
            {
                if (simplify_ds->is_eliminated(v) ||
                    !propagator.is_open(lit::positive_lit(v)))
                {
                    continue;
                }
                last_changed |= p_push_var(v);
            }
            p_add_additional_arcs();
            last_changed |= p_tarjan_scc_find_equalities();
            changed |= last_changed;
        } while (last_changed);
        if (changed) {
            simplify_ds->add_clauses(clauses, old_counts);
            simplify_ds->apply_fixes_and_equalities(
                equality_graph.compute_old_to_new_map());
        }
        return changed;
    }

    bool run_contrapositive_and_strengthening_iteration() {
        return p_add_missing_contrapositivities() ||
               p_strengthen_with_implications();
    }

  private:
    bool p_strengthen_with_implications() {
        auto clauses_with_lit = p_make_clauses_with_lit();
        auto& clauses = simplify_ds->clauses();
        auto is_implied = [&](Lit l) { return reachable.contains(l); };
        bool result = false;
        for (Lit l = 0, nl = 2 * simplify_ds->original_num_vars(); l < nl; ++l)
        {
            if (simplify_ds->is_eliminated(lit::var(l)))
                continue;
            if (clauses_with_lit[l].empty())
                continue;
            p_stamp_reachable_black(l, true);
            reachable.erase(l);
            for (CRef cr : clauses_with_lit[l]) {
                auto& clause = clauses[cr];
                if (std::any_of(clause.begin(), clause.end(), is_implied)) {
                    stats->clauses_strengthened_by_binary_resolution += 1;
                    auto iter = std::find(clause.begin(), clause.end(), l);
                    std::swap(*iter, clause.back());
                    clause.pop_back();
                    if (clause.size() <= 2) {
                        result = true;
                    }
                }
            }
        }
        return result;
    }

    void p_stamp_reachable_black(Lit v, bool clear_first) {
        tj_stack.clear();
        if (clear_first)
            reachable.clear();
        reachable.insert(v);
        tj_stack.push_back(v);
        while (!tj_stack.empty()) {
            Lit w = tj_stack.back();
            tj_stack.pop_back();
            for (Lit x : implications_of[w]) {
                if (reachable.check_insert(x)) {
                    tj_stack.push_back(x);
                }
            }
        }
    }

    bool p_add_missing_contrapositivities() {
        Lit rlast = NIL;
        auto is_last_pushed = [&](Lit w) { return w == rlast; };
        std::size_t added_clauses = 0;
        auto& clauses = simplify_ds->clauses();
        for (Lit v : tj_order) {
            auto& aio = additional_implications_of[v];
            if (aio.empty())
                continue;
            auto& io = implications_of[v];
            p_stamp_reachable_black(
                v, !std::any_of(io.begin(), io.end(), is_last_pushed));
            rlast = v;
            for (Lit w : aio) {
                if (!reachable.contains(w)) {
                    if (arc_set.emplace(v, w).second) {
                        Lit lits[2] = {lit::negate(v), w};
                        io.push_back(w);
                        clauses.emplace_back(+lits, lits + 2);
                        ++added_clauses;
                    }
                }
            }
            aio.clear();
        }
        stats->binary_contrapositivities_learned += added_clauses;
        return added_clauses != 0;
    }

    void p_tarjan_scc_strongconnect(Lit v) {
        auto& data = tj_dfs_data[v];
        data.index = tj_index;
        data.lowlink = tj_index++;
        tj_stack.push_back(v);
        tj_on_stack[v] = true;
        auto handle_neighbor = [&](Lit w) {
            if (tj_dfs_data[w].index == NIL) {
                p_tarjan_scc_strongconnect(w);
                data.lowlink = (std::min)(data.lowlink, tj_dfs_data[w].lowlink);
            } else if (tj_on_stack[w]) {
                data.lowlink = (std::min)(data.lowlink, tj_dfs_data[w].index);
            }
        };
        for (Lit w : implications_of[v]) {
            handle_neighbor(w);
        }
        for (Lit w : additional_implications_of[v]) {
            handle_neighbor(w);
        }
        if (data.lowlink == data.index) {
            // v is root of a strongly connected component
            Lit w;
            while ((w = tj_stack.back()) != v) {
                if (equality_graph.make_equal(v, w)) {
                    tj_changes += 1;
                    stats->variables_eliminated_by_equality += 1;
                }
                tj_stack.pop_back();
                tj_on_stack[w] = false;
            }
            tj_stack.pop_back();
            tj_on_stack[v] = false;
            tj_order.push_back(v);
        }
    }

    bool p_tarjan_scc_find_equalities() {
        std::size_t nl = 2 * simplify_ds->original_num_vars();
        tj_dfs_data.assign(nl, TarjanSCCDFSInfo{});
        tj_on_stack.assign(nl, false);
        tj_stack.clear();
        tj_order.clear();
        tj_index = tj_changes = 0;
        for (Lit v = 0; v < nl; ++v) {
            if (tj_dfs_data[v].index != NIL)
                continue;
            p_tarjan_scc_strongconnect(v);
        }
        return tj_changes != 0;
    }

    bool p_init_iteration() {
        arc_set.clear();
        auto& io = implications_of;
        auto& aio = additional_implications_of;
        auto clear_vector = [](std::vector<Lit>& v) { v.clear(); };
        std::for_each(io.begin(), io.end(), clear_vector);
        std::for_each(aio.begin(), aio.end(), clear_vector);
        for (auto [l1, l2] : clauses.binary_clauses()) {
            if (propagator.is_open(l1) && propagator.is_open(l2)) {
                p_add_double_arc(lit::negate(l1), l2);
            }
        }
        return detail::fix_unaries(propagator, equality_graph, stats);
    }

    bool p_push_var(Var v) {
        Lit p = lit::positive_lit(v), n = lit::negative_lit(v);
        if (!p_push_lit(p))
            return true;
        p_implications_to_reachable();
        p_analyze_trail(p);
        propagator.pop_level();
        if (!p_push_lit(n))
            return true;
        p_analyze_trail(n);
        bool result = p_infer_fixed_from_doubly_implied();
        propagator.pop_level();
        if (result)
            propagator.incorporate_or_throw();
        return result;
    }

    void p_implications_to_reachable() {
        reachable.clear();
        for (Lit l : detail::implied_literals(propagator)) {
            reachable.insert(l);
        }
    }

    bool p_infer_fixed_from_doubly_implied() {
        bool result = false;
        for (Lit l : detail::implied_literals(propagator)) {
            if (reachable.count(l)) {
                result = true;
                clauses.add_clause(&l, &l + 1);
                stats->variables_fixed += 1;
            }
        }
        return result;
    }

    bool p_push_lit(Lit l) {
        if (!propagator.push_level(l)) {
            propagator.resolve_or_throw();
            stats->conflict_clauses_learned += 1;
            stats->variables_fixed += 1;
            stats->failed_literals += 1;
            equality_graph.make_false(l);
            return false;
        }
        return true;
    }

    void p_analyze_trail(Lit pushed) {
        auto iter_lits = propagator.current_level_begin() + 1,
             end = propagator.get_trail().end();
        auto iter_reasons = propagator.current_level_reasons_begin() + 1;
        while (iter_lits < end) {
            Lit limplied = *iter_lits++;
            Reason rimplied = *iter_reasons++;
            assert(rimplied.reason_length >= 2);
            if (rimplied.reason_length == 2)
                continue;
            p_buffer_clause(rimplied.clause);
            if (clause_analysis_buffer.size() == 2) {
                p_add_double_arc(lit::negate(clause_analysis_buffer[0]),
                                 clause_analysis_buffer[1]);
            } else {
                p_add_arc(pushed, limplied);
            }
        }
    }

    void p_add_additional_arcs() {
        for (Lit l = 0, nl = 2 * simplify_ds->original_num_vars(); l < nl; ++l)
        {
            Lit lneg = lit::negate(l);
            for (Lit o : implications_of[l]) {
                Lit oneg = lit::negate(o);
                if (!arc_set.count(std::pair<Lit, Lit>{oneg, lneg})) {
                    additional_implications_of[oneg].push_back(lneg);
                }
            }
        }
    }

    void p_buffer_clause(CRef clause) {
        clause_analysis_buffer.clear();
        auto lits = clauses.lits_of(clause);
        std::copy_if(lits.begin(), lits.end(),
                     std::back_inserter(clause_analysis_buffer), [&](Lit l) {
                         return propagator.get_decision_level(l) != 0;
                     });
    }

    void p_add_arc(Lit from, Lit to) {
        if (arc_set.emplace(from, to).second) {
            implications_of[from].push_back(to);
        }
    }

    void p_add_double_arc(Lit from, Lit to) {
        if (arc_set.emplace(from, to).second) {
            implications_of[from].push_back(to);
        }
        Lit bfrom = lit::negate(to);
        Lit bto = lit::negate(from);
        if (arc_set.emplace(bfrom, bto).second) {
            implications_of[bfrom].push_back(bto);
        }
    }

    std::vector<std::vector<CRef>> p_make_clauses_with_lit() {
        Lit nl = 2 * simplify_ds->original_num_vars();
        std::vector<std::vector<CRef>> clauses_with_lit(nl);
        auto& clauses = simplify_ds->clauses();
        for (std::size_t i = 0, n = clauses.size(); i < n; ++i) {
            const auto& clause = clauses[i];
            if (clause.size() <= 2)
                continue;
            for (Lit l : clause) {
                clauses_with_lit[l].push_back(CRef(i));
            }
        }
        return clauses_with_lit;
    }

    struct TarjanSCCDFSInfo {
        Lit index;
        Lit lowlink;

        TarjanSCCDFSInfo() noexcept {
            index = NIL;
            lowlink = 0;
        }
    };

    SimplifyDatastructure* simplify_ds;
    ClauseDB clauses;
    ClauseCounts old_counts;
    SharedDBPropagator propagator;
    EqualityGraph equality_graph;
    std::vector<std::vector<Lit>> implications_of, additional_implications_of;
    SimplificationStats stats_buffer;
    SimplificationStats* stats;
    EdgeSet arc_set;
    std::vector<Lit> clause_analysis_buffer;
    std::vector<TarjanSCCDFSInfo> tj_dfs_data;
    std::vector<bool> tj_on_stack;
    Lit tj_index, tj_changes;
    std::vector<Lit> tj_stack;
    std::vector<Lit> tj_order;
    StampSet<Lit> reachable;
};

/**
 * A very fast method to detect literals fixed at level 0,
 * and to eliminate those literals from the formula.
 * This method is mainly used for testing & debugging purposes.
 */
inline bool fix_level0_literals(SimplifyDatastructure& simplify_ds,
                                SimplificationStats* stats = nullptr) {
    ClauseDB exported{simplify_ds.extract_clause_db()};
    SharedDBPropagator propagator{&exported};
    EqualityGraph eqgraph{exported.num_vars()};
    SimplificationStats stat_buffer;
    if (!stats)
        stats = &stat_buffer;
    std::size_t count = 0;
    for (Lit l : propagator.get_trail()) {
        eqgraph.make_true(l);
        count += 1;
    }
    stats->variables_fixed += count;
    if (count != 0) {
        simplify_ds.apply_fixes_and_equalities(
            eqgraph.compute_old_to_new_map());
        return true;
    }
    return false;
}

/**
 * Detect fixed, failed and equal literals (by analyzing consequences from
 * pushing and propagating each literal at least once); this identifies most
 * equalities between literals and literals implicitly fixed to true.
 */
inline bool
detect_failed_and_equal_literals(SimplifyDatastructure& simplify_ds,
                                 SimplificationStats* stats = nullptr) {
    ImplicationGraphBuilder detector{&simplify_ds};
    if (stats)
        detector.set_stats(stats);
    if (detector.run_fixing_and_equality_iteration()) {
        return true;
    }
    return detector.run_contrapositive_and_strengthening_iteration();
}

} // namespace sammy

#endif
