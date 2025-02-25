#ifndef SAMMY_COMPRESS_BINARIES_H_INCLUDED_
#define SAMMY_COMPRESS_BINARIES_H_INCLUDED_

#include "dynamic_bitset.h"
#include "literals.h"
#include "shared_db_propagator.h"
#include "simplify_datastructure.h"
#include "stamp_set.h"
#include "thread_group.h"
#include "detect_equalities.h"

namespace sammy {

class BinaryClauseCompressor {
  public:
    explicit BinaryClauseCompressor(SimplifyDatastructure* simplify_ds) :
        simplify_ds(simplify_ds),
        m_clauses(simplify_ds->extract_clause_db()),
        m_out_edges(2 * simplify_ds->original_num_vars()),
        m_visited(2 * simplify_ds->original_num_vars()),
        m_children(2 * simplify_ds->original_num_vars()),
        m_remove(2 * simplify_ds->original_num_vars())
    {}

    void compute_binary_clause_graph() {
        for(const auto& cl : simplify_ds->clauses()) {
            if(cl.size() != 2) continue;
            Lit l1 = cl[0], l2 = cl[1];
            m_out_edges[lit::negate(l1)].push_back(l2);
            m_out_edges[lit::negate(l2)].push_back(l1);
        }
    }

    void transitively_reduce_dag() {
        std::vector<std::pair<Lit,Lit>> binaries;
        const Lit nl = 2 * simplify_ds->original_num_vars();
        for(Lit l = 0; l < nl; ++l) {
            if(simplify_ds->is_eliminated(lit::var(l))) { ++l; continue; }
            m_visited.clear();
            m_children.clear();
            m_remove.clear();
            const auto& succs = m_out_edges[l];
            std::for_each(succs.begin(), succs.end(), [&] (Lit succ) {m_children.insert(succ);});
            for(Lit succ : succs) {
                if(!m_visited.check_insert(succ)) {
                    m_remove.insert(succ);
                    continue;
                }
                m_dfs_stack.push_back(succ);
                while(!m_dfs_stack.empty()) {
                    Lit curr = m_dfs_stack.back();
                    m_dfs_stack.pop_back();
                    for(Lit next : m_out_edges[curr]) {
                        if(!m_visited.check_insert(next)) {
                            if(m_children.count(next)) {
                                m_remove.insert(next);
                            }
                        } else {
                            m_dfs_stack.push_back(next);
                        }
                    }
                }
            }
            for(Lit succ : succs) {
                if(!m_remove.count(succ)) {
                    Lit l1 = lit::negate(l), l2 = succ;
                    binaries.emplace_back((std::min)(l1,l2), (std::max)(l1,l2));
                }
            }
        }
        std::sort(binaries.begin(), binaries.end());
        binaries.erase(std::unique(binaries.begin(), binaries.end()), binaries.end());
        simplify_ds->replace_all_binaries(binaries);
    }

    void compute_transitive_implication_graph(ThreadGroup<void>* tgroup) {
        SharedDBPropagator propagator_ctx{&m_clauses};
        tgroup->parallel_foreach_iterator(Lit(0), Lit(2 * simplify_ds->original_num_vars()), propagator_ctx,
            [&] (SharedDBPropagator& propagator, Lit l) {
                if(simplify_ds->is_eliminated(lit::var(l)) || !propagator.is_open(l)) {
                    return;
                }
                if(!propagator.push_level(l)) {
                    throw std::logic_error("BinaryClauseCompressor should only be used after simplification!");
                }
                auto& out_edges = m_out_edges[l];
                auto implied = detail::implied_literals(propagator);
                std::copy(implied.begin(), implied.end(), std::back_inserter(out_edges));
                propagator.pop_level();
            }
        );
    }

    void transitive_reduction() {
        std::vector<std::pair<Lit,Lit>> binaries;
        for(Lit l = 0, nl = 2 * simplify_ds->original_num_vars(); l != nl; ++l) {
            if(simplify_ds->is_eliminated(lit::var(l))) continue;
            m_visited.clear();
            for(Lit succ : m_out_edges[l]) {
                for(Lit succ2 : m_out_edges[succ]) {
                    m_visited.insert(succ2);
                }
            }
            for(Lit succ : m_out_edges[l]) {
                if(!m_visited.count(succ)) {
                    Lit l1 = lit::negate(l);
                    Lit l2 = succ;
                    binaries.emplace_back(std::min(l1, l2), std::max(l1, l2));
                }
            }
        }
        std::sort(binaries.begin(), binaries.end());
        binaries.erase(std::unique(binaries.begin(), binaries.end()), binaries.end());
        simplify_ds->replace_all_binaries(binaries);
    }

  private:
    using Row = DynamicBitset;
    using Matrix = std::vector<Row>;
    SimplifyDatastructure *simplify_ds;
    ClauseDB m_clauses;
    std::vector<std::vector<Lit>> m_out_edges;
    std::vector<Lit> m_dfs_stack;
    StampSet<Lit, std::uint16_t> m_visited;
    StampSet<Lit, std::uint16_t> m_children;
    StampSet<Lit, std::uint16_t> m_remove;
};

}

#endif
