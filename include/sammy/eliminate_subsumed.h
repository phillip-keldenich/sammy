#ifndef SAMMY_ELIMINATE_SUBSUMED_H_INCLUDED_
#define SAMMY_ELIMINATE_SUBSUMED_H_INCLUDED_

#include "error.h"
#include "literals.h"
#include "simplification_stats.h"
#include "stamp_set.h"

namespace sammy {

template<typename ClauseType>
class SubsumptionChecker {
  public:
    SubsumptionChecker(std::vector<ClauseType>& clauses, Var n_all) :
        m_nv(n_all),
        m_nl(2 * n_all),
        m_clauses(clauses),
        m_in_clause(m_nl),
        m_watching_clauses(m_nl),
        stats(&stats_buffer)
    {
        p_init_watches();
    }

    void set_stats(SimplificationStats* stats) noexcept {
        this->stats = stats;
    }

    void remove_subsumed() {
        for(CRef c = 0, n = m_clauses.size(); c < n; ++c) {
            p_empty_if_subsumed(c);
        }
        auto deleted_begin = std::remove_if(m_clauses.begin(), m_clauses.end(), 
                                            [] (const ClauseType& cl) { return cl.empty(); });
        stats->clauses_subsumed += std::size_t(m_clauses.end() - deleted_begin);
        m_clauses.erase(deleted_begin, m_clauses.end());
    }

  private:
    bool p_walk_watch_list(CRef index, Lit l) {
        auto& watch_list = m_watching_clauses[l];
        auto end = watch_list.end();
        auto out = watch_list.begin();
        bool subsumed = false;
        for(auto in = watch_list.begin(); in != end; ++in) {
            CRef cother = *in;
            // we cannot subsume ourself. stay in the watch list.
            if(cother == index) { *out++ = cother; continue; }
            const ClauseType& other_lits = m_clauses[cother];
            // subsumed clauses do not participate in subsumption anymore;
            // they are dropped from watch lists without replacement when we
            // encounter them here.
            if(other_lits.empty()) { continue; }
            // find replacement watch (must not be in the current clause).
            auto replacement = std::find_if(other_lits.begin(), other_lits.end(), [&] (Lit l) {
                return !m_in_clause.count(l);
            });
            if(replacement == other_lits.end()) {
                // cother subsumes us.
                subsumed = true;
                // copy remaining watching clauses.
                out = std::copy(in, end, out);
                break;
            } else {
                // cother does not subsume us.
                m_watching_clauses[*replacement].push_back(cother);
            }
        }
        // trim watch list
        watch_list.erase(out, end);
        return subsumed;
    }

    void p_empty_if_subsumed(CRef index) {
        ClauseType& clause = m_clauses[index];
        m_in_clause.assign(clause.begin(), clause.end());
        for(Lit l : clause) {
            if(p_walk_watch_list(index, l)) {
                clause.clear();
                return;
            }
        }
    }

    void p_init_watches() {
        for(std::size_t ci = 0, cn = m_clauses.size(); ci < cn; ++ci) {
            const auto& cl = m_clauses[ci];
            m_watching_clauses[cl[0]].push_back(CRef(ci));
        }
    }

    Var m_nv;
    Lit m_nl;
    std::vector<ClauseType>& m_clauses;
    StampSet<Lit, std::uint16_t> m_in_clause;
    std::vector<std::vector<CRef>> m_watching_clauses;
    SimplificationStats stats_buffer;
    SimplificationStats *stats;
};

template<typename ClauseType>
inline void eliminate_subsumed(std::vector<ClauseType>& clauses,
                               Var n_all, SimplificationStats* stats = nullptr) 
{
    SubsumptionChecker<ClauseType> subsumption_checker{clauses, n_all};
    if(stats) subsumption_checker.set_stats(stats);
    subsumption_checker.remove_subsumed();
}

}

#endif
