#ifndef SAMMY_ELIMINATE_SUBSUMED_H_INCLUDED_
#define SAMMY_ELIMINATE_SUBSUMED_H_INCLUDED_

#include "error.h"
#include "literals.h"
#include "simplification_stats.h"
#include "stamp_set.h"

namespace sammy {

template <typename ClauseType> class SubsumptionChecker {
  public:
    SubsumptionChecker(std::vector<ClauseType>& clauses, Var n_all)
        : m_nv(n_all), m_nl(2 * n_all), m_clauses(clauses), m_in_clause(m_nl),
          m_watching_clauses(m_nl), stats(&stats_buffer) {
        p_init_watches();
    }

    void set_stats(SimplificationStats* stats) noexcept { this->stats = stats; }

    void remove_subsumed() {
        // Will empty subsumed clauses and in a second step remove
        // all empty clauses.
        for (CRef c = 0, n = m_clauses.size(); c < n; ++c) {
            p_empty_if_subsumed(c);
        }
        auto deleted_begin =
            std::remove_if(m_clauses.begin(), m_clauses.end(),
                           [](const ClauseType& cl) { return cl.empty(); });
        stats->clauses_subsumed += std::size_t(m_clauses.end() - deleted_begin);
        m_clauses.erase(deleted_begin, m_clauses.end());
    }

  private:
    /**
     *  Walk the watch-list for literal `l` and:
     *    * detect if any watched clause D subsumes the current clause C,
     *      in which case we stop and return true;
     *    * otherwise, shift each watched clause D to a new watch literal
     *      (one not in C) or drop it if D is already empty.
     *
     *  This keeps each D/l pair visited exactly once across the entire
     *  elimination pass, giving near-linear overall behavior.
     */
    bool p_walk_watch_list(CRef index, Lit l) {
        auto& watch_list    = m_watching_clauses[l];
        auto  write_it      = watch_list.begin();
        const auto read_end = watch_list.end();
        bool  subsumed      = false;

        // 1) Scan every watcher D of literal l
        for (auto read_it = watch_list.begin(); read_it != read_end; ++read_it) {
            CRef watcher = *read_it;

            // 1a) never drop our own watch
            if (watcher == index) {
                *write_it++ = watcher;
                continue;
            }

            const auto& clauseD = m_clauses[watcher];
            // 1b) already-eliminated clauses vanish
            if (clauseD.empty())
                continue;

            // 1c) find a new literal in D that's not in the current clause
            auto replacement_it = std::find_if(
                clauseD.begin(), clauseD.end(),
                [&](Lit lit) { return !m_in_clause.count(lit); }
            );

            if (replacement_it == clauseD.end()) {
                // 1d) no replacement => D subset-of C => D subsumes C
                subsumed = true;
                // copy remaining watchers to the write position
                // this will overwrite all shifted or eliminated clauses
                // and leave the rest unchanged
                write_it = std::copy(read_it, read_end, write_it);
                break;
            } else {
                // 1e) shift D to watch the new literal
                m_watching_clauses[*replacement_it].push_back(watcher);
            }
        }

        // 2) erase everything we didn't rewrite
        watch_list.erase(write_it, read_end);
        return subsumed;
    }

    void p_empty_if_subsumed(CRef index) {
        ClauseType& clause = m_clauses[index];
        // this is super fast thanks to the stamp set.
        m_in_clause.assign(clause.begin(), clause.end());
        for (Lit l : clause) {
            const auto is_subsumed = p_walk_watch_list(index, l);
            if (is_subsumed) {
                // Normally, an empty clause indicates infeasibility, but here
                // we just use it to indicate that the clause has been subsumed
                // and will delete it later.
                clause.clear();
                return;
            }
        }
    }

    void p_init_watches() {
        for (std::size_t ci = 0, cn = m_clauses.size(); ci < cn; ++ci) {
            const auto& cl = m_clauses[ci];
            if (cl.empty()) {
                // An empty clause indicates a bad formula. We would not want to
                // continue, though it is likely that this would have already been
                // detected previously.
                throw std::logic_error("Empty clause in formula for subsumption.");
            }
            m_watching_clauses[cl[0]].push_back(CRef(ci));
        }
    }

    Var m_nv;
    Lit m_nl;
    std::vector<ClauseType>& m_clauses;
    StampSet<Lit, std::uint16_t> m_in_clause;
    std::vector<std::vector<CRef>> m_watching_clauses;
    SimplificationStats stats_buffer;
    SimplificationStats* stats;
};

template <typename ClauseType>
inline void eliminate_subsumed(std::vector<ClauseType>& clauses, Var n_all,
                               SimplificationStats* stats = nullptr) {
    SubsumptionChecker<ClauseType> subsumption_checker{clauses, n_all};
    if (stats)
        subsumption_checker.set_stats(stats);
    subsumption_checker.remove_subsumed();
}

} // namespace sammy

#endif
