#ifndef SAMMY_RANDSAT_H_INCLUDED_
#define SAMMY_RANDSAT_H_INCLUDED_

#include "shared_db_propagator.h"
#include "error.h"
#include "rng.h"
#include <chrono>

namespace sammy {

struct RandSATStats {
    double solve_seconds = 0.0;
    std::size_t propagations = 0;
    std::size_t conflicts = 0;
    std::size_t total_learned_clause_length = 0;
    std::size_t conflict_level_sum = 0;
    std::size_t total_decisions = 0;

    template<typename OutObject>
    void export_to(OutObject& object) const
    {
        object["solve_time"] = solve_seconds;
        object["propagations"] = propagations;
        object["conflicts"] = conflicts;
        object["total_learned_clause_length"] = total_learned_clause_length;
        object["conflict_level_sum"] = conflict_level_sum;
        object["total_decisions"] = total_decisions;
    }
};

template<typename RNGType>
class RandSAT {
  public:
    struct TrivialClauseHandler {
        bool new_unary(Lit) { return true; }
        bool new_binary(Lit, Lit) { return true; }
        bool new_clause(ClauseDB::Lits) { return true; }
    };

    struct RecordSizeClauseHandler {
        bool new_unary(Lit) {
            that->stats->total_learned_clause_length += 1; 
            return true;
        }

        bool new_binary(Lit, Lit) {
            that->stats->total_learned_clause_length += 2;
            return true;
        }
        
        bool new_clause(ClauseDB::Lits lits) {
            that->stats->total_learned_clause_length += std::distance(lits.begin(), lits.end());
            return true;
        }

        RandSAT* that;
    };

    explicit RandSAT(RNGType& rng, const ClauseDB& clauses) : 
        rng(&rng),
        clauses(clauses),
        view(&this->clauses),
        propagator(&this->clauses),
        stats(&stats_buffer),
        literal_dist(0, 2 * clauses.num_vars() - 1),
        handler{this}
    {
        TrivialClauseHandler trivial;
        view.handle_new_clauses(trivial);
    }

    void set_stats(RandSATStats* stats) noexcept {
        this->stats = stats;
    }

    std::vector<bool> solve() {
        auto before = std::chrono::steady_clock::now();
        while(p_push_next_choice()) { stats->total_decisions += 1; }
        auto after = std::chrono::steady_clock::now();
        stats->solve_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(after-before).count();
        std::vector<bool> result(clauses.num_vars(), false);
        for(Lit l : propagator.get_trail()) {
            result[lit::var(l)] = !lit::negative(l);
        }
        return result;
    }

  private:
    bool p_push_next_choice() {
        auto c = p_pick_next_choice();
        if(!c) return false;
        bool presult = propagator.push_level(*c);
        stats->propagations += (propagator.get_trail().end() - propagator.current_level_begin()) - 1;
        if(!presult) {
            stats->conflicts += 1;
            stats->conflict_level_sum += propagator.get_current_level();
            propagator.resolve_or_throw();
            view.handle_new_clauses(handler);
        }
        return true;
    }

    std::optional<Lit> p_pick_next_choice() {
        if(propagator.get_trail().size() == clauses.num_vars()) return std::nullopt;
        for(int i = 0; i < 32; ++i) {
            Lit l = literal_dist(*rng);
            if(propagator.is_open(l)) {
                return l;
            }
        }
        open_buffer.clear();
        for(Var v = 0, nv = clauses.num_vars(); v < nv; ++v) {
            if(propagator.is_open(lit::positive_lit(v))) open_buffer.push_back(v);
        }
        std::uniform_int_distribution<std::size_t> range{0, open_buffer.size() - 1};
        Var open_var = open_buffer[range(*rng)];
        return std::uniform_int_distribution<int>{0,1}(*rng) ? 
            lit::positive_lit(open_var) : lit::negative_lit(open_var);
    }

    RNGType* rng;
    ClauseDB clauses;
    ClauseDBView view;
    SharedDBPropagator propagator;
    RandSATStats stats_buffer;
    RandSATStats* stats;
    std::uniform_int_distribution<Lit> literal_dist;
    std::vector<Var> open_buffer;
    RecordSizeClauseHandler handler;
};

template<class RNGType> RandSAT(RNGType& rng, const ClauseDB& clauses) -> RandSAT<RNGType>;

inline
std::vector<bool> randsolve(const ClauseDB& clauses) {
    RandSAT solver{sammy::rng(), clauses};
    return solver.solve();
}

}

#endif
