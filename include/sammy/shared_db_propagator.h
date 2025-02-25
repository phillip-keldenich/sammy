#ifndef SAMMY_SHARED_DB_PROPAGATOR_H_INCLUDED_
#define SAMMY_SHARED_DB_PROPAGATOR_H_INCLUDED_

#include "clause_db.h"
#include "literals.h"
#include "error.h"
#include <boost/circular_buffer.hpp>
#include <cassert>
#include <exception>
#include <iostream>
#include <stdexcept>

namespace sammy {
template <typename T>
using RingBuffer = boost::circular_buffer_space_optimized<T>;

/**
 * @brief A reason for a propagated literal.
 * Its either a decision (reason_length == 0),
 * a unary clause (reason_length == 1, clause in literals[0]),
 * a binary clause (reason_length == 2, clause in literals),
 * or a longer clause (clause referred to by clause).
 */
struct Reason {
    struct Decision {};
    struct Unary {
        Lit lit;
    };
    struct Binary {
        Lit lit1, lit2;
    };
    struct Clause {
        std::uint32_t length;
        CRef clause;
    };

    /* implicit */ Reason(Decision) noexcept : reason_length(0) {}

    /* implicit */ Reason(Unary unary) noexcept : reason_length(1) {
        literals[0] = unary.lit;
    }

    /* implicit */ Reason(Binary b) noexcept : reason_length(2) {
        literals[0] = b.lit1;
        literals[1] = b.lit2;
    }

    /* implicit */ Reason(Clause c) noexcept : reason_length(c.length) {
        clause = c.clause;
    }

    std::uint32_t reason_length;
    union {
        CRef clause;
        Lit literals[2];
    };

    ClauseDB::Lits lits(const ClauseDB& db) const noexcept {
        switch (reason_length) {
        case 0:
            return {nullptr, nullptr};
        case 1:
            return {+literals, literals + 1};
        case 2:
            return {+literals, literals + 2};
        default:
            return db.lits_of(clause);
        }
    }
};

struct Watcher {
    Lit blocker;
    CRef watch_info_offset;
    CRef clause;
};

struct WatchInfo {
    Lit watched[2];
};

class VariableState {
    std::int32_t value_code{-1};
    std::uint32_t stamp{0};
    std::uint32_t trail_pos{NIL};

  public:
    std::uint32_t get_trail_pos() const noexcept { return trail_pos; }

    std::uint32_t get_stamp() const noexcept { return stamp; }

    void stamp_with(std::uint32_t v) noexcept { stamp = v; }

    void assign(std::uint32_t tpos, Lit ltrue, std::int32_t level) {
        value_code = (level << 1) + std::int32_t(ltrue & 1);
        trail_pos = tpos;
    }

    std::int32_t level() const noexcept { return value_code >> 1; }

    std::int32_t state(Lit l) const noexcept {
        return (value_code >> 31) | (1 & (~std::int32_t(l) ^ value_code));
    }

    void make_open() noexcept { value_code = -1; }

    bool is_open() const noexcept { return value_code < 0; }

    bool is_false() const noexcept { return (value_code << 31) & ~value_code; }

    bool is_true() const noexcept { return !(value_code & 1); }

    bool is_false(Lit literal) const noexcept {
        // literal is true if its last bit
        // matches the value code's last bit
        // and its not open
        return ~(value_code >> 31) & 1 & (std::int32_t(literal) ^ value_code);
    }

    bool is_true(Lit literal) const noexcept {
        // literal is true if its last bit
        // matches the value code's last bit
        // and its not open
        return ~(value_code >> 31) & 1 & (~std::int32_t(literal) ^ value_code);
    }

    bool is_open_or_true(Lit literal) const noexcept {
        return ((value_code >> 31) | (~std::int32_t(literal) ^ value_code)) & 1;
    }
};

class LevelInfo {
    std::uint32_t trail_pos;
    std::uint32_t stamp{0};

  public:
    std::uint32_t get_stamp() const noexcept { return stamp; }

    void stamp_with(std::uint32_t v) noexcept { stamp = v; }

    std::uint32_t level_begin() const noexcept { return trail_pos; }

    explicit LevelInfo(std::uint32_t trail_pos) noexcept
        : trail_pos(trail_pos) {}
};

class SharedDBPropagator {
    friend class ClauseDBView;

    // View that manages the part of the clauses
    // in the database that we have already 'seen'.
    // A clause that we have seen is incorporated
    // in the watch lists etc...
    ClauseDBView view;

    // Clauses that we have not yet 'seen' from
    // our view, but that we have added watches for
    // because we learnt them (and added them to the DB).
    // Stored in ascending order of references (older clauses first).
    // Cleared when we 'see' the clause from our view.
    RingBuffer<CRef> additional_watches;

    // Additional binary clauses that we have not yet 'seen' from
    // our view, but that we want to watch because we learnt them (and added
    // them to the DB). Stored in chronologically ascending order for each
    // literal. Cleared when we 'see' the clause from our view.
    std::vector<RingBuffer<Lit>> additional_binaries;

    // The state of our variables.
    std::vector<VariableState> variables;
    // For each literal, a list of watchers.
    std::vector<std::vector<Watcher>> watchers;
    // For each pair of watchers, we have a watch_info
    // structure with the watched literals, since we
    // cannot do the typical MiniSat trick and reorder clauses.
    std::vector<WatchInfo> watch_info;

    // The literals on the trail.
    std::vector<Lit> trail_lits;
    // The reasons on the trail.
    std::vector<Reason> trail_reasons;

    // The levels of the trail.
    std::vector<LevelInfo> levels;

    // A buffer for the literals of the clause we
    // are currently learning.
    std::vector<Lit> learn_buffer;

    // The index of the next literal to propagate on.
    std::size_t trail_queue_head{0};

    // The reason for the current conflict.
    Reason conflict_reason{Reason::Decision{}};
    // The conflict literal of the current conflict (its negation is in the
    // trail).
    Lit conflict_lit{NIL};
    // A stamp counter for marking variables (conflict resolution/redundancy
    // information).
    std::uint32_t stamp_counter{0};
    // Whether we have a current conflict.
    bool conflicting{false};
    // Buffer for supporting decisions of a literal.
    std::vector<std::pair<std::int32_t,Lit>> supporting_decision_buffer;

    /**
     * Assign the given literal to true at decision level 0.
     */
    bool p_assign_at_0(Lit forced_true) {
        VariableState& vstate = variables[lit::var(forced_true)];
        if (vstate.is_open()) {
            vstate.assign(trail_lits.size(), forced_true, 0);
            trail_lits.push_back(forced_true);
            trail_reasons.push_back(Reason::Unary{forced_true});
        } else {
            if (vstate.is_false(forced_true)) {
                conflicting = true;
                return false;
            }
        }
        return true;
    }

    /**
     * Assign the given literal to true at the given level,
     * using the given reason.
     */
    template <typename ReasonType>
    void p_assign_at(VariableState& vstate, std::int32_t level, Lit literal,
                     ReasonType&& reason) {
        vstate.assign(trail_lits.size(), literal, level);
        trail_lits.push_back(literal);
        trail_reasons.emplace_back(std::forward<ReasonType>(reason));
    }

    /**
     * Initialize level 0 with the unary clauses from the clause database.
     */
    void p_init_unaries() {
        struct NewUnaryOnConstruction {
            SharedDBPropagator* that;
            bool new_unary(Lit forced_true) {
                return that->p_assign_at_0(forced_true);
            }
        };
        NewUnaryOnConstruction handler{this};
        view.handle_new_unaries(handler);
    }

    /**
     * Initialize the watches for binary clauses (handle differently compared to
     * larger clauses).
     */
    void p_init_binary_watches() {
        struct NewBinaryClauseOnConstruction {
            SharedDBPropagator* that;
            bool new_binary(Lit l1, Lit l2) {
                VariableState& v1 = that->variables[lit::var(l1)];
                VariableState& v2 = that->variables[lit::var(l2)];
                if (v1.is_false(l1)) {
                    if (v2.is_open()) {
                        that->db().add_clause(&l2, &l2 + 1);
                    }
                    return that->p_assign_at_0(l2);
                } else if (v2.is_false(l2)) {
                    if (v1.is_open()) {
                        that->db().add_clause(&l1, &l1 + 1);
                    }
                    return that->p_assign_at_0(l1);
                }
                return true;
            }
        };
        NewBinaryClauseOnConstruction handler{this};
        view.handle_new_binary_clauses(handler);
    }

    /**
     * Initialize the watches for longer clauses (length > 2).
     */
    void p_init_watches() {
        watchers.resize(2 * db().num_vars());
        struct NewLongClauseOnConstruction {
            SharedDBPropagator* that;
            bool new_clause(ClauseDB::Lits literals) {
                std::int32_t nws = 0;
                Lit ws[2];
                for (Lit l : literals) {
                    VariableState& vstate = that->variables[lit::var(l)];
                    auto s = vstate.state(l);
                    if (s == -1) {
                        if (nws < 2) {
                            ws[nws++] = l;
                        }
                    } else if (s == 1) {
                        nws = -1;
                        break;
                    }
                }
                if (nws == -1)
                    return true; // satisfied at level 0
                if (nws == 0) {  // violated at level 0
                    that->conflicting = true;
                    return false;
                }
                if (nws == 1) { // forcing at level 0
                    that->db().add_clause(&ws[0], &ws[1]);
                    return that->p_assign_at_0(ws[0]);
                }
                WatchInfo winfo{{ws[0], ws[1]}};
                std::uint32_t wisize = that->watch_info.size();
                CRef cr = that->db().cref_of(literals);
                assert(cr != NIL);
                Watcher w1{ws[1], wisize, cr}, w2{ws[0], wisize, cr};
                that->watch_info.push_back(winfo);
                that->watchers[ws[0]].push_back(w1);
                that->watchers[ws[1]].push_back(w2);
                return true;
            }
        };
        NewLongClauseOnConstruction handler{this};
        view.handle_new_long_clauses(handler);
        if (!conflicting)
            p_init_binary_watches();
    }

    /**
     * Propagate using a binary clause.
     * Return false on conflict.
     */
    bool p_propagate_binary(Lit lfalse, Lit other, std::int32_t level) {
        Lit v = lit::var(other);
        VariableState& vs = variables[v];
        if (vs.is_open()) {
            p_assign_at(vs, level, other, Reason::Binary{lfalse, other});
        } else {
            if (vs.is_false(other)) {
                conflicting = true;
                conflict_reason = Reason::Binary{lfalse, other};
                conflict_lit = other;
                return false;
            }
        }
        return true;
    }

    bool p_propagate_through_binaries(Lit ltrue) {
        Lit lfalse = lit::negate(ltrue);
        auto level = std::int32_t(levels.size() - 1);
        for (Lit other : view.binary_partners_of(lfalse)) {
            if (!p_propagate_binary(lfalse, other, level))
                return false;
        }
        if (!additional_binaries.empty()) {
            for (Lit other : additional_binaries[lfalse]) {
                if (!p_propagate_binary(lfalse, other, level))
                    return false;
            }
        }
        return true;
    }

    template <typename WatcherIterator>
    bool p_has_true_blocker(WatcherIterator& watcher_in,
                            WatcherIterator& watcher_out) {
        Lit blocker = watcher_in->blocker;
        Lit bvar = lit::var(blocker);
        if (variables[bvar].is_true(blocker)) {
            *watcher_out = *watcher_in;
            ++watcher_out, ++watcher_in;
            return true;
        }
        return false;
    }

    Lit p_find_replacement(ClauseDB::Lits clause_lits, Lit other_watched) {
        for (Lit l : clause_lits) {
            if (l == other_watched)
                continue;
            Lit v = lit::var(l);
            VariableState& vs = variables[v];
            if (vs.is_open_or_true(l)) {
                return l;
            }
        }
        return NIL;
    }

    bool p_propagate_through_longer(Lit ltrue) {
        Lit lfalse = lit::negate(ltrue);
        auto level = std::int32_t(levels.size() - 1);
        std::vector<Watcher>& ws = watchers[lfalse];
        auto watcher_in = ws.begin(), watcher_out = ws.begin(),
             watcher_end = ws.end();
        while (watcher_in != watcher_end) {
            assert(watcher_in->clause != NIL);
            if (p_has_true_blocker(watcher_in, watcher_out))
                continue;
            WatchInfo& winfo = watch_info[watcher_in->watch_info_offset];
            if (winfo.watched[0] == lfalse) {
                std::swap(winfo.watched[0], winfo.watched[1]);
            }
            Lit new_blocker = winfo.watched[0];
            Watcher new_watcher{new_blocker, watcher_in->watch_info_offset,
                                watcher_in->clause};
            VariableState& nbstate = variables[lit::var(new_blocker)];
            if (new_blocker != watcher_in->blocker) {
                if (nbstate.is_true(new_blocker)) {
                    *watcher_out = new_watcher;
                    ++watcher_out, ++watcher_in;
                    continue;
                }
            }
            auto clause_lits = db().lits_of(watcher_in->clause);
            Lit replacement = p_find_replacement(clause_lits, new_blocker);
            if (replacement == NIL) {
                auto rs = Reason::Clause{
                    std::uint32_t(clause_lits.end() - clause_lits.begin()),
                    watcher_in->clause};
                if (nbstate.is_false(new_blocker)) {
                    conflict_lit = new_blocker;
                    conflict_reason = rs;
                    conflicting = true;
                    for (; watcher_in != watcher_end;
                         ++watcher_in, ++watcher_out)
                    {
                        *watcher_out = *watcher_in;
                    }
                    break;
                } else {
                    p_assign_at(nbstate, level, new_blocker, rs);
                    *watcher_out = new_watcher;
                    ++watcher_out, ++watcher_in;
                    continue;
                }
            }
            winfo.watched[1] = replacement;
            ++watcher_in;
            assert(new_watcher.clause != NIL);
            watchers[replacement].push_back(new_watcher);
        }
        ws.erase(watcher_out, watcher_end);
        return !conflicting;
    }

    bool p_propagate(Lit ltrue) {
        if (!p_propagate_through_binaries(ltrue))
            return false;
        return p_propagate_through_longer(ltrue);
    }

    std::uint32_t p_increase_stamp() noexcept {
        if (stamp_counter >= std::numeric_limits<std::uint32_t>::max() - 6) {
            for (VariableState& vs : variables) {
                vs.stamp_with(0);
            }
            for (LevelInfo& lvl : levels) {
                lvl.stamp_with(0);
            }
            stamp_counter = 0;
        }
        stamp_counter += 3;
        return stamp_counter;
    }

    void p_stamp_level(std::int32_t level) {
        LevelInfo& li = levels[level];
        if (li.get_stamp() < stamp_counter) {
            li.stamp_with(stamp_counter);
        } else {
            li.stamp_with(stamp_counter + 1);
        }
    }

    std::uint32_t p_stamp_and_count(std::int32_t level,
                                    ClauseDB::Lits literals) {
        std::uint32_t count = 0;
        for (Lit l : literals) {
            Lit v = lit::var(l);
            VariableState& vs = variables[v];
            std::int32_t vlvl = vs.level();
            if (vlvl >= level) {
                if (vs.get_stamp() >= stamp_counter)
                    continue;
                ++count;
                vs.stamp_with(stamp_counter);
            } else {
                if (vlvl <= 0)
                    continue;
                std::uint32_t vstamp = vs.get_stamp();
                if (vstamp < stamp_counter) {
                    p_stamp_level(vlvl);
                    learn_buffer.push_back(l);
                    vs.stamp_with(stamp_counter);
                }
            }
        }
        return count;
    }

    std::uint32_t p_stamp_and_count(std::int32_t level, Reason reason) {
        if (reason.reason_length <= 2) {
            auto lits = ClauseDB::Lits{&reason.literals[0],
                                       &reason.literals[reason.reason_length]};
            return p_stamp_and_count(level, lits);
        } else {
            return p_stamp_and_count(level, db().lits_of(reason.clause));
        }
    }

    void p_bfs_reasons(std::uint32_t current_stamp) {
        std::size_t lbpos = 0;
        while(lbpos < learn_buffer.size()) {
            Lit next = learn_buffer[lbpos++];
            std::size_t tindex = get_trail_index(next);
            if(trail_reasons[tindex].reason_length == 0) {
                supporting_decision_buffer.emplace_back(get_decision_level(next), next);
            } else {
                Reason reason = trail_reasons[tindex];
                for(Lit lr : reason.lits(db())) {
                    if(lr != next) {
                        Var v = lit::var(lr);
                        if(variables[v].get_stamp() != current_stamp) {
                            variables[v].stamp_with(current_stamp);
                            learn_buffer.push_back(lit::negate(lr));
                        }
                    }
                }
            }
        }
    }

    bool p_is_redundant(Lit v) {
        auto& vs = variables[v];
        auto s = vs.get_stamp();
        if (s == stamp_counter + 1)
            return true;
        if (s == stamp_counter + 2)
            return false;
        auto tloc = variables[v].get_trail_pos();
        const Reason& r = trail_reasons[tloc];
        if (r.reason_length == 0) {
            vs.stamp_with(stamp_counter + 2);
            return false;
        }
        auto reason_lits = r.lits(db());
        for (Lit rl : reason_lits) {
            Lit rv = lit::var(rl);
            if (rv == v)
                continue;
            auto rlvl = variables[rv].level();
            if (rlvl == 0)
                continue;
            auto rvs = variables[rv];
            auto rs = rvs.get_stamp();
            if (rs == stamp_counter + 2)
                return false;
            if (rs < stamp_counter) {
                if (levels[rlvl].get_stamp() < stamp_counter ||
                    !p_is_redundant(rv))
                {
                    rvs.stamp_with(stamp_counter + 2);
                    return false;
                }
            }
        }
        vs.stamp_with(stamp_counter + 1);
        return true;
    }

    void p_filter_redundancies() {
        // move the conflict-level literal to the front
        std::swap(learn_buffer.back(), learn_buffer.front());
        learn_buffer.erase(std::remove_if(learn_buffer.begin() + 1,
                                          learn_buffer.end(),
                                          [&](Lit l) {
                                              Lit v = lit::var(l);
                                              auto vlvl = variables[v].level();
                                              if (vlvl == 0)
                                                  return true;
                                              if (levels[vlvl].get_stamp() !=
                                                  stamp_counter + 1)
                                                  return false;
                                              return p_is_redundant(v);
                                          }),
                           learn_buffer.end());
    }

    std::pair<std::int32_t, Lit> p_target_level() {
        std::int32_t target_level = 0;
        Lit target_lit = learn_buffer.front();
        for (auto i = learn_buffer.begin() + 1, e = learn_buffer.end(); i != e;
             ++i)
        {
            Lit l = *i;
            Lit v = lit::var(l);
            auto lvl = variables[v].level();
            if (lvl > target_level) {
                target_level = lvl;
                target_lit = l;
            }
        }
        return {target_level, target_lit};
    }

    template <typename AssignmentHandler>
    void p_rollback_level(AssignmentHandler& handler, bool report) {
        if (levels.back().level_begin() == 0) {
            for (auto i = trail_lits.rbegin(), e = trail_lits.rend(); i != e;
                 ++i)
            {
                Lit l = *i;
                if (report)
                    handler.assignment_undone(l);
                variables[lit::var(l)].make_open();
            }
            trail_lits.clear();
            trail_reasons.clear();
        } else {
            auto current_end =
                trail_lits.begin() + (levels.back().level_begin() - 1);
            auto current_begin = trail_lits.end() - 1;
            for (; current_begin != current_end; --current_begin) {
                trail_reasons.pop_back();
                Lit l = *current_begin;
                if (report)
                    handler.assignment_undone(l);
                variables[lit::var(l)].make_open();
            }
            trail_lits.erase(current_end + 1, trail_lits.end());
        }
        levels.pop_back();
    }

    template <typename AssignmentHandler>
    std::pair<std::int32_t, Lit>
    p_jumpback_to_target(AssignmentHandler& handler) {
        auto [tlvl, tlit] = p_target_level();
        p_rollback_level(handler, false);
        while (levels.size() > std::size_t(tlvl + 1)) {
            p_rollback_level(handler, true);
        }
        trail_queue_head = trail_lits.size();
        return {tlvl, tlit};
    }

    template <typename AssignmentHandler>
    void p_handle_conflict_clause(CRef cref_if_long,
                                  AssignmentHandler& handler) {
        auto [tlvl, tlit] = p_jumpback_to_target(handler);
        Lit learned = learn_buffer.front();
        Lit lv = lit::var(learned);
        std::uint32_t len = learn_buffer.size();
        switch (len) {
        case 1:
            p_assign_at(variables[lv], tlvl, learned, Reason::Unary{learned});
            break;
        case 2: {
            Lit other = learn_buffer[1];
            p_assign_at(variables[lv], tlvl, learned,
                        Reason::Binary{learned, other});
            if (additional_binaries.empty()) {
                additional_binaries.resize(db().num_vars() * 2);
            }
            additional_binaries[learned].push_back(other);
            additional_binaries[other].push_back(other);
            break;
        }
        default: {
            p_assign_at(variables[lv], tlvl, learned,
                        Reason::Clause{len, cref_if_long});
            p_new_watch(learned, tlit, cref_if_long);
            additional_watches.push_back(cref_if_long);
            break;
        }
        }
        learn_buffer.clear();
    }

    void p_compute_conflict_clause() {
        p_increase_stamp();
        std::int32_t level(levels.size() - 1);
        std::uint32_t on_current_level =
            p_stamp_and_count(level, conflict_reason);
        auto trail_lit_iter = trail_lits.end() - 1;
        auto trail_reason_iter = trail_reasons.end() - 1;
        while (on_current_level > 1) {
            Lit l = *trail_lit_iter;
            Lit v = lit::var(l);
            if (variables[v].get_stamp() >= stamp_counter) {
                on_current_level +=
                    p_stamp_and_count(level, *trail_reason_iter);
                --on_current_level;
            }
            --trail_lit_iter;
            --trail_reason_iter;
        }
        for (;;) {
            Lit l = *trail_lit_iter;
            Lit v = lit::var(l);
            if (variables[v].get_stamp() >= stamp_counter)
                break;
            --trail_lit_iter;
            --trail_reason_iter;
        }
        learn_buffer.push_back(lit::negate(*trail_lit_iter));
        p_filter_redundancies();
    }

    void p_reset_conflict() noexcept {
        conflicting = false;
        conflict_lit = NIL;
        conflict_reason = Reason::Decision{};
    }

    void p_new_watch(Lit l1, Lit l2, CRef clause) {
        assert(clause != NIL);
        WatchInfo winfo;
        winfo.watched[0] = l1;
        winfo.watched[1] = l2;
        CRef woffs = watch_info.size();
        Watcher w1{l2, woffs, clause};
        Watcher w2{l1, woffs, clause};
        watch_info.push_back(winfo);
        watchers[l1].push_back(w1);
        watchers[l2].push_back(w2);
    }

    class L0ClauseIncorporationHandler {
        SharedDBPropagator* that;

      public:
        bool found_unsat = false;

        explicit L0ClauseIncorporationHandler(SharedDBPropagator* that) noexcept
            : that(that) {}

        bool new_unary(Lit l1) {
            if (!that->p_assign_at_0(l1)) {
                found_unsat = true;
                return false;
            }
            return true;
        }

        bool new_binary(Lit l1, Lit l2) {
            if (!that->additional_binaries.empty()) {
                auto& abin1 = that->additional_binaries[l1];
                auto& abin2 = that->additional_binaries[l2];
                if (!abin1.empty() && abin1.front() == l2) {
                    abin1.pop_front();
                    abin2.pop_front();
                    return true;
                }
            }
            if (that->is_false(l1)) {
                if (!that->p_assign_at_0(l2)) {
                    found_unsat = true;
                    return false;
                }
            } else if (that->is_false(l2)) {
                if (!that->p_assign_at_0(l1)) {
                    found_unsat = true;
                    return false;
                }
            }
            return true;
        }

        bool new_clause(ClauseDB::Lits lits) {
            CRef cref = that->db().cref_of(lits);
            auto& aw = that->additional_watches;
            if (!aw.empty() && aw.front() == cref) {
                // already watched
                aw.pop_front();
                return true;
            }
            for (Lit l : lits) {
                if (that->is_true(l)) {
                    // satisfied at level 0 (no need to watch)
                    return true;
                }
            }
            Lit lopen[2];
            int num_open = 0;
            for (Lit l : lits) {
                if (that->is_open(l)) {
                    lopen[num_open++] = l;
                    if (num_open == 2)
                        break;
                }
            }
            if (num_open == 0) {
                found_unsat = true;
                return false;
            }
            if (num_open == 1) {
                if (!that->p_assign_at_0(lopen[0])) {
                    found_unsat = true;
                    return false;
                }
                return true;
            }
            that->p_new_watch(lopen[0], lopen[1], cref);
            return true;
        }
    };



  public:
    /**
     * @brief Get the underlying database.
     *
     * @return ClauseDB&
     */
    ClauseDB& db() noexcept { return view.database(); }

    /**
     * @brief Get the underlying database.
     *
     * @return const ClauseDB&
     */
    const ClauseDB& db() const noexcept { return view.database(); }

    /**
     * @brief Construct a new SharedDBPropagator using a clause database.
     * On calls to resolve_conflicts, the propagator adds clauses to the
     * database. These clauses are implied (in non-obvious ways) by the clauses
     * already in the database.
     *
     * @param db
     */
    explicit SharedDBPropagator(ClauseDB* db)
        : view(db), variables(db->num_vars()), levels{{LevelInfo{0}}} {
        p_init_unaries();
        if (conflicting)
            return;
        p_init_watches();
        propagate();
    }

    /**
     * @brief Trigger propagation; it should not be necessary to call this
     * manually.
     *
     * @return true
     * @return false
     */
    bool propagate() {
        if (conflicting)
            return false;
        while (trail_queue_head < trail_lits.size()) {
            Lit prop = trail_lits[trail_queue_head++];
            if (!p_propagate(prop))
                return false;
        }
        return true;
    }

    /**
     * @brief Check if the given literal is assigned to true in the current
     * trail.
     *
     * @param literal
     * @return true
     * @return false
     */
    bool is_true(Lit literal) const noexcept {
        Lit v = lit::var(literal);
        return variables[v].is_true(literal);
    }

    /**
     * Incorporate a full assignment into this propagator.
     * Throws an exception if the assignment is infeasible.
     */
    template<typename FullAssignment>
    void incorporate_assignment(const FullAssignment& assignment)
    {
        for(Var v = 0, nv = db().num_vars(); v < nv; ++v) {
            Lit l = assignment[v] ? lit::positive_lit(v) : lit::negative_lit(v);
            if(is_false(l) || (is_open(l) && !push_level(l))) {
                throw std::logic_error("Assignment is infeasible!");
            }
        }
    }

    /**
     * Compute a list of [Level, Literal] pairs of decisions
     * that ultimately led to including l in the trail.
     */
    const std::vector<std::pair<std::int32_t,Lit>> &decisions_leading_to(Lit l) {
        if(conflicting) throw std::logic_error("decisions_leading_to called on propagator with conflict!");
        if(is_open(l)) throw std::logic_error("decisions_leading_to called with open literal!");
        supporting_decision_buffer.clear();
        std::size_t tindex = get_trail_index(l);
        if(trail_reasons[tindex].reason_length == 0) {
            supporting_decision_buffer.emplace_back(get_decision_level(l), l);
            return supporting_decision_buffer;
        }

        auto current = p_increase_stamp();
        Reason reason = trail_reasons[tindex];
        for(Lit lr : reason.lits(db())) {
            if(lr != l) {
                variables[lit::var(lr)].stamp_with(current);
                learn_buffer.push_back(lit::negate(lr));
            }
        }
        p_bfs_reasons(current);
        learn_buffer.clear();
        return supporting_decision_buffer;
    }

    /*
     * Compute a list of [Level, Literal] pairs of decisions
     * that ultimately led to the current conflict.
     */
    const std::vector<std::pair<std::int32_t,Lit>>& decisions_leading_to_conflict() {
        if(!conflicting) 
            throw std::logic_error("decisions_leading_to_conflict called on non-conflicting propagator!");
    
        supporting_decision_buffer.clear();
        auto current = p_increase_stamp();
        for(Lit lr : conflict_reason.lits(db())) {
            if(lr != conflict_lit) {
                variables[lit::var(lr)].stamp_with(current);
                learn_buffer.push_back(lit::negate(lr));
            }
        }
        variables[lit::var(conflict_lit)].stamp_with(current);
        Lit lc = lit::negate(conflict_lit);
        for(Lit lr : get_reason(lc).lits(db()))  {
            if(variables[lit::var(lr)].get_stamp() != current) {
                variables[lit::var(lr)].stamp_with(current);
                learn_buffer.push_back(lit::negate(lr));
            }
        }
        p_bfs_reasons(current);
        learn_buffer.clear();
        return supporting_decision_buffer;
    }

    /**
     * @brief Check if the given literal is assigned to false in the current
     * trail.
     *
     * @param literal
     * @return true
     * @return false
     */
    bool is_false(Lit literal) const noexcept {
        Lit v = lit::var(literal);
        return variables[v].is_false(literal);
    }

    /**
     * @brief Check if the given literal is unassigned/open in the current
     * trail.
     *
     * @param literal
     * @return true
     * @return false
     */
    bool is_open(Lit literal) const noexcept {
        Lit v = lit::var(literal);
        return variables[v].is_open();
    }

    /**
     * @brief Check if the given non-open literal was assigned as a decision.
     *
     * @param literal
     * @return true
     * @return false
     */
    bool is_decision(Lit literal) const noexcept {
        assert(!is_open(literal));
        auto tpos = variables[lit::var(literal)].get_trail_pos();
        return trail_reasons[tpos].reason_length == 0;
    }

    /**
     * @brief Get the decision level of the literal if it is in the trail.
     * @return The decision level, or a negative value if the literal is open.
     */
    std::int32_t get_decision_level(Lit literal) const noexcept {
        return variables[lit::var(literal)].level();
    }

    /**
     * @brief Get the reason for the literal if it is in the trail;
     *        otherwise, causes undefined behavior.
     */
    Reason get_reason(Lit literal) const noexcept {
        auto tpos = variables[lit::var(literal)].get_trail_pos();
        return trail_reasons[tpos];
    }

    /**
     * @brief Push a decision literal and create a new decision level.
     * Automatically propagates the decision and all consequences.
     * Throws an error if the decision literal is already assigned (true or
     * false).
     *
     * @param decision
     * @return true If the decision did not result in a conflict.
     * @return false If the decision resulted in a conflict.
     */
    bool push_level(Lit decision) {
        Lit dvar = lit::var(decision);
        VariableState& vstate = variables[dvar];
        if (!vstate.is_open()) {
            throw std::invalid_argument(
                "The given decision literal was already assigned!");
        }
        std::uint32_t tpos = trail_lits.size();
        std::int32_t new_level = levels.size();
        levels.emplace_back(tpos);
        p_assign_at(vstate, new_level, decision, Reason::Decision{});
        return propagate();
    }

    /**
     * @brief Pop the highest decision level without learning.
     * There is no need for a conflict to use this method.
     * On conflict, it should be preferred to call resolve_conflicts instead,
     * but pop_level will also clear the conflict.
     */
    void pop_level() {
        if (levels.size() == 1) {
            throw std::invalid_argument(
                "Trying to pop level from propagator at level 0!");
        }
        struct TrivialHandler {
            void assignment_undone(Lit) {}
        } handler;
        p_rollback_level(handler, false);
        trail_queue_head = trail_lits.size();
        if (conflicting)
            p_reset_conflict();
    }

    std::int32_t get_current_level() const noexcept {
        return std::int32_t(levels.size() - 1);
    }

    std::vector<Lit>::const_iterator current_level_begin() const noexcept {
        return trail_lits.begin() + levels.back().level_begin();
    }

    std::vector<Reason>::const_iterator current_level_reasons_begin() const noexcept {
        return trail_reasons.begin() + levels.back().level_begin();
    }

    /**
     * @brief Get an iterator to the beginning of the given level in the trail.
     *
     * @param level
     * @return std::vector<Lit>::const_iterator
     */
    std::vector<Lit>::const_iterator level_begin(std::uint32_t level) const {
        return trail_lits.begin() + levels[level].level_begin();
    }

    /**
     * @brief Get an iterator to the end of the given level in the trail.
     *
     * @param level
     * @return std::vector<Lit>::const_iterator
     */
    std::vector<Lit>::const_iterator level_end(std::uint32_t level) const {
        if (level >= levels.size() - 1)
            return trail_lits.end();
        return trail_lits.begin() + levels[level + 1].level_begin();
    }

    std::vector<Lit>::const_iterator literal_position(Lit lit) const {
        return trail_lits.begin() + variables[lit::var(lit)].get_trail_pos();
    }

    /**
     * @brief Get the list of literals that are currently assigned to true.
     *
     * @return const std::vector<Lit>&
     */
    const std::vector<Lit>& get_trail() const noexcept { return trail_lits; }

    /**
     * @brief Get the list of reasons for the literals that are currently assigned to true.
     */
    const std::vector<Reason>&  get_reasons() const noexcept { return trail_reasons; }

    /**
     * @brief Get a list of all decision literals on the trail.
     *        Creates a new vector (unlike get_trail, which returns
     *        a reference to an internal structure).
     *
     * @return std::vector<Lit>
     */
    std::vector<Lit> get_decisions() const {
        std::vector<Lit> result;
        result.reserve(levels.size() - 1);
        for (auto ilvl = levels.begin() + 1, iend = levels.end(); ilvl != iend;
             ++ilvl)
        {
            result.push_back(trail_lits[ilvl->level_begin()]);
        }
        return result;
    }

    /**
     * @brief Check whether we currently have a conflict.
     *
     * @return true
     * @return false
     */
    bool is_conflicting() const noexcept { return conflicting; }

    /**
     * @brief Get the conflict literal and reason.
     * @return std::pair<Lit, Reason>
     */
    std::pair<Lit, Reason> get_conflict() const noexcept {
        return {conflict_lit, conflict_reason};
    }

    /**
     * @brief Resolve a conflict by learning a clause and jumping back
     * to the appropriate decision level (at least one level down).
     * At least one assignment is forced on the target decision level
     * from the conflict clause.
     *
     * In any case, all assignments on the current level are undone because
     * its decision led to a conflict; this is *NOT* reported to the given
     * AssignmentHandler.
     *
     * However, all assignments on lower levels that are undone or forced by
     * this action *ARE* reported. After learning a conflict clause and jumping
     * back, we continue propagation. It is possible that another conflict
     * occurs. This conflict is also handled recursively (and now all undone
     * assignments are reported). This is repeated until we reach a state
     * without conflicts (we return true), or we reach a conflict on level 0 (we
     * return false and the formula is UNSAT).
     *
     * @tparam AssignmentHandler A type that implements methods
     * assignment_undone(Lit) and assignment_forced(Lit).
     * @param assignments The AssignmentHandler that is notified of changes.
     * @return true if at some level, we ended in a non-conflicting state.
     * @return false if we encountered a conflict at level 0, indicating
     * infeasibility.
     */
    template <typename AssignmentHandler>
    bool resolve_conflicts(AssignmentHandler& assignments) {
        if (!conflicting)
            return true;
        if (levels.size() == 1)
            return false;
        p_compute_conflict_clause();
        assert(learn_buffer.size() > 0);
        CRef cc = db().add_clause(learn_buffer.data(),
                                  learn_buffer.data() + learn_buffer.size());
        assert(learn_buffer.size() <= 2 || cc != NIL);
        p_handle_conflict_clause(cc, assignments);
        p_reset_conflict();
        std::size_t tsize = trail_queue_head;
        std::size_t lbegin = levels.back().level_begin();
        if (!propagate()) {
            for (std::size_t cpos = tsize - 1; cpos != lbegin - 1; --cpos) {
                assignments.assignment_undone(trail_lits[cpos]);
            }
            return resolve_conflicts(assignments);
        } else {
            for (auto i = trail_lits.begin() + tsize, e = trail_lits.end();
                 i != e; ++i)
            {
                assignments.assignment_forced(*i);
            }
            return true;
        }
    }

    /**
     * Compute the total memory used by this propagator, in bytes.
     */
    std::size_t total_memory_usage() const noexcept {
        std::size_t total_watcher_size =
            watchers.capacity() * sizeof(std::vector<Watcher>);
        for (const auto& w : watchers) {
            total_watcher_size += w.capacity() * sizeof(Watcher);
        }
        std::size_t additional_binary_size =
            additional_binaries.capacity() * sizeof(RingBuffer<Lit>);
        for (const auto& a : additional_binaries) {
            additional_binary_size += a.capacity() * sizeof(Lit);
        }
        return view.total_memory_usage() +
               additional_watches.capacity() * sizeof(CRef) +
               variables.capacity() * sizeof(VariableState) +
               watch_info.capacity() * sizeof(WatchInfo) +
               trail_lits.capacity() * sizeof(Lit) +
               trail_reasons.capacity() * sizeof(Reason) +
               levels.capacity() * sizeof(LevelInfo) +
               learn_buffer.capacity() * sizeof(Lit) + total_watcher_size +
               additional_binary_size + sizeof(SharedDBPropagator);
    }

    /**
     * @brief Resolve conflicts without handler.
     *
     * @return true
     * @return false
     */
    bool resolve_conflicts() {
        struct TrivialAssignmentHandler {
            void assignment_undone(Lit) const noexcept {}
            void assignment_forced(Lit) const noexcept {}
        };
        TrivialAssignmentHandler handler;
        return resolve_conflicts(handler);
    }

    /**
     * @brief Reset the propagator to level 0. Incorporates new clauses.
     * @throws UNSATError if this results in a conflict at level 0 (the formula is UNSAT).
     */
    void reset_or_throw() {
        if(is_conflicting()) resolve_or_throw();
        while(get_current_level() > 0) {
            pop_level();
        }
        incorporate_or_throw();
    }

    /**
     * @brief Resolve conflicts without handler.
     * @throws UNSATError if the conflict cannot be resolved (the formula is UNSAT).
     */
    void resolve_or_throw() {
        if(!resolve_conflicts()) throw UNSATError();
    }

    /**
     * Reset propagator to level 0.
     */
    void reset_to_zero() noexcept {
        while(get_current_level() > 0) {
            pop_level();
        }
    }

    /**
     * @brief Incorporate all new clauses.
     *        The propagator must be at level 0.
     *
     * @return true if we did not find unsatisfiability.
     * @return false if we found unsatisfiability due to the new clauses.
     */
    bool incorporate_new_clauses_at_level_0() {
        assert(levels.size() == 1);
        L0ClauseIncorporationHandler handler{this};
        view.handle_new_clauses(handler);
        conflicting = handler.found_unsat;
        if (!conflicting) {
            additional_binaries.clear();
            propagate();
        }
        return !conflicting;
    }

    /**
     * @brief Like incorporate_new_clauses_at_level_0.
     * @throws UNSATError if the formula becomes UNSAT by the added clauses.
     */
    void incorporate_or_throw() {
        if(!incorporate_new_clauses_at_level_0()) throw UNSATError();
    }

    /**
     * @brief Extract an assignment as bit-vector, where result[i] == true means
     *        that variable i (internal 0-based indexing) is set to true.
     */
    std::vector<bool> extract_assignment() const {
        const Var nv = db().num_vars();
        if(get_trail().size() != nv) {
            throw std::logic_error("Trail incomplete in extract_assignment!");
        }
        std::vector<bool> result(nv, false);
        for(Lit l : get_trail()) {
            if(!lit::negative(l)) {
                result[lit::var(l)] = true;
            }
        }
        return result;
    }

    /**
     * Get the index in the trail of the given literal.
     * Undefined behaviour if lit is open.
     */
    std::size_t get_trail_index(Lit lit) const noexcept {
        return variables[lit::var(lit)].get_trail_pos();
    }
};

} // namespace hs

#endif
