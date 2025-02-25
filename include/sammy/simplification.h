#ifndef SAMMY_SIMPLIFICATION_H_INCLUDED_
#define SAMMY_SIMPLIFICATION_H_INCLUDED_

#include "eliminate_subsumed.h"
#include "simplify_datastructure.h"
#include "detect_equalities.h"
#include "bounded_variable_elimination.h"
#include "vivify.h"
#include "compress_binaries.h"

namespace sammy {

inline
SimplifiedInstance run_simplifier(SimplifyDatastructure& simplifier, SimplificationStats& stats) {
    bool changed;
    Var nv = simplifier.original_num_vars();
    auto before = Clock::now();
    eliminate_subsumed(simplifier.clauses(), nv, &stats);
    assert(!simplifier.has_duplicate_binary_clause());
    while(detect_failed_and_equal_literals(simplifier, &stats));
    eliminate_subsumed(simplifier.clauses(), nv, &stats);
    assert(!simplifier.has_duplicate_binary_clause());
    do {
        changed = bounded_variable_elimination(simplifier, 20, &stats);
        if(changed) eliminate_subsumed(simplifier.clauses(), nv, &stats);
        assert(!simplifier.has_duplicate_binary_clause());
        changed |= vivify(simplifier, &stats);
        changed |= detect_failed_and_equal_literals(simplifier, &stats);
        if(changed) eliminate_subsumed(simplifier.clauses(), nv, &stats);
        assert(!simplifier.has_duplicate_binary_clause());
        stats.simplification_rounds += 1;
    } while(changed);
    BinaryClauseCompressor comp{&simplifier};
    comp.compute_binary_clause_graph();
    comp.transitively_reduce_dag();
    SimplifiedInstance simplified = simplifier.compress();
    auto after = Clock::now();
    stats.simplification_time = seconds_between(before, after);
    return simplified;
}

/**
 * Produce a new version of the given simplified instance
 * with all subsumed clauses removed.
 */
inline
SimplifiedInstance remove_subsumed(const SimplifiedInstance& instance) {
    const auto n_all = instance.formula.num_vars();
    auto clause_list = to_clause_list<CVec>(instance.formula);
    eliminate_subsumed(clause_list, n_all);
    return SimplifiedInstance{
        instance.new_to_old,
        ClauseDB(n_all, clause_list), 
        instance.num_concrete
    };
}

}

#endif
