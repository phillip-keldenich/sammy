#ifndef SAMMY_SIMPLIFICATION_STATS_H_INCLUDED_
#define SAMMY_SIMPLIFICATION_STATS_H_INCLUDED_

#include "clause_db.h"

#include <cstddef>
#include <iostream>
#include <string>
#include <unordered_map>

namespace sammy {

struct SimplificationStats {
    double simplification_time = 0.0;
    std::size_t simplification_rounds = 0;
    std::size_t variables_eliminated_by_equality = 0;
    std::size_t variables_fixed = 0;
    std::size_t variables_eliminated_by_resolution = 0;
    std::size_t pure_literals_eliminated = 0;
    std::size_t failed_literals = 0;
    std::size_t clauses_strengthened_by_binary_resolution = 0;
    std::size_t clauses_strengthened_by_vivification = 0;
    std::size_t clauses_subsumed = 0;
    std::size_t conflict_clauses_learned = 0;
    std::size_t binary_contrapositivities_learned = 0;
    Var n_all_before = 0;
    Var n_all_after = 0;
    Var n_concrete_before = 0;
    Var n_concrete_after = 0;
    CRef n_clause_before = 0;
    CRef n_clause_after = 0;
    CRef n_binary_before = 0;
    CRef n_binary_after = 0;
    CRef n_long_before = 0;
    CRef n_long_after = 0;
    CRef total_clause_size_before = 0;
    CRef total_clause_size_after = 0;

    void capture_before(const ClauseDB& input, Var n_concrete) noexcept {
        n_all_before = input.num_vars();
        n_concrete_before = n_concrete;
        n_clause_before = input.num_clauses();
        n_binary_before = input.num_binaries();
        n_long_before = n_clause_before - n_binary_before - input.num_unaries();
        total_clause_size_before = 2 * n_binary_before + input.num_unaries() +
                                   input.literal_db_size() - n_long_before;
    }

    void capture_after(const ClauseDB& simplified, Var n_concrete) noexcept {
        n_all_after = simplified.num_vars();
        n_concrete_after = n_concrete;
        n_clause_after = simplified.num_clauses();
        n_binary_after = simplified.num_binaries();
        n_long_after =
            n_clause_after - n_binary_after - simplified.num_unaries();
        total_clause_size_after = 2 * n_binary_after +
                                  simplified.num_unaries() +
                                  simplified.literal_db_size() - n_long_after;
    }
};

template <typename JSONType>
inline void add_simplification_stats(JSONType& output,
                                     const SimplificationStats& stats) {
    std::unordered_map<std::string, std::size_t> stat_section;
    stat_section["variables_eliminated_by_equality"] =
        stats.variables_eliminated_by_equality;
    stat_section["variables_fixed"] = stats.variables_fixed;
    stat_section["failed_literals"] = stats.failed_literals;
    stat_section["variables_eliminated_by_resolution"] =
        stats.variables_eliminated_by_resolution;
    stat_section["pure_literals_eliminated"] = stats.pure_literals_eliminated;
    stat_section["clauses_strengthened_by_binary_resolution"] =
        stats.clauses_strengthened_by_binary_resolution;
    stat_section["clauses_strengthened_by_vivification"] =
        stats.clauses_strengthened_by_vivification;
    stat_section["clauses_subsumed"] = stats.clauses_subsumed;
    stat_section["conflict_clauses_learned"] = stats.conflict_clauses_learned;
    stat_section["binary_contrapositivities_learned"] =
        stats.binary_contrapositivities_learned;
    stat_section["simplification_rounds"] = stats.simplification_rounds;
    stat_section["variables_before"] = stats.n_all_before;
    stat_section["variables_after"] = stats.n_all_after;
    stat_section["concrete_before"] = stats.n_concrete_before;
    stat_section["concrete_after"] = stats.n_concrete_after;
    stat_section["clauses_before"] = stats.n_clause_before;
    stat_section["clauses_after"] = stats.n_clause_after;
    stat_section["binaries_before"] = stats.n_binary_before;
    stat_section["binaries_after"] = stats.n_binary_after;
    stat_section["long_clauses_before"] = stats.n_long_before;
    stat_section["long_clauses_after"] = stats.n_long_after;
    stat_section["formula_length_before"] = stats.total_clause_size_before;
    stat_section["formula_length_after"] = stats.total_clause_size_after;
    output["simplification_stats"] = stat_section;
    output["simplification_stats"]["simplification_time"] =
        stats.simplification_time;
}

} // namespace sammy

#endif
