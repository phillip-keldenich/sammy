#ifndef SAMMY_INITIAL_PHASE_RESULT_H_INCLUDED_
#define SAMMY_INITIAL_PHASE_RESULT_H_INCLUDED_

#include "literals.h"
#include "output.h"
#include "pair_infeasibility_map.h"
#include <sstream>

namespace sammy {

/**
 * Store the result of the initial phase, i.e.,
 * the phase that runs the heuristic to find an
 * initial solution and an initial mutually
 * exclusive set.
 */
struct InitialPhaseResult {
    PairInfeasibilityMap inf_map;
    std::vector<std::vector<bool>> best_solution;
    std::vector<Vertex> best_spawners;
    std::vector<Vertex> best_mutually_exclusive;
    std::vector<Vertex> all_spawners;
    std::vector<Vertex> coloring_order;
    std::size_t universe_size;

    static InitialPhaseResult import_from_output(const OutputObject& obj) {
        auto internalize = [](const auto& obj, const char* key) {
            return lit::internalize(
                obj.at(key).template get<std::vector<ExternalVertex>>());
        };
        return InitialPhaseResult{
            PairInfeasibilityMap::import_bits(obj.at("infeasibility_map")),
            obj.at("best_solution").get<std::vector<std::vector<bool>>>(),
            internalize(obj, "best_spawners"),
            internalize(obj, "best_mutually_exclusive"),
            internalize(obj, "all_spawners"),
            internalize(obj, "coloring_order"),
            obj.at("universe_size").get<std::size_t>()};
    }

    OutputObject export_to_output() const {
        return OutputObject{
            {"infeasibility_map", inf_map.export_bits()},
            {"best_solution", best_solution},
            {"best_spawners", lit::externalize(best_spawners)},
            {"best_mutually_exclusive",
             lit::externalize(best_mutually_exclusive)},
            {"all_spawners", lit::externalize(all_spawners)},
            {"coloring_order", lit::externalize(coloring_order)},
            {"universe_size", universe_size}};
    }
};

/**
 * Generate an OutputObject from a set of clauses
 * (processed or not) and an initial phase result.
 */
inline OutputObject
export_initial_phase_result(const ClauseDB& clauses,
                            const InitialPhaseResult& result) {
    std::stringstream clause_export_stream;
    clauses.export_to(clause_export_stream);
    return OutputObject{{"clauses", clause_export_stream.str()},
                        {"initial_phase_result", result.export_to_output()}};
}

inline std::pair<ClauseDB, InitialPhaseResult>
import_initial_phase_result(const OutputObject& obj) {
    std::stringstream clause_import_stream(
        obj.at("clauses").get<std::string>());
    return {
        ClauseDB::import_from(clause_import_stream),
        InitialPhaseResult::import_from_output(obj.at("initial_phase_result"))};
}

} // namespace sammy

#endif
