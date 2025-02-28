#ifndef SAMMY_IO_H_INCLUDED_
#define SAMMY_IO_H_INCLUDED_

#include <sammy/clause_db.h>
#include <sammy/output.h>
#include <sammy/barrage.h>

namespace sammy {

struct ProblemInput {
    ClauseDB clauses;
    InitialPhaseResult initial_phase;
    nlohmann::json raw_input;
};

inline ProblemInput interpret_problem_json(nlohmann::json input_data) {
    auto num_vars = input_data.at("num_variables").get<std::uint32_t>();
    auto clauses = input_data.at("clauses").get<std::vector<ExternalClause>>();
    return {ClauseDB{num_vars, clauses},
            InitialPhaseResult::import_from_output(input_data),
            std::move(input_data)};
}

inline ProblemInput
load_problem_input(const std::filesystem::path& formula_file) {
    if (!std::filesystem::is_regular_file(formula_file)) {
        throw std::runtime_error("The given formula file does not exist!");
    }
    std::ifstream input;
    input.exceptions(std::ios::failbit | std::ios::badbit);
    input.open(formula_file, std::ios::in);
    nlohmann::json input_data;
    input >> input_data;
    return interpret_problem_json(std::move(input_data));
}

struct SubproblemInput {
    std::vector<Vertex> uncovered_vertices;
    std::vector<Vertex> best_local_mes;
    std::size_t num_nonremoved_configs;
    std::size_t best_global_mes;
    std::size_t best_global_lb;
    std::vector<DynamicBitset> covering_assignments;
    std::optional<DynamicBitset> nonremoved_config;
    nlohmann::json raw_input;
};

inline SubproblemInput interpret_subproblem_json(nlohmann::json input_data) {
    auto uncovered_vertices = lit::internalize(
        input_data.at("uncovered").get<std::vector<ExternalVertex>>());
    auto best_local_mes =
        lit::internalize(input_data.at("initial_uncovered_mes")
                             .get<std::vector<ExternalVertex>>());
    auto num_nonremoved_configs =
        input_data.at("num_nonremoved_configs").get<std::size_t>();
    auto best_global_mes =
        input_data.at("global_best_mes_size").get<std::size_t>();
    auto best_global_lb =
        input_data.at("global_best_lower_bound").get<std::size_t>();
    std::vector<DynamicBitset> covering_assignments;
    for (const auto& e : input_data.at("removed_configs")) {
        covering_assignments.emplace_back(e.get<std::vector<bool>>());
    }
    std::optional<DynamicBitset> remaining{std::nullopt};
    const auto& rem = input_data.at("remaining_config");
    if (!rem.is_null()) {
        remaining.emplace(rem.get<std::vector<bool>>());
    }
    return {std::move(uncovered_vertices),
            std::move(best_local_mes),
            num_nonremoved_configs,
            best_global_mes,
            best_global_lb,
            std::move(covering_assignments),
            std::move(remaining),
            std::move(input_data)};
}

inline SubproblemInput
load_subproblem_input(const std::filesystem::path& subproblem_file) {
    if (!std::filesystem::is_regular_file(subproblem_file)) {
        throw std::runtime_error("The given subproblem file does not exist!");
    }
    std::ifstream input;
    input.exceptions(std::ios::failbit | std::ios::badbit);
    input.open(subproblem_file, std::ios::in);
    nlohmann::json input_data;
    input >> input_data;
    return interpret_subproblem_json(std::move(input_data));
}

}

#endif
