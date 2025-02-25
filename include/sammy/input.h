#ifndef SAMMY_INPUT_H_INCLUDED_
#define SAMMY_INPUT_H_INCLUDED_

#include <nlohmann/json.hpp>
#include "clause_db.h"
#include <filesystem>
#include <iostream>
#include <fstream>
#include <string>

namespace sammy {

struct InputData {
    ClauseDB formula;
    Var num_concrete;
    std::string name;
};

inline InputData read_input(const std::filesystem::path& path) {
    std::ifstream input;
    input.exceptions(std::ios::badbit | std::ios::failbit);
    input.open(path, std::ios::in);
    nlohmann::json json_data;
    input >> json_data;
    if(json_data.at("type") != "software configuration model") {
        throw std::runtime_error("JSON data missing the 'software configuration model' type flag!");
    }
    auto clauses = json_data.at("cnf_clauses").get<std::vector<ExternalClause>>();
    auto n_all = json_data.at("num_variables").get<Var>();
    auto n_concrete = json_data.at("num_concrete_features").get<Var>();
    return InputData{ClauseDB{n_all, clauses}, n_concrete, json_data.at("name").get<std::string>()};
}

}

#endif
