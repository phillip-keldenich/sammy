#ifndef SAMMY_IO_H_INCLUDED_
#define SAMMY_IO_H_INCLUDED_

#include <sammy/clause_db.h>
#include <sammy/output.h>
#include <sammy/barrage.h>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/exception/diagnostic_information.hpp>

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

inline nlohmann::json 
read_json_path(const std::filesystem::path& file)
{
    auto ext = file.extension().string();
    std::ifstream raw_input;
    boost::iostreams::filtering_istream input;
    auto input_flags = std::ios::in;
    if(ext == ".xz" || ext == ".XZ" ||
       ext == ".lzma" || ext == ".LZMA")
    {
        input_flags |= std::ios::binary;
        input.push(boost::iostreams::lzma_decompressor());
    } else if(ext == ".gz" || ext == ".GZ") {
        input_flags |= std::ios::binary;
        input.push(boost::iostreams::gzip_decompressor());
    } else if(ext == ".bz2" || ext == ".BZ2" || 
              ext == ".bzip2" || ext == ".BZIP2") 
    {
        input_flags |= std::ios::binary;
        input.push(boost::iostreams::bzip2_decompressor());
    }
    try {
        raw_input.exceptions(std::ios::failbit | std::ios::badbit);
        raw_input.open(file, input_flags);
        raw_input.exceptions(std::ios::badbit);
    } catch (const std::exception& e) {
        throw std::runtime_error("Failed to open the input file: " + 
                                 file.string());
    }
    try {
        input.push(raw_input);
        return nlohmann::json::parse(input);
    } catch(const boost::iostreams::lzma_error& e) {
        std::cerr << "Failed to read JSON data from LZMA file "
                  << file << ": " << e.error() << " - " 
                  << e.code().message() << std::endl;
        std::exit(1);
    } catch(const boost::exception& e) {
        std::cerr << "Failed to read JSON data from file "
                  << file << ": " << boost::diagnostic_information(e, true)
                  << std::endl;
        std::exit(1);
    } catch(const std::exception& e) {
        std::cerr << "Failed to read JSON data from file "
                  << file << ": " << e.what() << std::endl;
        std::exit(1);
    }
}

inline void
write_json_path(const std::filesystem::path& file, const nlohmann::json& d) {
    namespace io = boost::iostreams;
    auto ext = file.extension().string();
    std::ofstream raw_output;
    boost::iostreams::filtering_ostream output;
    auto output_flags = std::ios::out | std::ios::trunc;

    if(ext == ".xz" || ext == ".XZ" ||
       ext == ".lzma" || ext == ".LZMA")
    {
        output_flags |= std::ios::binary;
        auto level = io::lzma::best_compression;
        output.push(io::lzma_compressor(io::lzma_params{level, 1}));
    } else if(ext == ".gz" || ext == ".GZ" ||
              ext == ".gzip" || ext == ".GZIP")
    {
        output_flags |= std::ios::binary;
        auto level = io::gzip::best_compression;
        output.push(io::gzip_compressor(io::gzip_params{level}));
    } else if(ext == ".bz2" || ext == ".BZ2" ||
              ext == ".bzip2" || ext == ".BZIP2")
    {
        output_flags |= std::ios::binary;
        output.push(io::bzip2_compressor());
    }
    raw_output.exceptions(std::ios::failbit | std::ios::badbit);
    raw_output.open(file, output_flags);
    output.push(raw_output);
    output << d.dump(2);
}

inline ProblemInput
load_problem_input(const std::filesystem::path& formula_file) {
    if (!std::filesystem::is_regular_file(formula_file)) {
        throw std::runtime_error("The given formula file does not exist!");
    }
    return interpret_problem_json(read_json_path(formula_file));
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
    return interpret_subproblem_json(read_json_path(subproblem_file));
}

}

#endif
