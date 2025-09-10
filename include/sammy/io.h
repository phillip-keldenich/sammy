#ifndef SAMMY_JSON_IO_H_INCLUDED_
#define SAMMY_JSON_IO_H_INCLUDED_

#include <iostream>
#include <exception>
#include <stdexcept>
#include <filesystem>
#include <fstream>
#include <utility>
#include <algorithm>
#include <string>
#include <tuple>
#include <cstdlib>
#include <boost/exception/diagnostic_information.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/lzma.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <nlohmann/json.hpp>

namespace sammy {

inline nlohmann::json read_json_path(const std::filesystem::path& file) {
    auto ext = file.extension().string();
    std::ifstream raw_input;
    boost::iostreams::filtering_istream input;
    auto input_flags = std::ios::in;
    if (ext == ".xz" || ext == ".XZ" || ext == ".lzma" || ext == ".LZMA") {
        input_flags |= std::ios::binary;
        input.push(boost::iostreams::lzma_decompressor());
    } else if (ext == ".gz" || ext == ".GZ") {
        input_flags |= std::ios::binary;
        input.push(boost::iostreams::gzip_decompressor());
    } else if (ext == ".bz2" || ext == ".BZ2" || ext == ".bzip2" ||
               ext == ".BZIP2")
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
    } catch (const boost::iostreams::lzma_error& e) {
        std::cerr << "Failed to read JSON data from LZMA file " << file << ": "
                  << e.error() << " - " << e.code().message() << std::endl;
        std::exit(1);
    } catch (const boost::exception& e) {
        std::cerr << "Failed to read JSON data from file " << file << ": "
                  << boost::diagnostic_information(e, true) << std::endl;
        std::exit(1);
    } catch (const std::exception& e) {
        std::cerr << "Failed to read JSON data from file " << file << ": "
                  << e.what() << std::endl;
        std::exit(1);
    }
}

inline void write_json_path(const std::filesystem::path& file,
                            const nlohmann::json& d) {
    namespace io = boost::iostreams;
    auto ext = file.extension().string();
    std::ofstream raw_output;
    boost::iostreams::filtering_ostream output;
    auto output_flags = std::ios::out | std::ios::trunc;

    if (ext == ".xz" || ext == ".XZ" || ext == ".lzma" || ext == ".LZMA") {
        output_flags |= std::ios::binary;
        auto level = io::lzma::best_compression;
        output.push(io::lzma_compressor(io::lzma_params{level, 1}));
    } else if (ext == ".gz" || ext == ".GZ" || ext == ".gzip" || ext == ".GZIP")
    {
        output_flags |= std::ios::binary;
        auto level = io::gzip::best_compression;
        output.push(io::gzip_compressor(io::gzip_params{level}));
    } else if (ext == ".bz2" || ext == ".BZ2" || ext == ".bzip2" ||
               ext == ".BZIP2")
    {
        output_flags |= std::ios::binary;
        output.push(io::bzip2_compressor());
    }
    raw_output.exceptions(std::ios::failbit | std::ios::badbit);
    raw_output.open(file, output_flags);
    output.push(raw_output);
    output << d.dump(2);
}

} // namespace sammy

#endif
