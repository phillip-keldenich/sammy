#ifndef SAMMY_OUTPUT_H_INCLUDED_
#define SAMMY_OUTPUT_H_INCLUDED_

#include "dynamic_bitset.h"
#include "literals.h"
#include "simplification_stats.h"
#include "simplify_datastructure.h"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <numeric>
#include <stdexcept>

#include <nlohmann/json.hpp>

#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>

namespace sammy {

using OutputObject = nlohmann::json;

inline std::size_t estimate_bytes_used(const OutputObject& obj) noexcept {
    if (obj.is_string()) {
        return sizeof(obj) + sizeof(nlohmann::json::string_t) +
               obj.get_ref<const nlohmann::json::string_t&>().capacity();
    }
    if (obj.is_binary()) {
        return sizeof(obj) + sizeof(nlohmann::json::binary_t) +
               obj.get_ref<const nlohmann::json::binary_t&>().capacity();
    }
    if (obj.is_array()) {
        const auto& ref = obj.get_ref<const nlohmann::json::array_t&>();
        std::size_t initial = sizeof(nlohmann::json::array_t) +
                              sizeof(obj) * (1 + ref.capacity() - ref.size());
        return std::accumulate(ref.begin(), ref.end(), initial,
                               [](std::size_t acc, const OutputObject& o) {
                                   return acc + estimate_bytes_used(o);
                               });
    }
    if (obj.is_object()) {
        const auto& ref = obj.get_ref<const nlohmann::json::object_t&>();
        std::size_t initial = sizeof(obj) + sizeof(nlohmann::json::object_t) +
                              (3 * sizeof(void*) +
                               sizeof(nlohmann::json::string_t) + sizeof(obj)) *
                                  ref.size();
        return std::accumulate(ref.begin(), ref.end(), initial,
                               [](std::size_t acc, const auto& p) {
                                   return acc + estimate_bytes_used(p.second) +
                                          p.first.capacity();
                               });
    }
    return sizeof(obj);
}

inline OutputObject bitset_to_json(const DynamicBitset& b) {
    OutputObject result;
    for (std::size_t i = 0; i < b.size(); ++i) {
        result.push_back(bool(b[i]));
    }
    return result;
}

inline OutputObject bitsets_to_json(const std::vector<DynamicBitset>& b) {
    OutputObject result;
    for (const DynamicBitset& bs : b) {
        result.emplace_back(bitset_to_json(bs));
    }
    return result;
}

struct StoredEvent {
    std::string type;
    double time;
    OutputObject data;
};

class EventRecorder {
  public:
    explicit EventRecorder() : m_begin_time(std::chrono::steady_clock::now()) {}

    double now() const noexcept {
        auto t = std::chrono::steady_clock::now();
        return std::chrono::duration_cast<Seconds>(t - m_begin_time).count();
    }

    void store_event_silent(std::string type) {
        m_events.emplace_back(
            StoredEvent{std::move(type), now(), OutputObject{}});
    }

    void store_event(std::string type) {
        m_events.emplace_back(
            StoredEvent{std::move(type), now(), OutputObject{}});
        if (print_events)
            p_print();
    }

    template <typename... PrintArgs>
    void store_event(std::string type, OutputObject data,
                     PrintArgs&&... p_args) {
        if (!data.is_null() && !data.is_object()) {
            throw std::runtime_error("Event data must be an object or null!");
        }
        m_events.emplace_back(
            StoredEvent{std::move(type), now(), std::move(data)});
        if (print_events)
            p_print(std::forward<PrintArgs>(p_args)...);
    }

    void store_event_silent(std::string type, OutputObject data) {
        m_events.emplace_back(
            StoredEvent{std::move(type), now(), std::move(data)});
    }

    const std::vector<StoredEvent>& events() const noexcept { return m_events; }

    void set_print_events(bool pe) noexcept { print_events = pe; }

    void merge_events(const EventRecorder& other,
                      const std::string& worker_id) {
        double delta_t = seconds_between(m_begin_time, other.m_begin_time);
        std::size_t old_size = m_events.size();
        for (const StoredEvent& event : other.m_events) {
            m_events.emplace_back(event);
            m_events.back().time += delta_t;
            try {
                m_events.back().data["worker_id"] = worker_id;
            } catch (...) {
                std::cerr << "Failed to add worker_id to event!" << std::endl;
                std::cerr << "EVENT: " << m_events.back().type << std::endl;
                throw;
            }
        }
        std::inplace_merge(m_events.begin(), m_events.begin() + old_size,
                           m_events.end(),
                           [](const StoredEvent& a, const StoredEvent& b) {
                               return a.time < b.time;
                           });
    }

    /**
     * Synchronize the time stamp of this recorder with another,
     * pretending that *this was started at the same time as other.
     */
    void synchronize_with(const EventRecorder& other) {
        m_begin_time = other.m_begin_time;
    }

    /**
     * Do some estimation of the number of bytes used
     * by this event recorder.
     */
    std::size_t bytes_used() const noexcept {
        std::size_t result = sizeof(*this);
        result +=
            (sizeof(StoredEvent) - sizeof(OutputObject)) * m_events.capacity();
        return std::accumulate(m_events.begin(), m_events.end(), result,
                               [](std::size_t acc, const StoredEvent& e) {
                                   return acc + estimate_bytes_used(e.data);
                               });
    }

  private:
    void p_print_args(std::ostream& outbuffer) { outbuffer << " }"; }

    template <typename PArg1, typename... PArgs>
    void p_print_args(std::ostream& outbuffer, PArg1&& parg1,
                      PArgs&&... pargs) {
        outbuffer << std::forward<PArg1>(parg1) << ": "
                  << m_events.back().data[std::forward<PArg1>(parg1)];
        if (sizeof...(pargs)) {
            outbuffer << ", ";
        }
        p_print_args(outbuffer, std::forward<PArgs>(pargs)...);
    }

    template <typename... PrintArgs> void p_print(PrintArgs&&... p_args) {
        std::ostringstream outbuffer;
        char buffer[32];
        auto& event = m_events.back();
        std::snprintf(buffer, 32, "[%11.4f] ", event.time);
        outbuffer << buffer << event.type;
        if (sizeof...(PrintArgs)) {
            outbuffer << " { ";
            p_print_args(outbuffer, std::forward<PrintArgs>(p_args)...);
        }
        outbuffer << std::endl;
        std::cout << outbuffer.str();
    }

    using Seconds = std::chrono::duration<double>;
    std::chrono::steady_clock::time_point m_begin_time;
    std::vector<StoredEvent> m_events;
    bool print_events = false;
};

inline void export_events(OutputObject& output,
                          const std::vector<StoredEvent>& events) {
    nlohmann::json& out = output["events"];
    for (const StoredEvent& event : events) {
        out.push_back(event.data);
        nlohmann::json& current = out.back();
        current["time"] = event.time;
        current["type"] = event.type;
    }
}

inline void export_events(OutputObject& output, const EventRecorder& recorder) {
    export_events(output, recorder.events());
}

inline void export_simplification(
    OutputObject& output, const SimplifyDatastructure& simplifier,
    const SimplifiedInstance& compressed,
    const SimplificationStats* simplification_stats = nullptr) {
    nlohmann::json& json = output["simplification"];
    json["simplified_formula"] = compressed.formula.export_all_clauses();
    std::vector<ExternalLit> new_to_old_external;
    std::transform(compressed.new_to_old.begin(), compressed.new_to_old.end(),
                   std::back_inserter(new_to_old_external), [](Var o) {
                       return lit::externalize(lit::positive_lit(o));
                   });
    std::vector<ExternalClause> reconstruction_stack_external;
    std::transform(simplifier.get_reconstruction_stack().begin(),
                   simplifier.get_reconstruction_stack().end(),
                   std::back_inserter(reconstruction_stack_external),
                   [](const auto& c) { return lit::externalize(c); });
    json["compressed_variable_to_original"] = new_to_old_external;
    json["reconstruction_stack"] = reconstruction_stack_external;
    if (simplification_stats) {
        add_simplification_stats(json, *simplification_stats);
    }
}

inline void export_solution(OutputObject& output,
                            const std::vector<std::vector<bool>>& solution,
                            const std::string& tag,
                            const OutputObject& extra_info) {
    nlohmann::json& json = output[tag];
    json["solution"] = extra_info;
    json["solution"]["configurations"] = solution;
}

inline void export_bound(OutputObject& output, const std::vector<Vertex>& lb,
                         const std::string& tag,
                         const OutputObject& extra_info) {
    std::vector<std::pair<ExternalLit, ExternalLit>> external;
    external.reserve(lb.size());
    auto externalize_vertex = [](Vertex v) {
        return std::pair<ExternalLit, ExternalLit>{lit::externalize(v.first),
                                                   lit::externalize(v.second)};
    };
    std::transform(lb.begin(), lb.end(), std::back_inserter(external),
                   externalize_vertex);
    nlohmann::json& json = output[tag];
    json["lb"] = extra_info;
    json["lb"]["mutually_exclusive_set"] = lb;
}

inline std::string tolower(const std::string& mixed) {
    std::string result;
    result.reserve(mixed.size());
    std::transform(mixed.begin(), mixed.end(), std::back_inserter(result),
                   [](char c) { return char(std::tolower(c)); });
    return result;
}

inline bool recognized_output_extension(const std::filesystem::path& path) {
    auto lowext = tolower(path.extension().string());
    return lowext == ".json" || lowext == ".jsonl";
}

inline void output_data(const OutputObject& output,
                        const std::filesystem::path& path) {
    auto lowext = tolower(path.extension().string());
    if (lowext == ".json") {
        std::ofstream output_stream;
        output_stream.exceptions(std::ios::failbit | std::ios::badbit);
        output_stream.open(path, std::ios::out | std::ios::trunc);
        output_stream << output;
    } else if (lowext == ".jsonl") {
        std::ofstream output_stream;
        output_stream.exceptions(std::ios::failbit | std::ios::badbit);
        output_stream.open(path, std::ios::out | std::ios::app);
        output_stream << output << std::endl;
    } else {
        throw std::runtime_error("Unrecognized output extension '" + lowext +
                                 "'!");
    }
}

} // namespace sammy

#endif
