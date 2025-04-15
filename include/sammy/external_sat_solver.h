#ifndef SAMMY_EXTERNAL_SAT_SOLVER_H_INCLUDED_
#define SAMMY_EXTERNAL_SAT_SOLVER_H_INCLUDED_

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <vector>

#include <boost/dll/runtime_symbol_info.hpp>

#include <sammy/cadical_solver.h>
#include <sammy/cmsat5_solver.h>
#include <sammy/kissat_solver.h>
#include <sammy/lingeling_solver.h>

#include "memory_usage.h"
#include "process.h"
#include "range.h"
#include "time.h"

namespace sammy {

enum class ExternalSolverType { KISSAT, CADICAL, LINGELING, CRYPTOMINISAT };

inline const char* get_solver_name(ExternalSolverType e) noexcept {
    switch (e) {
	default:
    case ExternalSolverType::KISSAT:
        return "kissat";
    case ExternalSolverType::CADICAL:
        return "CaDiCaL";
    case ExternalSolverType::LINGELING:
        return "Lingeling";
    case ExternalSolverType::CRYPTOMINISAT:
        return "cryptominisat5";
    }
}

inline std::optional<ExternalSolverType> is_sat_only_call(int argc,
                                                          char** argv) {
    if (argc != 3 || argv[1] != std::string{"--sat-solver-only"}) {
        return std::nullopt;
    }
    std::string solver_name(argv[2]);
    if (solver_name == "kissat") {
        return ExternalSolverType::KISSAT;
    } else if (solver_name == "CaDiCaL") {
        return ExternalSolverType::CADICAL;
    } else if (solver_name == "Lingeling") {
        return ExternalSolverType::LINGELING;
    } else if (solver_name == "cryptominisat5") {
        return ExternalSolverType::CRYPTOMINISAT;
    } else {
        return std::nullopt;
    }
}

inline std::size_t raw_read_stdin(char* buffer, std::size_t size) {
    std::size_t read_bytes = 0;
    while (read_bytes < size) {
        auto rb = ::read(STDIN_FILENO, buffer + read_bytes, size - read_bytes);
        if (rb < 0) {
            throw std::runtime_error("Failed to read from stdin.");
        }
        if (rb == 0) {
            break;
        }
        read_bytes += static_cast<std::size_t>(rb);
    }
    return read_bytes;
}

template <typename Callable> void sat_only_read_input(Callable&& callable) {
    std::vector<std::int32_t> buffer(262144, 0);
    char* buffer_data = reinterpret_cast<char*>(buffer.data());
    std::size_t block_size = buffer.size() * sizeof(std::int32_t);
    bool read_more = true;
    while (read_more) {
        std::size_t read_bytes = raw_read_stdin(buffer_data, block_size);
        if (read_bytes < block_size) {
            read_more = false;
        }
        if (read_bytes % sizeof(std::int32_t) != 0) {
            throw std::runtime_error("Got wrong input size!");
        }
        std::size_t read_lits = read_bytes / sizeof(std::int32_t);
        std::invoke(std::forward<Callable>(callable),
                    IteratorRange(buffer.begin(), buffer.begin() + read_lits));
    }
}

inline std::int32_t sat_only_main_prepare() {
    std::ios::sync_with_stdio(false);
    std::int32_t num_vars;
    if (raw_read_stdin(reinterpret_cast<char*>(&num_vars), sizeof(num_vars)) !=
            sizeof(num_vars) ||
        num_vars < 0)
    {
        std::exit(1);
    }
    return num_vars;
}

inline void sat_only_handle_result(std::optional<bool> result) {
    if (!result) {
        std::exit(1);
    }
    std::size_t mem_peak = current_peak_rss();
    std::cout.write(*result ? "S" : "U", 1);
    std::cout.write(reinterpret_cast<const char*>(&mem_peak), sizeof(mem_peak));
}

template <typename RawMapType>
inline void sat_only_write_raw_map(const RawMapType& model_raw) {
    std::vector<std::uint8_t> model_out;
    std::transform(model_raw.begin(), model_raw.end(),
                   std::back_inserter(model_out),
                   [](bool b) -> std::uint8_t { return b ? 1 : 0; });
    std::cout.write(reinterpret_cast<const char*>(model_out.data()),
                    sizeof(std::uint8_t) * model_out.size());
}

template <typename SatSolver,
          std::enable_if_t<std::is_same_v<SatSolver, CMSAT5Solver>, int> = 0>
int sat_only_main() {
    std::int32_t num_vars = sat_only_main_prepare();
    SatSolver solver;
    solver.new_vars(static_cast<std::size_t>(num_vars));
    sat_only_read_input([&](auto&& lit_range) {
        for (std::int32_t l : lit_range) {
            if (l == 0) {
                solver.finish_clause();
            } else {
                solver.add_literal(solver.lit_from_dimacs_int(l));
            }
        }
    });
    auto result = solver.solve();
    if (!result)
        return 1;
    sat_only_handle_result(result);
    if (*result) {
        auto model_map = solver.get_model();
        sat_only_write_raw_map(model_map.raw());
    }
    return 0;
}

template <
    typename SatSolver,
    std::enable_if_t<std::is_integral_v<typename SatSolver::Lit>, int> = 0>
int sat_only_main() {
    std::int32_t num_vars = sat_only_main_prepare();
    SatSolver solver;
    solver.new_vars(num_vars);
    sat_only_read_input([&](auto&& lit_range) {
        for (std::int32_t l : lit_range) {
            solver.add_literal(l);
        }
    });
    auto result = solver.solve();
    sat_only_handle_result(result);
    if (*result) {
        auto model_map = solver.get_model();
        sat_only_write_raw_map(model_map.raw());
    }
    return 0;
}

inline int sat_only_entry_point(ExternalSolverType solver_type) {
    switch (solver_type) {
	default:
    case ExternalSolverType::KISSAT:
        return sat_only_main<KissatSolver>();
    case ExternalSolverType::CADICAL:
        return sat_only_main<CadicalSolver>();
    case ExternalSolverType::LINGELING:
        return sat_only_main<LingelingSolver>();
    case ExternalSolverType::CRYPTOMINISAT:
        return sat_only_main<CMSAT5Solver>();
    }
}

/**
 * External non-incremental SAT solver
 * running in another process; mostly
 * useful because, unlike other approaches,
 * it can relatively precisely measure its
 * memory requirements.
 */
template <ExternalSolverType EST> class ExternalNonIncrementalSAT {
  public:
    using Var = std::int32_t;
    using Lit = std::int32_t;

    ~ExternalNonIncrementalSAT() {
        std::unique_lock l{m_process_mutex};
        if (m_process.valid()) {
            m_process.reset();
        }
    }

    class ModelMap {
      public:
        bool operator[](Lit l) const noexcept {
            if (l < 0) {
                auto index = -(l + 1);
                assert(index < m_raw.size());
                return !m_raw[index];
            } else {
                auto index = l - 1;
                assert(index < m_raw.size());
                return m_raw[index];
            }
        }

        const std::vector<std::uint8_t>& raw() const noexcept { return m_raw; }

        std::vector<std::uint8_t>& raw() noexcept { return m_raw; }

      private:
        friend class ExternalNonIncrementalSAT;

        explicit ModelMap(std::vector<std::uint8_t> raw) noexcept
            : m_raw(std::move(raw)) {}

        std::vector<std::uint8_t> m_raw;
    };

    Lit new_var() { return ++m_max_var; }

    Lit num_vars() { return m_max_var; }

    void add_literal(Lit l) {
        if (std::abs(l) > m_max_var) {
            m_max_var = std::abs(l);
        }
        if (!m_chunks.empty()) {
            auto& cnk = m_chunks.back();
            if (cnk.size() < cnk.capacity()) {
                cnk.push_back(l);
                return;
            }
        }
        p_add_literal_slow_path(l);
    }

    ModelMap get_model() { return m_model.value(); }

    static const char* name() noexcept { return get_solver_name(EST); }

    template <typename... Lits> void add_literals(Lits&&... lits) {
        (add_literal(lits), ...);
    }

    void finish_clause() { add_literal(0); }

    template <typename... Lits> void add_short_clause(Lits&&... lits) {
        (add_literal(lits), ...);
        finish_clause();
    }

    template <typename Iterator> void add_clause(Iterator begin, Iterator end) {
        std::for_each(begin, end, [this](Lit l) { add_literal(l); });
        finish_clause();
    }

    std::optional<bool> solve() {
        static boost::filesystem::path location{boost::dll::program_location()};
        {
            std::unique_lock l{m_process_mutex};
            if (m_process.valid()) {
                throw std::logic_error("solve() reentered during solve!");
            }
            if (m_dont_start) {
                return std::nullopt;
            }
            m_process =
                Process(location, {"--sat-solver-only", get_solver_name(EST)});
            p_write_formula();
        }
        char result_byte;
        auto do_read = [&](auto* buffer, std::size_t size = 1) {
            char* cbuffer = reinterpret_cast<char*>(buffer);
            std::size_t read_bytes =
                size * sizeof(std::remove_pointer_t<decltype(buffer)>);
            if (m_process.read_from_process(cbuffer, read_bytes) != read_bytes)
            {
                return false;
            }
            return true;
        };
        if (!do_read(&result_byte) || !do_read(&m_solver_total_memory)) {
            // EOF on output indicates this wait won't be long.
            std::unique_lock l{m_process_mutex};
            m_process.wait();
            return std::nullopt;
        }
        /*std::stringstream msg;
        msg << "Formula size: "
            << m_total_size * sizeof(std::int32_t) / (1024.0*1024.0) << " MiB, "
            << "solver used " << m_solver_total_memory / (1024.0*1024.0)
            << " MiB\n";
        std::cerr << msg.str();*/
        if (result_byte == 'U') {
            std::unique_lock l{m_process_mutex};
            m_process.wait();
            return false;
        }
        if (result_byte != 'S') {
            throw std::runtime_error("Unexpected result from solver");
        }
        std::vector<std::uint8_t> model_raw(m_max_var);
        if (!do_read(model_raw.data(), model_raw.size())) {
            std::unique_lock l{m_process_mutex};
            m_process.wait();
            return std::nullopt;
        } else {
            std::unique_lock l{m_process_mutex};
            m_model.emplace(ModelMap(std::move(model_raw)));
            m_process.wait();
            return true;
        }
    }

    void terminate() {
        std::unique_lock l{m_process_mutex};
        if (m_process.valid()) {
            if (!m_process.have_waited()) {
                m_process.terminate();
            }
        } else {
            m_dont_start = true;
        }
    }

  private:
    using Chunk = std::vector<Lit>;

    void p_add_literal_slow_path(Lit l) {
        std::size_t chunk_capacity = (1 << 15);
        if (!m_chunks.empty()) {
            chunk_capacity = 2 * m_chunks.back().capacity();
        }
        m_chunks.emplace_back();
        auto& cnk = m_chunks.back();
        cnk.reserve(chunk_capacity);
        cnk.push_back(l);
    }

    void p_write_formula() {
		m_total_size = 0;
		std::for_each(m_chunks.begin(), m_chunks.end(), 
				      [&] (const Chunk& c) { m_total_size += c.size(); });
        m_process.write_to_process(reinterpret_cast<const char*>(&m_max_var),
                                   sizeof(m_max_var));
        for (const auto& c : m_chunks) {
            m_process.write_to_process(reinterpret_cast<const char*>(c.data()),
                                       c.size() * sizeof(Lit));
        }
        m_chunks.clear();
        m_process.send_eof();
    }

    Var m_max_var{0};
    std::vector<Chunk> m_chunks;
    std::size_t m_total_size{0};
    std::size_t m_solver_total_memory{0};
    std::mutex m_process_mutex;
    Process m_process;
    std::optional<ModelMap> m_model;
    bool m_dont_start{false};
};

} // namespace sammy

#endif
