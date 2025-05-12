==> ./equality_graph.h <==
#ifndef SAMMY_EQUALITY_GRAPH_H_INCLUDED_
#define SAMMY_EQUALITY_GRAPH_H_INCLUDED_

#include "error.h"
#include "literals.h"
#include <algorithm>
#include <utility>

namespace sammy {

/**
 * Basically a union-find datastructure
 * that allows making literals equal to other literals,
 * or to a true 'dummy' variable's literals
 * (i.e., fixing it to true or false).
 */
class EqualityGraph {
    std::vector<Var> m_component_root;
    std::vector<bool> m_component_root_negated;
    std::vector<std::uint8_t> m_rank;
    std::vector<std::pair<Var, bool>> m_path_buffer;
    Var m_dummy; // we have a dummy variable 'True' at num_vars + 1

  public:
    explicit EqualityGraph(Var num_vars)
        : m_component_root(num_vars + 1, NIL),
          m_component_root_negated(num_vars + 1, false),
          m_rank(num_vars + 1, std::uint8_t(0)), m_dummy(num_vars) {
        m_path_buffer.reserve(64);
        // make sure the dummy will always be picked as representative
        m_rank[m_dummy] = std::numeric_limits<std::uint8_t>::max();
        for (Var v = 0; v <= num_vars; ++v) {
            m_component_root[v] = v;
        }
    }

    Var dummy() const noexcept { return m_dummy; }

    Lit find(Lit l) noexcept {
        Var v = lit::var(l);
        bool negated = lit::negative(l);
        Var r;
        bool rn;
        std::tie(r, rn) = p_find_pcompress(v);
        return (negated ^ rn) ? lit::negative_lit(r) : lit::positive_lit(r);
    }

    bool make_equal(Lit x, Lit y) {
        Var v1 = lit::var(x);
        Var v2 = lit::var(y);
        auto [r1, r1n] = p_find_pcompress(v1);
        auto [r2, r2n] = p_find_pcompress(v2);
        r1n ^= lit::negative(x);
        r2n ^= lit::negative(y);
        bool negated = r1n ^ r2n;
        if (r1 == r2) {
            if (negated)
                throw UNSATError();
            return false;
        }
        auto k1 = m_rank[v1];
        auto k2 = m_rank[v2];
        if (k1 < k2) {
            m_component_root[r1] = r2;
            m_component_root_negated[r1] = negated;
        } else {
            if (k1 == k2) {
                m_rank[r1] += 1;
            }
            m_component_root[r2] = r1;
            m_component_root_negated[r2] = negated;
        }
        return true;
    }

    bool make_true(Lit l) { return make_equal(l, lit::positive_lit(dummy())); }

    bool make_false(Lit l) { return make_equal(l, lit::negative_lit(dummy())); }

    std::vector<Lit> compute_old_to_new_map() {
        Var nv = m_component_root.size() - 1;
        std::vector<Lit> old_to_new(2 * nv, NIL);
        for (Var v = 0; v < nv; ++v) {
            Lit p = lit::positive_lit(v);
            Lit m = find(p);
            if (lit::var(m) == dummy()) {
                if (lit::negative(m)) {
                    old_to_new[p] = simplify::fixed_negative();
                    old_to_new[lit::negate(p)] = simplify::fixed_positive();
                } else {
                    old_to_new[p] = simplify::fixed_positive();
                    old_to_new[lit::negate(p)] = simplify::fixed_negative();
                }
            } else {
                old_to_new[p] = m;
                old_to_new[lit::negate(p)] = lit::negate(m);
            }
        }
        return old_to_new;
    }

  private:
    std::pair<Var, bool> p_find_pcompress(Var v) {
        Var p;
        bool pneg, rneg = false;
        for (;;) {
            p = m_component_root[v];
            pneg = m_component_root_negated[v];
            rneg ^= pneg;
            if (p == v) {
                bool cneg = false;
                for (auto i = m_path_buffer.rbegin(), e = m_path_buffer.rend();
                     i != e; ++i)
                {
                    auto [v_i, n_i] = *i;
                    cneg ^= n_i;
                    m_component_root[v_i] = p;
                    m_component_root_negated[v_i] = cneg;
                }
                m_path_buffer.clear();
                return {p, rneg};
            }
            m_path_buffer.emplace_back(v, pneg);
            v = p;
        }
    }
};

} // namespace sammy

#endif
==> ./vertex_operations.h <==
#ifndef SAMMY_VERTEX_OPERATIONS_H_INCLUDED_
#define SAMMY_VERTEX_OPERATIONS_H_INCLUDED_

#include "literals.h"
#include "shared_db_propagator.h"

namespace sammy {

/**
 * @brief Reset a propagator and push a vertex to it.
 *        Does not use conflict resolution and
 *        throws an exception if the vertex cannot
 *        be pushed without conflict.
 *
 * @param propagator
 * @param vertex
 */
inline void reset_and_push_noresolve(SharedDBPropagator& propagator,
                                     Vertex vertex) {
    propagator.reset_or_throw();
    Lit lmin = vertex.first;
    Lit lmax = vertex.second;
    if (propagator.is_false(lmin) || propagator.is_false(lmax)) {
        throw std::logic_error("Infeasible vertex pushed!");
    }
    if (propagator.is_open(lmin)) {
        if (!propagator.push_level(lmin) || propagator.is_false(lmax)) {
            throw std::logic_error("Infeasible vertex pushed!");
        }
    }
    if (propagator.is_open(lmax)) {
        if (!propagator.push_level(lmax)) {
            throw std::logic_error("Infeasible vertex pushed!");
        }
    }
}

/**
 * @brief Tests whether the given vertex can be pushed onto the given
 * propagator. Does not learn from conflicts; returns the propagator into the
 * previous state.
 */
inline bool can_push(SharedDBPropagator& propagator, Vertex vertex) {
    Lit lmin = vertex.first;
    Lit lmax = vertex.second;
    if (propagator.is_false(lmin) || propagator.is_false(lmax)) {
        return false;
    }
    int push_count = 0;
    if (propagator.is_open(lmin)) {
        if (!propagator.push_level(lmin)) {
            propagator.pop_level();
            return false;
        }
        if (propagator.is_false(lmax)) {
            propagator.pop_level();
            return false;
        }
        ++push_count;
    }
    bool result = true;
    if (propagator.is_open(lmax)) {
        result = propagator.push_level(lmax);
        ++push_count;
    }
    for (int i = 0; i < push_count; ++i) {
        propagator.pop_level();
    }
    return result;
}

inline int push_vertex(SharedDBPropagator& propagator, Vertex vertex) {
    Lit lmin = vertex.first;
    Lit lmax = vertex.second;
    if (propagator.is_false(lmin) || propagator.is_false(lmax)) {
        return -1;
    }
    int pushed = 0;
    if (propagator.is_open(lmin)) {
        if (!propagator.push_level(lmin) || propagator.is_false(lmax)) {
            propagator.pop_level();
            return -1;
        }
        pushed = 1;
    }
    if (propagator.is_open(lmax)) {
        if (!propagator.push_level(lmax)) {
            propagator.pop_level();
            if (pushed)
                propagator.pop_level();
            return -1;
        }
        ++pushed;
    }
    return pushed;
}

inline bool push_vertex_pair(SharedDBPropagator& propagator, Vertex v1,
                             Vertex v2) {
    int pc1 = push_vertex(propagator, v1);
    if (pc1 < 0)
        return false;
    if (push_vertex(propagator, v2) < 0) {
        for (int i = 0; i < pc1; ++i)
            propagator.pop_level();
        return false;
    }
    return true;
}

} // namespace sammy

#endif
==> ./error.h <==
#ifndef SAMMY_ERROR_H_INCLUDED_
#define SAMMY_ERROR_H_INCLUDED_

#include <exception>
#include <stdexcept>
#include <string>

namespace sammy {

class UNSATError : public std::runtime_error {
  public:
    UNSATError() : std::runtime_error("The given model is unsatisfiable") {}
};

} // namespace sammy

#endif
==> ./verify.h <==
#ifndef SAMMY_VERIFY_H_INCLUDED_
#define SAMMY_VERIFY_H_INCLUDED_

#include "clause_db.h"
#include "greedysat.h"
#include "logging.h"
#include "range.h"
#include "shared_db_propagator.h"
#include "vertex_operations.h"
#include <optional>
#include <sstream>
#include <string>

namespace sammy {

template <typename AssignmentType>
inline std::optional<std::string>
solution_has_error(const ClauseDB& clauses, const AssignmentType& solution,
                   Reason* reason = nullptr) {
    if (solution.size() != clauses.num_vars()) {
        return "Error: Solution length (= " + std::to_string(solution.size()) +
               ") does not match number of variables (= " +
               std::to_string(clauses.num_vars()) + ")!";
    }
    auto literal_satisfied = [&](Lit l) -> bool {
        Var v = lit::var(l);
        bool val = solution[v];
        return val ^ lit::negative(l);
    };
    auto clause_safisfied = [&](const Lit* beg, const Lit* end) {
        return std::any_of(beg, end, literal_satisfied);
    };
    auto vio_unary = std::find_if(clauses.unary_literals().begin(),
                                  clauses.unary_literals().end(),
                                  [&](Lit l) { return !literal_satisfied(l); });
    if (vio_unary != clauses.unary_literals().end()) {
        if (reason) {
            *reason = Reason::Unary{*vio_unary};
        }
        return "Error: Unary clause (" +
               std::to_string(lit::externalize(*vio_unary)) + ") is violated!";
    }
    auto vio_binary =
        std::find_if(clauses.binary_clauses().begin(),
                     clauses.binary_clauses().end(), [&](auto cl) {
                         return !literal_satisfied(cl.first) &&
                                !literal_satisfied(cl.second);
                     });
    if (vio_binary != clauses.binary_clauses().end()) {
        if (reason) {
            *reason = Reason::Binary{vio_binary->first, vio_binary->second};
        }
        return "Error: Binary clause (" +
               std::to_string(lit::externalize(vio_binary->first)) + " " +
               std::to_string(lit::externalize(vio_binary->second)) +
               ") is violated!";
    }
    for (CRef c = 1, n = clauses.literal_db_size(); c < n;
         c = clauses.next_clause(c))
    {
        auto lits = clauses.lits_of(c);
        if (!clause_safisfied(lits.begin(), lits.end())) {
            if (reason) {
                *reason = Reason::Clause{clauses.clause_length(c), c};
            }
            std::ostringstream buffer;
            buffer << "Error: Clause (";
            print_clause_external(buffer, lits);
            buffer << ") is violated!";
            return buffer.str();
        }
    }
    return std::nullopt;
}

inline std::string vertex_to_string(Vertex v) {
    return "(" + std::to_string(lit::externalize(v.first)) + ", " +
           std::to_string(lit::externalize(v.second)) + ")";
}

inline std::optional<std::string>
mutually_exclusive_set_has_error(ClauseDB& clauses, Var n_concrete,
                                 const std::vector<Vertex>& vertices_) {
    // normalize vertices and look for invalid ones
    auto normalize = [](Vertex v) -> Vertex {
        return {(std::min)(v.first, v.second), (std::max)(v.first, v.second)};
    };
    auto is_bad_vertex = [](Vertex v) {
        return lit::var(v.first) == lit::var(v.second);
    };
    auto is_oor_vertex = [&](Vertex v) {
        return lit::var(v.second) >= n_concrete;
    };
    std::vector<Vertex> vertices;
    vertices.reserve(vertices_.size());
    std::transform(vertices_.begin(), vertices_.end(),
                   std::back_inserter(vertices), normalize);
    auto bad_vertex =
        std::find_if(vertices.begin(), vertices.end(), is_bad_vertex);
    if (bad_vertex != vertices.end()) {
        return "Error: Vertex " + vertex_to_string(*bad_vertex) +
               " in mutually exclusive set!";
    }
    auto oor_vertex =
        std::find_if(vertices.begin(), vertices.end(), is_oor_vertex);
    if (oor_vertex != vertices.end()) {
        return "Error: Vertex " + vertex_to_string(*oor_vertex) +
               " (n_concrete = " + std::to_string(n_concrete) +
               ") in mutually exclusive set!";
    }

    // check for duplicates
    EdgeSet present_vertices;
    for (Vertex v : vertices) {
        if (!present_vertices.emplace(v).second) {
            return "Error: Duplicate vertex " + vertex_to_string(v) +
                   " in mutually exclusive set!";
        }
    }

    // verify mutual exclusiveness
    SharedDBPropagator propagator{&clauses};
    std::vector<std::pair<Vertex, Vertex>> explicitly_check_vertices;
    for (std::size_t ind_in_vertices = 0, n_vertices = vertices.size();
         ind_in_vertices < n_vertices; ++ind_in_vertices)
    {
        Vertex v = vertices[ind_in_vertices];
        reset_and_push_noresolve(propagator, v);
        for (std::size_t j = ind_in_vertices + 1; j < n_vertices; ++j) {
            Vertex w = vertices[j];
            if (propagator.is_false(w.first) || propagator.is_false(w.second)) {
                continue;
            }
            bool pushed_w1 = false;
            if (propagator.is_open(w.first)) {
                if (!propagator.push_level(w.first) ||
                    propagator.is_false(w.second))
                {
                    propagator.pop_level();
                    continue;
                }
                pushed_w1 = true;
            }
            // propagator has w.first on trail and w.second isn't false
            if (propagator.is_open(w.second)) {
                bool pres = propagator.push_level(w.second);
                propagator.pop_level();
                if (pushed_w1)
                    propagator.pop_level();
                if (!pres) {
                    continue;
                }
            } else if (pushed_w1) {
                propagator.pop_level();
            }
            // now we're back at the level we had after pushing v
            explicitly_check_vertices.emplace_back(v, w);
        }
    }
    for (auto [v, w] : explicitly_check_vertices) {
        Lit literals[4] = {v.first, v.second, w.first, w.second};
        ClauseDB new_clauses{clauses};
        new_clauses.add_clause(&literals[0], &literals[1]);
        new_clauses.add_clause(&literals[1], &literals[2]);
        new_clauses.add_clause(&literals[2], &literals[3]);
        new_clauses.add_clause(&literals[3], literals + 4);
        SharedDBPropagator pp{&new_clauses};
        GreedySAT<PreferFalse> gsolver{&pp, n_concrete, PreferFalse{}};
        try {
            if (gsolver.solve()) {
                return "Error: Vertices " + vertex_to_string(v) + " and " +
                       vertex_to_string(w) + " are not mutually exclusive!\n";
            }
        } catch (UNSATError&) {
            continue;
        }
    }

    return std::nullopt;
}

inline bool verify_solution(const ClauseDB& clauses,
                            const std::vector<bool>& solution) {
    return !solution_has_error(clauses, solution);
}

} // namespace sammy

#endif
==> ./experiment_flags.h <==
#ifndef SAMMY_EXPERIMENT_FLAGS_H_INCLUDED_
#define SAMMY_EXPERIMENT_FLAGS_H_INCLUDED_

#include "basic_cuda_iteration.h"
#include "lower_bound_mip_settings.h"
#include "output.h"
#include "run_initial.h"
#include <boost/program_options.hpp>
#include <cstdlib>

namespace sammy {

template <typename ValueType> inline auto value_with_default(ValueType& vref) {
    ValueType default_value{vref};
    return boost::program_options::value(&vref)->default_value(default_value);
}

template <typename ValueType> inline auto required_value(ValueType& vref) {
    return boost::program_options::value(&vref)->required();
}

inline auto bool_switch(bool& bref) {
    return boost::program_options::bool_switch(&bref)->default_value(false);
}

struct ExperimentFlagsConfig {
    RunInitialConfig initial_heuristic_config;
    std::string output_file, input_file;
    std::string cuda_mode{"auto"};
    bool dont_simplify;
    bool initial_progress;
    bool print_events;
    LowerBoundMIPConfig lb_mip_config;
    SATColoringConfig sat_coloring_config;

    std::optional<std::string> process_parsed_args() {
        initial_heuristic_config.simplify = !dont_simplify;
        initial_heuristic_config.quiet = !initial_progress;
        if (!output_file.empty() && !recognized_output_extension(output_file)) {
            return "Error: Unrecognized output file extension (use .json or "
                   ".jsonl)!";
        }
        if (!std::filesystem::is_regular_file(input_file)) {
            return "Error: The given input file '" + input_file +
                   "' is not an existing regular file!";
        }
        if (cuda_mode == "off") {
            set_cuda_mode(CUDA_USAGE_DISABLED);
        } else if (cuda_mode == "forced") {
            set_cuda_mode(CUDA_USAGE_FORCED);
        } else {
            set_cuda_mode(CUDA_USAGE_IF_AVAILABLE);
        }
        return std::nullopt;
    }
};

inline void add_gurobi_experiment_flags(
    boost::program_options::options_description& description,
    ExperimentFlagsConfig& config) {
    auto& lcfg = config.lb_mip_config;
    description.add_options()(
        "lb-exactly-solve-full-mip", bool_switch(lcfg.solve_exactly_full_mip),
        "Try to solve the clique problem on the universe conflict graph "
        "exactly, "
        "using a MIP with additional cutting planes applied at the root node.")(
        "lb-exactly-branch-and-price",
        bool_switch(lcfg.solve_exactly_cut_and_price),
        "Try to solve the clique problem on the universe conflict graph "
        "exactly, "
        "using not a MIP but cutting planes and pricing for upper bounds and "
        "'accidental' integrality, rounding and other heuristics for lower "
        "bounds.")("lb-mip-stop-cuts-iteration-window",
                   value_with_default(lcfg.stop_cuts_iteration_window),
                   "Specify the length (in iterations) of the window to "
                   "consider when deciding "
                   "whether to stop adding cutting planes.")(
        "lb-mip-stop-cuts-below-relative-improvement",
        value_with_default(lcfg.stop_cuts_below_relative_improvement),
        "Specify the required relative improvement of the upper bound "
        "needed accross an iteration window below which generating cutting "
        "planes "
        "shall be stopped.")(
        "separation-rounds-before-pricing",
        value_with_default(lcfg.separation_rounds_before_pricing),
        "Specify the number of separation rounds to do between pricing rounds "
        "when faced with a non-optimal-for-subgraph LP relaxation.")(
        "lb-mip-timeout", value_with_default(lcfg.total_mip_timeout),
        "Specify the maximum amount of time to spend on lower bound MIP "
        "approach.")(
        "lb-mip-gurobi-quiet", bool_switch(lcfg.quiet_gurobi),
        "Tell Gurobi to be quiet while solving the lower bound MIPs/LPs.");
}

inline void
add_sat_coloring_flags(boost::program_options::options_description& description,
                       ExperimentFlagsConfig& config) {
    auto& scfg = config.sat_coloring_config;
    description.add_options()(
        "initial-sat-coloring", bool_switch(scfg.initial_sat_coloring),
        "Use SAT coloring for some time to initialize upper and lower bounds.")(
        "initial-sat-coloring-time-limit",
        value_with_default(scfg.initial_sat_coloring_timeout),
        "Time limit for initial SAT coloring.");
}

inline void
add_experiment_flags(boost::program_options::options_description& description,
                     ExperimentFlagsConfig& config,
                     bool enable_gurobi_clique_flags = false,
                     bool enable_sat_coloring_flags = false) {
    auto& icfg = config.initial_heuristic_config;
    description.add_options()("help,h", "produce help output")(
        "dont-simplify", bool_switch(config.dont_simplify),
        "do not apply simplification")("initial-time-limit",
                                       value_with_default(icfg.max_time),
                                       "Time limit for the initial heuristic.")(
        "cuda-mode", value_with_default(config.cuda_mode),
        "select CUDA mode: either 'off', 'force', or 'auto'. "
        "If 'off', no CUDA will be used. If 'force', CUDA will be used; "
        "an error is raised if this cannot be done. "
        "If 'auto' (or any other value), CUDA will be used if possible, "
        "falling back to CPU otherwise.")(
        "initial-iteration-goal", value_with_default(icfg.goal_iterations),
        "Iterations of the initial heuristic to run at least (unless the time "
        "limit is exceeded).")(
        "initial-goal-time", value_with_default(icfg.min_time),
        "Time (in seconds) to spend on the initial heuristic at least "
        "(unless the time limit is exceeded or optimality is proved).")(
        "initial-random-clique-restarts",
        value_with_default(icfg.random_clique_restarts_per_iteration),
        "Number of repeats (per thread) to run the randomized greedy clique "
        "algorithm for per iteration of the initial heuristic.")(
        "initial-heuristic-progress", bool_switch(config.initial_progress),
        "Print progress information during the initial heuristic")(
        "print-events", bool_switch(config.print_events),
        "Print time-stamped event information")(
        "output-file,o", value_with_default(config.output_file));
    if (enable_gurobi_clique_flags) {
        add_gurobi_experiment_flags(description, config);
    }
    if (enable_sat_coloring_flags) {
        add_sat_coloring_flags(description, config);
    }
}

inline void add_input_file_flag(
    boost::program_options::options_description& description,
    boost::program_options::positional_options_description& positional,
    ExperimentFlagsConfig& config) {
    description.add_options()(
        "input-file",
        boost::program_options::value<std::string>(&config.input_file)
            ->required(),
        "The input file to read.");
    positional.add("input-file", 1);
}

inline bool parse_cmdline(
    const boost::program_options::options_description& all_options,
    const boost::program_options::positional_options_description& positional,
    boost::program_options::variables_map& vmap, int argc, char** argv) {
    auto cmdline_parser = boost::program_options::command_line_parser(
        argc, const_cast<const char* const*>(argv));
    boost::program_options::store(
        cmdline_parser.options(all_options).positional(positional).run(), vmap);
    vmap.notify();
    return !vmap.count("help");
}

inline void get_config_or_exit(
    int argc, char** argv, ExperimentFlagsConfig& defaults_and_out,
    const boost::program_options::options_description& extra_options = {},
    bool enable_gurobi_clique_flags = false,
    bool enable_sat_coloring_flags = false) {
    boost::program_options::options_description visible("Allowed options");
    boost::program_options::options_description hidden("Hidden options");
    boost::program_options::positional_options_description positional;
    visible.add(extra_options);
    add_experiment_flags(visible, defaults_and_out, enable_gurobi_clique_flags,
                         enable_sat_coloring_flags);
    add_input_file_flag(hidden, positional, defaults_and_out);
    boost::program_options::options_description all;
    all.add(visible).add(hidden);
    boost::program_options::variables_map vmap;
    try {
        if (!parse_cmdline(all, positional, vmap, argc, argv)) {
            std::cout << visible << std::endl;
            std::exit(0);
        }
    } catch (const boost::program_options::error& err) {
        int ecode = 1;
        if (!vmap.count("help")) {
            std::cerr << "Invalid command line options!" << std::endl;
            std::cerr << err.what() << std::endl;
            ecode = 0;
        }
        std::cerr << visible << std::endl;
        std::exit(ecode);
    }
    defaults_and_out.process_parsed_args();
}

inline ExperimentFlagsConfig get_config_or_exit(
    int argc, char** argv,
    const boost::program_options::options_description& extra_options = {},
    bool enable_gurobi_clique_flags = false,
    bool enable_sat_coloring_flags = false) {
    ExperimentFlagsConfig config;
    get_config_or_exit(argc, argv, config, extra_options,
                       enable_gurobi_clique_flags, enable_sat_coloring_flags);
    return config;
}

} // namespace sammy

#endif
==> ./barrage.h <==
#ifndef SAMMY_BARRAGE_H_INCLUDED_
#define SAMMY_BARRAGE_H_INCLUDED_

#include "barrage_lns_subproblem.h"
#include "clique_storage.h"
#include "experiment_flags.h"
#include "implied_vertices.h"
#include "initial_coloring_heuristic.h"
#include "initial_phase_result.h"
#include "output.h"
#include "pair_infeasibility_map.h"
#include "simplification.h"
#include "thread_clauses.h"
#include <any>
#include <filesystem>
#include <optional>
#include <queue>
#include <utility>

namespace sammy {

using EventMask = std::uint32_t;

enum class PortfolioEvent : EventMask {
    TIMEOUT = 1,
    BETTER_LOWER_BOUND = 2,
    BETTER_UPPER_BOUND = 4,
    BETTER_MES = 8,
    ALARM = 16,
    OPTIMALITY = 32
};

struct InterruptionCheckInfo {
    std::size_t best_mes;
    std::size_t best_lower_bound;
    std::size_t best_upper_bound;
};

class PortfolioSolver;
struct LNSTimeAndSuccessInfo {
    struct RemovedClassesInfo {
        std::size_t complete_try_successes = 0;
        std::size_t complete_tries_since_last_success = 0;
        std::size_t complete_tries_total = 0;
        double complete_tries_total_time = 0.0;
    };

    inline LNSTimeAndSuccessInfo(PortfolioSolver* solver) noexcept;

    std::map<std::size_t, RemovedClassesInfo> removed_classes_info;
    double current_goal_time = 30.0;
    double num_failure_threshold = 20.0;
    std::size_t universe_size;
    std::size_t min_num_tries_for_time = 5;

    void report_success(std::size_t removed_classes, double time) {
        auto& info = removed_classes_info[removed_classes];
        info.complete_tries_since_last_success = 0;
        info.complete_tries_total++;
        info.complete_try_successes++;
        info.complete_tries_total_time += time;
    }

    void report_failure(std::size_t removed_classes, double time) {
        auto& info = removed_classes_info[removed_classes];
        info.complete_tries_since_last_success++;
        info.complete_tries_total++;
        info.complete_tries_total_time += time;
    }

    std::size_t select_goal_num_removed() {
        std::size_t r = p_select_goal_num_removed();
        return r;
    }

  private:
    std::size_t p_select_goal_num_removed() {
        if (removed_classes_info.empty()) {
            return p_empty_select_num_removed();
        }
        auto last_usable = removed_classes_info.end();
        for (auto it = removed_classes_info.begin(),
                  e = removed_classes_info.end();
             it != e; ++it)
        {
            const auto& info = it->second;
            if (info.complete_tries_since_last_success == 0 &&
                info.complete_tries_total > 0)
            {
                last_usable = it;
                break;
            }
            if (info.complete_tries_since_last_success >= num_failure_threshold)
            {
                continue;
            }
            if (info.complete_tries_total < min_num_tries_for_time) {
                last_usable = it;
                continue;
            }
            double avg_time_per_try =
                info.complete_tries_total_time / info.complete_tries_total;
            if (avg_time_per_try >= current_goal_time) {
                last_usable = it;
            }
        }
        if (last_usable == removed_classes_info.end()) {
            return p_enlarge_num_removed();
        }
        return last_usable->first;
    }

    std::size_t p_enlarge_num_removed() {
        return removed_classes_info.rbegin()->first + 1;
    }

    std::size_t p_empty_select_num_removed() {
        if (universe_size < 100'000)
            return 9;
        if (universe_size < 500'000)
            return 7;
        if (universe_size < 1'000'000)
            return 5;
        if (universe_size < 5'000'000)
            return 4;
        return 3;
    }
};

class PortfolioElement {
  public:
    explicit PortfolioElement(PortfolioSolver* solver) : solver(solver) {}

    /**
     * Destroy the portfolio element.
     * Joins if the worker has not yet been joined on.
     */
    virtual ~PortfolioElement() { join(); }

    /**
     * Wait for the worker thread to end.
     * This does not notify the worker at all,
     * so the worker must be notified beforehand.
     */
    void join() {
        if (m_worker.joinable()) {
            m_worker.join();
        }
    }

    /**
     * Get a thread-local copy of the clauses.
     */
    inline ClauseDB& get_clauses() const;

    /**
     * Get the infeasibility map.
     */
    inline PairInfeasibilityMap& get_infeasibility_map() noexcept;

    /*
     * Routine that is called by the event
     * delivery mechanism (coordinator thread)
     * to check if this worker should be interrupted,
     * and interrupt it if necessary.
     * Called with the mutex held.
     */
    virtual void interrupt_if_necessary(const InterruptionCheckInfo& /*info*/) {
    }

    /**
     * Get an event recorder (if this element has any).
     */
    virtual const EventRecorder* get_recorder() const { return nullptr; }

    /**
     * If this element has any, synchronize the event recorder.
     */
    virtual void synchronize_recorder(const EventRecorder& /*other*/) {}

    /**
     * If this element has any, set the quiet flag of the event recorder.
     */
    virtual void set_recorder_quiet(bool quiet) {}

    /**
     * Get a description of this element (it it has any).
     */
    virtual std::string get_description() const { return ""; }

    /**
     * Called by the coordinator to
     * deliver new events to this worker.
     */
    inline void deliver_events(EventMask new_events,
                               const InterruptionCheckInfo& info);

    /**
     * Called by the coordinator to deliver alarm events.
     */
    inline void deliver_alarm(std::size_t alarm_id,
                              const InterruptionCheckInfo& info);

    /**
     * Set an alarm (this causes the coordinator to notify us,
     * generating an ALARM event at the requested time).
     * Automatically discards the currently set alarm (if there is any).
     */
    inline void set_alarm(double in_seconds);

    /**
     * Discard the currently set alarm.
     */
    inline void discard_alarm();

    /**
     * Can be called to obtain the raised events
     * and reset them to their unraised state.
     */
    EventMask consume_events() {
        std::unique_lock l{mutex};
        EventMask result = events;
        events = 0;
        return result;
    }

    /**
     * Called to start the worker thread.
     */
    void start() {
        if (m_worker.joinable())
            throw std::logic_error("Trying to start running worker!");
        m_worker = std::thread{[&]() { this->main(); }};
    }

  protected:
    /**
     * Method implementing the actual work.
     */
    virtual void main() = 0;

    /**
     * Mutex that locks the 'events' field,
     * and any data that can be accessed by
     * multiple threads, e.g., data accessed
     * in the 'interrupt_if_necessary' method.
     */
    std::mutex mutex;

    /**
     * Condition variable that can be
     * waited on for events.
     */
    std::condition_variable condition;

    /**
     * The coordinating solver.
     */
    PortfolioSolver* solver;

    /**
     * The events currently raised/pending.
     */
    EventMask events = 0;

    /**
     * An atomic flag that we can check to see if
     * we should terminate ASAP. This flag is
     * already set on delivery (not consumption)
     * of the TIMEOUT event.
     */
    std::atomic<bool> should_terminate{false};

  private:
    std::thread m_worker; //< the worker thread
    std::optional<std::size_t> m_alarm{std::nullopt};
    std::size_t m_alarm_id_counter{0};
};

class PortfolioSolver {
  public:
    PortfolioSolver(ClausesTicket clauses, EventRecorder* global_recorder,
                    InitialPhaseResult&& initial_phase)
        : m_inf_map(std::move(initial_phase.inf_map)),
          m_implied_cache(&m_inf_map, initial_phase.universe_size),
          m_best_spawners(std::move(initial_phase.best_spawners)),
          m_all_spawners(std::move(initial_phase.all_spawners)),
          m_coloring_order(std::move(initial_phase.coloring_order)),
          m_universe_size(initial_phase.universe_size), m_ticket(clauses),
          m_global_recorder(global_recorder),
          m_best_mes(std::move(initial_phase.best_mutually_exclusive)),
          m_best_lb_vertex_set(m_best_mes), m_lower_bound(m_best_mes.size()),
          m_best_solution(local_clauses(clauses).num_vars(), &m_inf_map,
                          initial_phase.best_solution.begin(),
                          initial_phase.best_solution.end()),
          m_lns_info(this) {}

    void limited_reduce_universe(double time_limit) {
        if (!m_implied_cache.have_reduced_universe()) {
            m_global_recorder->store_event(
                "BEGIN_LIMITED_IMPLIED_VERTEX_ELIMINATION");
            m_implied_cache.limited_reduce_universe(local_clauses(m_ticket),
                                                    time_limit);
            p_post_reduce_universe();
        }
    }

    void reduce_universe() {
        if (!m_implied_cache.have_reduced_universe()) {
            m_global_recorder->store_event("BEGIN_IMPLIED_VERTEX_ELIMINATION");
            m_implied_cache.reduce_universe(local_clauses(m_ticket));
            p_post_reduce_universe();
        }
    }

    const ImpliedVertexCache& implied_cache() const noexcept {
        return m_implied_cache;
    }

    PartialSolution get_best_solution() const {
        std::unique_lock l{m_mutex};
        PartialSolution partial = m_best_solution;
        return partial;
    }

    const std::vector<Vertex>& get_best_spawners() const noexcept {
        return m_best_spawners;
    }

    const std::vector<Vertex>& get_all_spawners() const noexcept {
        return m_all_spawners;
    }

    const std::vector<Vertex>& get_coloring_order() const noexcept {
        return m_coloring_order;
    }

    std::size_t get_universe_size() const noexcept { return m_universe_size; }

    void set_alarm(PortfolioElement* element, double in_seconds,
                   std::size_t alarm_id) {
        std::unique_lock l{m_mutex};
        AlarmTime time =
            Clock::now() + std::chrono::duration<double>(in_seconds);
        m_alarm_requests.emplace_back(AlarmRequest{time, alarm_id, element});
        if (m_alarm_requests.size() == 1) {
            m_condition.notify_one();
        }
    }

    ClauseDB& get_clauses() const { return local_clauses(m_ticket); }

    PairInfeasibilityMap& get_infeasibility_map() noexcept { return m_inf_map; }

    void run(double time_limit = 1000.0 * 365.0 * 24.0 * 60.0 * 60.0) {
        set_alarm(nullptr, time_limit, 0);
        p_main_loop();
    }

    void merge_events() {
        std::size_t worker_id = 0;
        for (const auto& e : m_elements) {
            const auto* events = e->get_recorder();
            if (events) {
                std::string description = e->get_description();
                if (description.empty())
                    description = std::to_string(worker_id);
                m_global_recorder->merge_events(*events, description);
            }
            ++worker_id;
        }
    }

    void add_element(std::unique_ptr<PortfolioElement> element) {
        std::unique_lock l{m_mutex};
        m_waiting_elements.emplace_back(std::move(element));
        if (m_waiting_elements.size() == 1) {
            m_condition.notify_one();
        }
    }

    template <typename ConcreteType, typename... Args>
    ConcreteType& emplace_element(Args&&... args) {
        auto element =
            std::make_unique<ConcreteType>(std::forward<Args>(args)...);
        ConcreteType& ref = *element;
        add_element(std::move(element));
        return ref;
    }

    /**
     * @brief Record a new event in the global event recorder.
     */
    template <typename Tag, typename... EventArgs>
    void report_global(Tag&& tag, OutputObject o, EventArgs&&... p) {
        std::unique_lock l{m_mutex};
        m_global_recorder->store_event(std::forward<Tag>(tag), std::move(o),
                                       std::forward<EventArgs>(p)...);
    }

    /**
     * @brief Obtain a clique cache view.
     */
    CliqueStorage::StorageView clique_cache_view() const noexcept {
        return m_clique_cache.obtain_view();
    }

    /**
     * @brief Report a new mutually exclusive set.
     */
    void report_mes(const std::vector<Vertex>& vertices, const char* source) {
        std::unique_lock l{m_mutex};
        if (m_best_mes.size() >= vertices.size())
            return;
        m_best_mes = vertices;
        if (m_implied_cache.have_reduced_universe()) {
            m_implied_cache.replace_implied(m_best_mes);
            std::size_t old_size = m_best_mes.size();
            std::sort(m_best_mes.begin(), m_best_mes.end());
            m_best_mes.erase(std::unique(m_best_mes.begin(), m_best_mes.end()),
                             m_best_mes.end());
            if (m_best_mes.size() != old_size) {
                throw std::logic_error(
                    "MES size changed during implied vertex elimination!");
            }
        }
        EventMask event = static_cast<EventMask>(PortfolioEvent::BETTER_MES);
        m_global_recorder->store_event("IMPROVED_MES",
                                       {{"size", m_best_mes.size()},
                                        {"vertices", m_best_mes},
                                        {"source", source}},
                                       "size", "source");
        if (m_best_mes.size() > m_lower_bound) {
            m_lower_bound = m_best_mes.size();
            m_best_lb_vertex_set = m_best_mes;
            m_global_recorder->store_event("IMPROVED_LB",
                                           {{"lb", m_best_mes.size()},
                                            {"vertices", m_best_mes},
                                            {"source", source}},
                                           "lb", "source");
            event |= static_cast<EventMask>(PortfolioEvent::BETTER_LOWER_BOUND);
            if (m_lower_bound >= m_best_solution.size()) {
                m_global_recorder->store_event("OPTIMALITY_REACHED");
                event |= static_cast<EventMask>(PortfolioEvent::OPTIMALITY);
            }
        }
        p_raise_events(event);
    }

    /**
     * @brief Called by workers to report a new solution.
     * Returns true if the reported solution is better than the
     * current best and was stored in the portfolio solver.
     */
    bool report_solution(const PartialSolution& solution, const char* source) {
        std::unique_lock l{m_mutex};
        if (m_best_solution.size() <= solution.size())
            return false;
        m_best_solution = solution;
        m_global_recorder->store_event("IMPROVED_SOLUTION",
                                       {{"size", m_best_solution.size()},
                                        {"source", source},
                                        {"lb", m_lower_bound}},
                                       "size", "lb", "source");
        EventMask events =
            static_cast<EventMask>(PortfolioEvent::BETTER_UPPER_BOUND);
        if (m_best_solution.size() <= m_lower_bound) {
            m_global_recorder->store_event("OPTIMALITY_REACHED");
            events |= static_cast<EventMask>(PortfolioEvent::OPTIMALITY);
        }
        p_raise_events(events);
        return true;
    }

    /**
     * @brief Called by workers to report new lower bounds.
     */
    void report_lower_bound(std::size_t lower_bound,
                            const std::vector<Vertex>& subgraph,
                            const char* source) {
        std::unique_lock l{m_mutex};
        if (lower_bound <= m_lower_bound)
            return;
        m_best_lb_vertex_set = subgraph;
        m_lower_bound = lower_bound;
        m_global_recorder->store_event("IMPROVED_LB",
                                       {{"lb", m_lower_bound},
                                        {"source", source},
                                        {"vertices", m_best_lb_vertex_set},
                                        {"ub", m_best_solution.size()}},
                                       "lb", "source", "ub");
        EventMask events =
            static_cast<EventMask>(PortfolioEvent::BETTER_LOWER_BOUND);
        if (m_best_solution.size() <= m_lower_bound) {
            m_global_recorder->store_event("OPTIMALITY_REACHED");
            events |= static_cast<EventMask>(PortfolioEvent::OPTIMALITY);
        }
        p_raise_events(events);
    }

    /**
     * Add a clique to the cache.
     */
    void add_clique(const std::vector<Vertex>& vertices) {
        m_clique_cache.push_clique(vertices.begin(), vertices.end());
    }

    /**
     * Add a clique to the cache.
     */
    template <typename Iterator> void add_clique(Iterator begin, Iterator end) {
        m_clique_cache.push_clique(begin, end);
    }

    /**
     * Mark a clique as used.
     */
    void clique_was_used(CliqueStorage::StorageView& view,
                         CliqueStorage::StorageView::Iterator iter) {
        m_clique_cache.used_clique(view, iter);
    }

    /**
     * Mark a clique as used.
     */
    void clique_was_used(CliqueStorage::StorageView& view, Index index) {
        m_clique_cache.used_clique(view, index);
    }

    /**
     * `Get the best mutually exclusive set currently available.
     */
    std::vector<Vertex> get_best_mes() const {
        std::unique_lock l{m_mutex};
        std::vector<Vertex> result = m_best_mes;
        return result;
    }

    /**
     * Get the best mutually exclusive set size currently available.
     */
    std::size_t get_best_mes_size() const noexcept {
        std::unique_lock l{m_mutex};
        return m_best_mes.size();
    }

    /**
     * Get the best lower bound currently available.
     */
    std::size_t get_best_lower_bound() const noexcept {
        std::unique_lock l{m_mutex};
        return m_lower_bound;
    }

    std::vector<Vertex> get_best_lower_bound_vertex_set() const {
        std::unique_lock l{m_mutex};
        return m_best_lb_vertex_set;
    }

    /**
     * Get the best upper bound currently available.
     */
    std::size_t get_best_solution_size() const noexcept {
        std::unique_lock l{m_mutex};
        return m_best_solution.size();
    }

    std::size_t lns_select_removed_class_count() {
        std::unique_lock l{m_mutex};
        return m_lns_info.select_goal_num_removed();
    }

    void lns_report_failure(std::size_t num_removed, double time) {
        std::unique_lock l{m_mutex};
        m_lns_info.report_failure(num_removed, time);
    }

    void lns_report_success(std::size_t num_removed, double time) {
        std::unique_lock l{m_mutex};
        m_lns_info.report_success(num_removed, time);
    }

    void enable_subproblem_reporting(const std::filesystem::path& out_dir) {
        if (m_subproblem_dir) {
            throw std::logic_error("Subproblem reporting already enabled!");
        }
        if (!exists(out_dir)) {
            create_directory(out_dir);
        }
        m_subproblem_dir = out_dir;
        auto universe_and_clauses_file = out_dir / "universe_and_clauses.json";
        nlohmann::json universe_and_clauses;
        auto universe = lit::externalize(m_implied_cache.get_universe());
        const auto& clauses = local_clauses(m_ticket);
        universe_and_clauses["infeasibility_map"] = m_inf_map.export_bits();
        universe_and_clauses["best_solution"] =
            m_best_solution.assignments_as<std::vector<bool>>();
        universe_and_clauses["best_spawners"] =
            lit::externalize(m_best_spawners);
        universe_and_clauses["best_mutually_exclusive"] =
            lit::externalize(m_best_mes);
        universe_and_clauses["all_spawners"] = lit::externalize(m_all_spawners);
        universe_and_clauses["coloring_order"] =
            lit::externalize(m_coloring_order);
        universe_and_clauses["universe_size"] = m_universe_size;
        universe_and_clauses["universe"] = std::move(universe);
        universe_and_clauses["clauses"] = clauses.export_all_clauses();
        universe_and_clauses["num_variables"] = clauses.num_vars();
        universe_and_clauses["num_concrete"] = m_inf_map.get_n_concrete();
        std::ofstream outfile;
        outfile.exceptions(std::ios::badbit | std::ios::failbit);
        outfile.open(universe_and_clauses_file,
                     std::ios::out | std::ios::trunc);
        outfile << universe_and_clauses;
    }

    bool subproblem_reporting_enabled() const noexcept {
        return m_subproblem_dir.has_value();
    }

    std::uint64_t new_subproblem_id() noexcept {
        return m_subproblem_counter++;
    }

    void report_subproblem(const LNSSubproblem& subproblem,
                           const PartialSolution& remaining_configs,
                           std::size_t global_best_mes_size,
                           std::size_t global_best_lower_bound,
                           const char* subproblem_type_) {
        if (!subproblem_reporting_enabled())
            return;
        auto id = new_subproblem_id();
        std::ostringstream name_fmt;
        std::string subproblem_type(subproblem_type_);
        auto name_filter = [](const std::string& str) {
            std::string result;
            for (char c : str) {
                if (std::isalnum(c)) {
                    result.push_back(c);
                } else {
                    result.push_back('_');
                }
            }
            return result;
        };
        name_fmt << "subproblem-" << name_filter(subproblem_type) << "-"
                 << std::setw(5) << std::setfill('0') << id << ".json";
        auto subproblem_file = *m_subproblem_dir / name_fmt.str();
        nlohmann::json output;
        output["id"] = id;
        output["subproblem_type"] = subproblem_type;
        output["uncovered"] = lit::externalize(subproblem.uncovered_universe);
        output["initial_uncovered_mes"] =
            lit::externalize(subproblem.mutually_exclusive_set);
        output["num_nonremoved_configs"] = remaining_configs.size();
        output["num_removed_configs"] =
            subproblem.removed_configurations.size();
        output["lns_info"] = p_lns_info_to_json();
        output["removed_configs"] =
            bitsets_to_json(subproblem.removed_configurations);
        output["global_best_mes_size"] = global_best_mes_size;
        output["global_best_lower_bound"] = global_best_lower_bound;
        if (!remaining_configs.empty()) {
            output["remaining_config"] =
                std::vector<bool>(remaining_configs.get_assignment(0));
        } else {
            output["remaining_config"] = nlohmann::json{};
        }
        std::ofstream outfile;
        outfile.exceptions(std::ios::badbit | std::ios::failbit);
        outfile.open(subproblem_file, std::ios::out | std::ios::trunc);
        outfile << output;
    }

    void report_subproblem(const std::vector<Vertex>& uncovered,
                           const std::vector<Vertex>& initial_uncovered_mes,
                           const PartialSolution& remaining_configs,
                           const std::vector<DynamicBitset>& removed_configs,
                           std::size_t global_best_mes_size,
                           std::size_t global_best_lower_bound,
                           const char* subproblem_type) {
        if (!subproblem_reporting_enabled())
            return;
        auto id = new_subproblem_id();
        std::ostringstream name_fmt;
        name_fmt << "subproblem-" << subproblem_type << "-" << std::setw(5)
                 << std::setfill('0') << id << ".json";
        auto subproblem_file = *m_subproblem_dir / name_fmt.str();
        nlohmann::json subproblem;
        subproblem["id"] = id;
        subproblem["subproblem_type"] = subproblem_type;
        subproblem["uncovered"] = lit::externalize(uncovered);
        subproblem["initial_uncovered_mes"] =
            lit::externalize(initial_uncovered_mes);
        subproblem["num_nonremoved_configs"] = remaining_configs.size();
        subproblem["num_removed_configs"] = removed_configs.size();
        subproblem["lns_info"] = p_lns_info_to_json();
        subproblem["removed_configs"] = bitsets_to_json(removed_configs);
        subproblem["global_best_mes_size"] = global_best_mes_size;
        subproblem["global_best_lower_bound"] = global_best_lower_bound;
        if (!remaining_configs.empty()) {
            subproblem["remaining_config"] =
                std::vector<bool>(remaining_configs.get_assignment(0));
        } else {
            subproblem["remaining_config"] = nlohmann::json{};
        }
        std::ofstream outfile;
        outfile.exceptions(std::ios::badbit | std::ios::failbit);
        outfile.open(subproblem_file, std::ios::out | std::ios::trunc);
        outfile << subproblem;
    }

  private:
    void p_raise_events(EventMask events) {
        EventMask prev = m_raised;
        if ((prev | events) == prev)
            return;
        m_raised |= events;
        m_condition.notify_one();
    }

    using AlarmTime =
        decltype(Clock::now() + std::chrono::duration<double>(1.0));
    struct AlarmRequest {
        AlarmTime time;
        std::size_t id;
        PortfolioElement* element;

        bool operator<(const AlarmRequest& other) const noexcept {
            return other.time < time;
        }
    };

    void p_post_reduce_universe() {
        m_global_recorder->store_event(
            "DONE_IMPLIED_VERTEX_ELIMINATION",
            {{"original_size", m_implied_cache.original_universe_size()},
             {"reduced_size", m_implied_cache.reduced_universe_size()}},
            "original_size", "reduced_size");
        m_implied_cache.replace_implied(m_best_spawners);
        sort_unique(m_best_spawners);
        m_implied_cache.replace_implied(m_all_spawners);
        sort_unique(m_all_spawners);
        m_implied_cache.replace_implied(m_coloring_order);
        nosort_unique<PairHashSet<Vertex>>(m_coloring_order);
        m_implied_cache.replace_implied(m_best_mes);
        std::size_t old_mes_size = m_best_mes.size();
        sort_unique(m_best_mes);
        if (old_mes_size != m_best_mes.size()) {
            throw std::logic_error(
                "MES size changed during implied vertex elimination!");
        }
    }

    void p_update_alarms() {
        for (const auto& r : m_alarm_requests) {
            m_alarm_queue.push(r);
        }
        m_alarm_requests.clear();
    }

    bool p_should_wake() const {
        return !m_alarm_requests.empty() || !m_waiting_elements.empty() ||
               m_raised != 0;
    }

    void p_extract_due_alarms() {
        auto now = Clock::now();
        while (!m_alarm_queue.empty() && m_alarm_queue.top().time <= now) {
            m_due_alarms.push_back(m_alarm_queue.top());
            m_alarm_queue.pop();
        }
    }

    bool p_handle_due_alarms(const InterruptionCheckInfo& info) {
        for (const AlarmRequest& a : m_due_alarms) {
            if (!a.element) {
                p_raise_timeout();
                return true;
            } else {
                a.element->deliver_alarm(a.id, info);
            }
        }
        m_due_alarms.clear();
        return false;
    }

    void p_raise_timeout() {
        InterruptionCheckInfo info;
        {
            std::unique_lock l{m_mutex};
            m_shutting_down = true;
            info = p_get_info();
            m_global_recorder->store_event("PORTFOLIO_TIMEOUT");
        }
        p_deliver_events(static_cast<EventMask>(PortfolioEvent::TIMEOUT), info);
    }

    void p_deliver_events(EventMask events, const InterruptionCheckInfo& info) {
        for (auto& e : m_elements) {
            e->deliver_events(events, info);
        }
    }

    void p_update_elements() {
        if (!m_shutting_down) {
            for (auto& e : m_waiting_elements) {
                e->start();
                m_elements.emplace_back(std::move(e));
            }
            m_waiting_elements.clear();
        }
    }

    InterruptionCheckInfo p_main_locked(EventMask& new_events) {
        std::unique_lock l{m_mutex};
        p_update_elements();
        p_update_alarms();
        if (m_alarm_queue.empty()) {
            m_condition.wait(l, [&]() { return p_should_wake(); });
        } else {
            m_condition.wait_until(l, m_alarm_queue.top().time,
                                   [&]() { return p_should_wake(); });
        }
        new_events = m_raised;
        m_raised = 0;
        p_update_elements();
        p_update_alarms();
        return p_get_info();
    }

    OutputObject p_lns_info_to_json() const {
        std::unique_lock l{m_mutex};
        OutputObject result{
            {"current_goal_time", m_lns_info.current_goal_time},
            {"num_failure_threshold", m_lns_info.num_failure_threshold},
            {"min_num_tries_for_time", m_lns_info.min_num_tries_for_time},
            {"removed_classes_info", nlohmann::json{}}};
        auto& removed_classes_info = result["removed_classes_info"];
        for (const auto& [num_removed, info] : m_lns_info.removed_classes_info)
        {
            nlohmann::json info_json{
                {"num_removed", num_removed},
                {"complete_try_successes", info.complete_try_successes},
                {"complete_tries_since_last_success",
                 info.complete_tries_since_last_success},
                {"complete_tries_total", info.complete_tries_total},
                {"complete_tries_total_time", info.complete_tries_total_time}};
            removed_classes_info.push_back(std::move(info_json));
        }
        return result;
    }

    bool p_main_unlocked(const InterruptionCheckInfo& info,
                         EventMask new_events) {
        p_extract_due_alarms();
        if (p_handle_due_alarms(info)) {
            return true;
        }
        if (new_events) {
            p_deliver_events(new_events, info);
            if (new_events & static_cast<EventMask>(PortfolioEvent::OPTIMALITY))
            {
                std::unique_lock l{m_mutex};
                m_shutting_down = true;
                return true;
            }
        }
        return false;
    }

    void p_await_completion() {
        for (auto& e : m_elements) {
            e->join();
        }
    }

    void p_main_loop() {
        for (;;) {
            EventMask new_events = 0;
            auto info = p_main_locked(new_events);
            if (p_main_unlocked(info, new_events)) {
                break;
            }
        }
        p_await_completion();
    }

    InterruptionCheckInfo p_get_info() const {
        return InterruptionCheckInfo{m_best_mes.size(), m_lower_bound,
                                     m_best_solution.size()};
    }

    // --- read-only elements ---
    /**
     * Pre-filled pair infeasibility map.
     */
    PairInfeasibilityMap m_inf_map;

    /**
     * ImpliedVertexCache for reducing the universe size.
     */
    ImpliedVertexCache m_implied_cache;

    /**
     * Class spawners for best sample from initial phase.
     */
    std::vector<Vertex> m_best_spawners;

    /**
     * Class spawners for all samples from the initial phase.
     */
    std::vector<Vertex> m_all_spawners;

    /**
     * Coloring order for the last heuristic run of the initial phase.
     */
    std::vector<Vertex> m_coloring_order;

    /**
     * The number of valid interactions.
     */
    std::size_t m_universe_size;

    /**
     * Ticket to access thread-local
     * mutable versions of the clause set.
     */
    ClausesTicket m_ticket;

    /**
     * Directory to which the subproblems are to be reported.
     */
    std::optional<std::filesystem::path> m_subproblem_dir;

    /**
     * Atomic counter for subproblem ids.
     */
    std::atomic<std::uint64_t> m_subproblem_counter{0};

    // --- events and coordination ---
    /**
     * Mutex for all global data and events.
     */
    mutable std::mutex m_mutex;

    /**
     * A flag that tracks if we are currently shutting down.
     */
    bool m_shutting_down{false};

    /**
     * Event recorder for output/logging purposes.
     */
    EventRecorder* m_global_recorder;

    /**
     * Condition variable to wait for events, alarm requests, ...
     */
    mutable std::condition_variable m_condition;

    /**
     * Pending alarm requests (new alarm requests go here).
     */
    std::vector<AlarmRequest> m_alarm_requests;

    /**
     * Unlocked container for alarms that have expired,
     * but have not yet been dispatched.
     */
    std::vector<AlarmRequest> m_due_alarms;

    /**
     * Priority queue of alarms. Contains a 'pseudo-alarm'
     * elements that signifies timeout expiry.
     */
    std::priority_queue<AlarmRequest> m_alarm_queue;

    /**
     * The set of raised events that have to be dispatched.
     */
    EventMask m_raised = 0;

    /**
     * The owning container of PortfolioElements that
     * are currently active.
     * It is not allowed to delete from this set while
     * the algorithm is running.
     */
    std::vector<std::unique_ptr<PortfolioElement>> m_elements;

    /**
     * List of elements that have been created and are waiting
     * to activate (and be started).
     */
    std::vector<std::unique_ptr<PortfolioElement>> m_waiting_elements;

    // --- best solution, LB vertex set, MES ---
    /**
     * Best mutually exclusive set.
     */
    std::vector<Vertex> m_best_mes;

    /**
     * Vertex set of best LB (these vertices require at least m_lower_bound
     * configurations). Need not be equal to m_best_mes.
     */
    std::vector<Vertex> m_best_lb_vertex_set;

    /**
     * Numerical LB value.
     */
    std::size_t m_lower_bound;

    /**
     * Best full solution.
     */
    PartialSolution m_best_solution;

    /**
     * Storage cache for cliques.
     * They are used for the destroy heuristic.
     */
    CliqueStorage m_clique_cache;

    /**
     * Information for deciding how many configurations
     * to remove in the destroy operation.
     */
    LNSTimeAndSuccessInfo m_lns_info;
};

void PortfolioElement::set_alarm(double in_seconds) {
    std::size_t alarm_id;
    {
        std::unique_lock l{mutex};
        if (m_alarm) {
            m_alarm.reset();
        }
        events &= ~static_cast<EventMask>(PortfolioEvent::ALARM);
        alarm_id = m_alarm_id_counter++;
        m_alarm = alarm_id;
    }
    solver->set_alarm(this, in_seconds, alarm_id);
}

void PortfolioElement::discard_alarm() {
    std::unique_lock l{mutex};
    if (m_alarm) {
        m_alarm.reset();
    }
    events &= ~static_cast<EventMask>(PortfolioEvent::ALARM);
}

void PortfolioElement::deliver_alarm(std::size_t alarm_id,
                                     const InterruptionCheckInfo& info) {
    std::unique_lock l{mutex};
    if (!m_alarm || alarm_id != *m_alarm) {
        return;
    }
    bool was_empty = (events == 0);
    events |= static_cast<EventMask>(PortfolioEvent::ALARM);
    interrupt_if_necessary(info);
    if (was_empty) {
        condition.notify_one();
    }
}

void PortfolioElement::deliver_events(EventMask new_events,
                                      const InterruptionCheckInfo& info) {
    std::unique_lock l{mutex};
    bool was_empty = (events == 0);
    events |= new_events;
    if (new_events & (static_cast<EventMask>(PortfolioEvent::TIMEOUT) |
                      static_cast<EventMask>(PortfolioEvent::OPTIMALITY)))
    {
        should_terminate.store(true);
    }
    interrupt_if_necessary(info);
    if (was_empty) {
        condition.notify_one();
    }
}

PairInfeasibilityMap& PortfolioElement::get_infeasibility_map() noexcept {
    return solver->get_infeasibility_map();
}

ClauseDB& PortfolioElement::get_clauses() const {
    return solver->get_clauses();
}

LNSTimeAndSuccessInfo::LNSTimeAndSuccessInfo(PortfolioSolver* solver) noexcept
    : universe_size(solver->get_universe_size()) {}

} // namespace sammy

#endif
==> ./universe_subgraph.h <==
#ifndef SAMMY_UNIVERSE_SUBGRAPH_H_INCLUDED_
#define SAMMY_UNIVERSE_SUBGRAPH_H_INCLUDED_

#include "clique_or_indset_builder.h"
#include "dynamic_bitset.h"
#include "pair_infeasibility_map.h"
#include "parallel_bit_filter.h"
#include "shared_db_propagator.h"
#include "thread_group.h"
#include "thread_interrupt.h"
#include "vertex_operations.h"

namespace sammy {

/**
 * UniverseSubgraph represents a growable subgraph of
 * the universe of interactions with an edge between
 * two interactions if they are mutually exclusive.
 */
class UniverseSubgraph {
  public:
    UniverseSubgraph(const UniverseSubgraph& o)
        : propagator(o.propagator), infeasibility_map(o.infeasibility_map),
          vertices(o.vertices), matrix(o.matrix),
          vertices_with_literal(o.vertices_with_literal),
          vertex_index_map(o.vertex_index_map), degree(o.degree),
          parallel_bits(&o.parallel_bits.thread_group()) {}

    UniverseSubgraph& operator=(const UniverseSubgraph& o) {
        propagator = o.propagator;
        infeasibility_map = o.infeasibility_map;
        vertices = o.vertices;
        matrix = o.matrix;
        vertices_with_literal = o.vertices_with_literal;
        vertex_index_map = o.vertex_index_map;
        degree = o.degree;
        return *this;
    }

    UniverseSubgraph(UniverseSubgraph&& o) noexcept
        : propagator(std::move(o.propagator)),
          infeasibility_map(o.infeasibility_map),
          vertices(std::move(o.vertices)), matrix(std::move(o.matrix)),
          vertices_with_literal(std::move(o.vertices_with_literal)),
          vertex_index_map(std::move(o.vertex_index_map)),
          degree(std::move(o.degree)),
          parallel_bits(&o.parallel_bits.thread_group()) {}

    UniverseSubgraph& operator=(UniverseSubgraph&& o) noexcept {
        std::swap(propagator, o.propagator);
        std::swap(infeasibility_map, o.infeasibility_map);
        std::swap(vertices, o.vertices);
        std::swap(matrix, o.matrix);
        std::swap(vertices_with_literal, o.vertices_with_literal);
        std::swap(vertex_index_map, o.vertex_index_map);
        std::swap(degree, o.degree);
        return *this;
    }

    bool operator==(const UniverseSubgraph& o) const noexcept {
        return vertices == o.vertices && matrix == o.matrix;
    }

    bool operator!=(const UniverseSubgraph& o) const noexcept {
        return !(*this == o);
    }

    UniverseSubgraph(ClauseDB* all_clauses, ThreadGroup<void>* thread_pool,
                     const PairInfeasibilityMap* infeasibility_map,
                     std::vector<Vertex> vertices)
        : propagator(all_clauses), infeasibility_map(infeasibility_map),
          vertices(std::move(vertices)),
          matrix(this->vertices.size(),
                 DynamicBitset(this->vertices.size(), false)),
          vertices_with_literal(2 * all_clauses->num_vars(),
                                DynamicBitset(this->vertices.size(), false)),
          parallel_bits(thread_pool) {
        p_build_vertex_index_map();
        p_build_vertices_with_literal();
        p_build_matrix();
    }

    SharedDBPropagator& get_propagator() noexcept { return propagator; }

    const SharedDBPropagator& get_propagator() const noexcept {
        return propagator;
    }

    DynamicBitset& matrix_row(std::size_t index) noexcept {
        return matrix[index];
    }

    const DynamicBitset& matrix_row(std::size_t index) const noexcept {
        return matrix[index];
    }

    std::size_t n() const noexcept { return vertices.size(); }

    const std::vector<std::size_t>& get_degrees() const noexcept {
        return degree;
    }

    const std::vector<Vertex>& vertex_set() const noexcept { return vertices; }

    Vertex vertex(std::size_t index) const noexcept { return vertices[index]; }

    bool has_vertex(Vertex v) const noexcept {
        return vertex_index_map.count(v);
    }

    std::size_t vertex_index(Vertex v) const { return vertex_index_map.at(v); }

    using CliqueBuilder = sammy::CliqueBuilder<UniverseSubgraph>;
    using IndependentSetBuilder =
        sammy::IndependentSetBuilder<UniverseSubgraph>;

    CliqueBuilder clique_builder() {
        return CliqueBuilder{this, &parallel_bits};
    }

    IndependentSetBuilder independent_set_builder() {
        return IndependentSetBuilder{this, &parallel_bits};
    }

    bool is_edge(Vertex v, Vertex w) const noexcept {
        return matrix[vertex_index(v)][vertex_index(w)];
    }

    bool is_edge(std::size_t i1, std::size_t i2) const noexcept {
        return matrix[i1][i2];
    }

    UniverseSubgraph restrict_to_degree_at_least(std::size_t min_degree) {
        std::vector<std::size_t> new_degrees = degree;
        std::vector<std::size_t> removal_queue;
        for (std::size_t i = 0, n = this->n(); i < n; ++i) {
            if (new_degrees[i] < min_degree) {
                removal_queue.push_back(i);
            }
        }
        std::size_t queue_pos = 0;
        while (queue_pos < removal_queue.size()) {
            std::size_t r = removal_queue[queue_pos++];
            for (std::size_t neighbor : matrix[r].ones()) {
                if (new_degrees[neighbor]-- == min_degree) {
                    removal_queue.push_back(neighbor);
                }
            }
        }
        std::sort(removal_queue.begin(), removal_queue.end());
        std::vector<Vertex> new_vertices;
        std::vector<DynamicBitset> new_matrix;
        std::vector<DynamicBitset> new_vertices_with_literal;
        removal_queue.push_back(this->n());
        p_walk_removals(removal_queue, new_vertices, new_matrix,
                        new_vertices_with_literal);
        return UniverseSubgraph(propagator, infeasibility_map,
                                std::move(new_vertices), std::move(new_matrix),
                                std::move(new_vertices_with_literal),
                                parallel_bits.thread_group());
    }

    /**
     * @brief Add a bulk of new vertices.
     *
     * @param vertices
     */
    void add_vertices(const std::vector<Vertex>& new_vertices) {
        std::size_t old_n = vertices.size();
        std::size_t new_n = old_n + new_vertices.size();
        std::size_t res_size = 1;
        while (res_size <= new_n) {
            res_size <<= 1;
        }
        for (auto& row : matrix) {
            row.reserve(res_size);
            row.resize(new_n, false);
        }
        matrix.reserve(res_size);
        matrix.resize(new_n, DynamicBitset(new_n, false));
        for (auto& vwl : vertices_with_literal) {
            vwl.reserve(res_size);
            vwl.resize(new_n, false);
        }
        vertices.insert(vertices.end(), new_vertices.begin(),
                        new_vertices.end());
        for (std::size_t i = old_n; i < new_n; ++i) {
            vertex_index_map.try_emplace(vertices[i], i);
        }
        degree.resize(vertices.size(), 0);
        p_extend_vertices_with_literal(old_n, new_n);
        p_extend_matrix(old_n, new_n);
    }

    Var get_n_concrete() const noexcept {
        return infeasibility_map->get_n_concrete();
    }

    const PairInfeasibilityMap* get_infeasibility_map() const noexcept {
        return infeasibility_map;
    }

    void nonedge_to_edge(std::size_t index1, std::size_t index2) {
        matrix[index1][index2] = true;
        matrix[index2][index1] = true;
        ++degree[index1];
        ++degree[index2];
    }

    template <typename InputIterator>
    void add_new_vertices(InputIterator begin, InputIterator end) {
        std::vector<Vertex> vnew;
        std::copy_if(begin, end, std::back_inserter(vnew),
                     [&](Vertex v) { return !has_vertex(v); });
        add_vertices(vnew);
    }

    template <typename VertexContainer>
    std::vector<std::size_t>
    to_indices(const VertexContainer& vertices) const noexcept {
        std::vector<std::size_t> result;
        result.reserve(vertices.size());
        std::transform(vertices.begin(), vertices.end(),
                       std::back_inserter(result),
                       [&](Vertex v) { return vertex_index(v); });
        return result;
    }

    /**
     * @brief Beside the edges implied by conflicts found in any case,
     *        we may also try to extend the information by
     *        checking, for each pair that is not already excluded,
     *        whether both vertices can be pushed simultaneously.
     *        This is MUCH more expensive than building the basic matrix.
     * @return The number of additional edges found.
     */
    std::size_t extend_matrix_by_propagation(bool interruptible = false) {
        if (is_extended) {
            return 0;
        }
        std::size_t count = 0, check_count = 0;
        for (std::size_t i = 1, nv = n(); i < nv; ++i) {
            Vertex v = vertices[i];
            reset_and_push_noresolve(propagator, v);
            DynamicBitset& row_i = matrix[i];
            for (std::size_t j = 0; j < i; ++j) {
                if (!row_i[j]) {
                    Vertex w = vertices[j];
                    if (!can_push(propagator, w)) {
                        row_i[j] = true;
                        matrix[j][i] = true;
                        ++degree[i];
                        ++degree[j];
                        ++count;
                    }
                }
            }
            if (interruptible && ++check_count % 1024 == 0) {
                throw_if_interrupted();
            }
        }
        is_extended = true;
        return count;
    }

    bool is_extended_by_propagation() const noexcept { return is_extended; }

    ThreadGroup<void>& thread_group() noexcept {
        return parallel_bits.thread_group();
    }

    std::size_t get_degree(std::size_t vertex_index) const noexcept {
        return degree[vertex_index];
    }

    std::size_t get_degree(Vertex v) const noexcept {
        return degree[vertex_index(v)];
    }

    BitsetOperationsBuffer* get_parallel_bits() noexcept {
        return &parallel_bits;
    }

  private:
    UniverseSubgraph(const SharedDBPropagator& propagator,
                     const PairInfeasibilityMap* infeasibility_map,
                     std::vector<Vertex> vertices,
                     std::vector<DynamicBitset> matrix,
                     std::vector<DynamicBitset> vertices_with_literal,
                     ThreadGroup<void>& thread_pool)
        : propagator(propagator), infeasibility_map(infeasibility_map),
          vertices(std::move(vertices)), matrix(std::move(matrix)),
          vertices_with_literal(std::move(vertices_with_literal)),
          degree(this->vertices.size(), 0), parallel_bits(&thread_pool) {
        this->propagator.reset_or_throw();
        p_build_vertex_index_map();
        for (std::size_t i = 0, n = this->n(); i < n; ++i) {
            degree[i] = this->matrix[i].count();
        }
    }

    void p_walk_removals_vertices(const std::vector<std::size_t>& removed,
                                  std::vector<Vertex>& nvertices) {
        std::size_t new_n = vertices.size() - removed.size();
        nvertices.reserve(new_n);
        auto rem_iter = removed.begin();
        for (std::size_t i = 0, nold = n(); i < nold; ++i) {
            if (i == *rem_iter) {
                ++rem_iter;
                continue;
            }
            nvertices.push_back(vertices[i]);
        }
    }

    void p_walk_removals_matrix(const std::vector<std::size_t>& removed,
                                std::vector<DynamicBitset>& nmatrix) {
        auto out_rem_iter = removed.begin();
        const std::size_t nold = n();
        std::size_t out_new_i = 0;
        for (std::size_t out_i = 0; out_i < nold; ++out_i) {
            if (out_i == *out_rem_iter) {
                ++out_rem_iter;
                continue;
            }
            auto& row = nmatrix[out_new_i];
            const auto& old_row = matrix[out_i];
            auto in_rem_iter = removed.begin();
            std::size_t in_new_i = 0;
            for (std::size_t in_i = 0; in_i < nold; ++in_i) {
                if (in_i == *in_rem_iter) {
                    ++in_rem_iter;
                    continue;
                }
                if (old_row[in_i]) {
                    row[in_new_i].set();
                }
                ++in_new_i;
            }
            ++out_new_i;
        }
    }

    void p_walk_removals_vwl(const std::vector<std::size_t>& removed,
                             std::vector<DynamicBitset>& nv_with_lit) {
        Lit nalit = 2 * propagator.db().num_vars();
        for (Lit l = 0; l < nalit; ++l) {
            auto rem_iter = removed.begin();
            std::size_t new_i = 0;
            const auto& old_row = vertices_with_literal[l];
            auto& new_row = nv_with_lit[l];
            for (std::size_t i = 0, n = this->n(); i < n; ++i) {
                if (i == *rem_iter) {
                    ++rem_iter;
                    continue;
                }
                if (old_row[i]) {
                    new_row[new_i].set();
                }
                ++new_i;
            }
        }
    }

    void p_walk_removals(const std::vector<std::size_t>& removed,
                         std::vector<Vertex>& nvertices,
                         std::vector<DynamicBitset>& nmatrix,
                         std::vector<DynamicBitset>& nv_with_lit) {
        p_walk_removals_vertices(removed, nvertices);
        auto nalit = 2 * propagator.db().num_vars();
        std::size_t new_n = nvertices.size();
        nmatrix.assign(new_n, DynamicBitset(new_n, false));
        nv_with_lit.assign(nalit, DynamicBitset(new_n, false));
        p_walk_removals_matrix(removed, nmatrix);
        p_walk_removals_vwl(removed, nv_with_lit);
    }

    void p_build_vertex_index_map() {
        vertex_index_map.reserve(vertices.size());
        for (std::size_t vi = 0, n = vertices.size(); vi < n; ++vi) {
            if (!vertex_index_map.try_emplace(vertices[vi], vi).second) {
                throw std::logic_error(
                    "Duplicate vertex given to UniverseSubgraph!");
            }
        }
    }

    void p_assert_symmetry() {
        for (std::size_t i = 1, n = this->n(); i < n; ++i) {
            const auto& row = matrix[i];
            for (std::size_t j = 0; j < i; ++j) {
                if (bool(row[j]) != bool(matrix[j][i])) {
                    std::cerr << "Asymmetry: matrix[" << i << "," << j
                              << "] != matrix[" << j << "," << i << "]!\n"
                              << std::flush;
                    std::abort();
                }
            }
        }
    }

    void p_build_matrix_row_range(std::size_t begin, std::size_t end,
                                  std::size_t old_n) {
        SharedDBPropagator local_prop = propagator;
        const std::size_t n = vertices.size();
        for (std::size_t vi = begin; vi != end; ++vi) {
            Vertex v = vertices[vi];
            auto& row = matrix[vi];
            reset_and_push_noresolve(local_prop, v);
            for (Lit lpos : local_prop.get_trail()) {
                Lit lneg = lit::negate(lpos);
                row |= vertices_with_literal[lneg];
            }
            if (is_extended) {
                // already checked this for old rows
                // during extension of existing rows
                for (std::size_t vj = 0; vj != old_n; ++vj) {
                    if (!row[vj]) {
                        row[vj] = matrix[vj][vi];
                    }
                }
                for (std::size_t vj = old_n; vj != vi; ++vj) {
                    if (row[vj])
                        continue;
                    Vertex w = vertices[vj];
                    if (!can_push(local_prop, w))
                        row[vj] = true;
                }
                for (std::size_t vj = vi + 1; vj < n; ++vj) {
                    if (row[vj])
                        continue;
                    Vertex w = vertices[vj];
                    if (!can_push(local_prop, w))
                        row[vj] = true;
                }
            }
            degree[vi] = row.count();
        }
    }

    void p_extend_matrix_row_range(std::size_t begin, std::size_t end,
                                   std::size_t start_from) {
        SharedDBPropagator local_prop = propagator;
        for (std::size_t vi = begin; vi != end; ++vi) {
            Vertex v = vertices[vi];
            reset_and_push_noresolve(local_prop, v);
            for (Lit lpos : local_prop.get_trail()) {
                Lit lneg = lit::negate(lpos);
                matrix[vi].binary_or(vertices_with_literal[lneg], start_from);
            }
            if (is_extended) {
                auto& row = matrix[vi];
                for (std::size_t vj = start_from, n = vertices.size(); vj != n;
                     ++vj)
                {
                    if (row[vj])
                        continue;
                    Vertex w = vertices[vj];
                    if (!can_push(local_prop, w)) {
                        row[vj] = true;
                    }
                }
            }
            degree[vi] = matrix[vi].count();
        }
    }

    void p_build_matrix(std::size_t row_begin, std::size_t row_end) {
        propagator.reset_or_throw();
        std::size_t num_threads = std::thread::hardware_concurrency();
        std::size_t tp_threads = parallel_bits.thread_group().num_threads() + 1;
        std::size_t num_vertices = row_end - row_begin;
        num_threads = (std::min)(num_threads, tp_threads);
        num_threads = (std::min)(num_threads, num_vertices);
        if (num_threads == 1) {
            p_build_matrix_row_range(row_begin, row_end, row_begin);
            return;
        }
        std::unique_ptr<std::thread[]> builders =
            std::make_unique<std::thread[]>(num_threads);
        std::size_t vertices_per_thread = num_vertices / num_threads;
        for (std::size_t i = 0, cur = row_begin; i < num_threads;
             ++i, cur += vertices_per_thread)
        {
            std::size_t cend = cur + vertices_per_thread;
            if (i == num_threads - 1)
                cend = row_end;
            builders[i] = std::thread(
                [this, &row_begin](std::size_t b, std::size_t e) {
                    p_build_matrix_row_range(b, e, row_begin);
                },
                cur, cend);
        }
        for (std::size_t i = 0; i < num_threads; ++i) {
            builders[i].join();
        }
        static_cast<void>(&UniverseSubgraph::p_assert_symmetry);
        assert((p_assert_symmetry(),
                "Graph/adjacency matrix should be symmetric!"));
    }

    void p_build_matrix() {
        degree.resize(vertices.size(), 0);
        p_build_matrix(0, vertices.size());
    }

    void p_extend_existing_rows(std::size_t existing_rows,
                                std::size_t begin_new) {
        std::size_t num_threads = std::thread::hardware_concurrency();
        std::size_t tp_threads = parallel_bits.thread_group().num_threads() + 1;
        num_threads = (std::min)(num_threads, tp_threads);
        num_threads = (std::min)(num_threads, existing_rows);
        if (num_threads == 1) {
            p_extend_matrix_row_range(0, existing_rows, begin_new);
            return;
        }
        std::unique_ptr<std::thread[]> builders =
            std::make_unique<std::thread[]>(num_threads);
        std::size_t vertices_per_thread = existing_rows / num_threads;
        for (std::size_t i = 0, cur = 0; i < num_threads;
             ++i, cur += vertices_per_thread)
        {
            std::size_t cend = cur + vertices_per_thread;
            if (i == num_threads - 1)
                cend = existing_rows;
            builders[i] = std::thread(
                [&](std::size_t b, std::size_t e) {
                    p_extend_matrix_row_range(b, e, begin_new);
                },
                cur, cend);
        }
        for (std::size_t i = 0; i < num_threads; ++i) {
            builders[i].join();
        }
    }

    void p_extend_matrix(std::size_t old_n, std::size_t new_n) {
        propagator.reset_or_throw();
        p_extend_existing_rows(old_n, old_n);
        p_build_matrix(old_n, new_n);
    }

    void p_extend_vertices_with_literal(std::size_t vindex_begin,
                                        std::size_t vindex_end) {
        for (std::size_t vi = vindex_begin; vi != vindex_end; ++vi) {
            Vertex v = vertices[vi];
            reset_and_push_noresolve(propagator, v);
            for (Lit l : propagator.get_trail()) {
                vertices_with_literal[l][vi].set();
            }
        }
    }

    void p_build_vertices_with_literal() {
        p_extend_vertices_with_literal(0, vertices.size());
    }

    SharedDBPropagator propagator;
    const PairInfeasibilityMap* infeasibility_map;
    std::vector<Vertex> vertices;
    std::vector<DynamicBitset> matrix;
    std::vector<DynamicBitset> vertices_with_literal;
    VertexMapTo<std::size_t> vertex_index_map;
    std::vector<std::size_t> degree;
    mutable BitsetOperationsBuffer parallel_bits;
    bool is_extended = false;
};

} // namespace sammy

#endif
==> ./compress_binaries.h <==
#ifndef SAMMY_COMPRESS_BINARIES_H_INCLUDED_
#define SAMMY_COMPRESS_BINARIES_H_INCLUDED_

#include "detect_equalities.h"
#include "dynamic_bitset.h"
#include "literals.h"
#include "shared_db_propagator.h"
#include "simplify_datastructure.h"
#include "stamp_set.h"
#include "thread_group.h"

namespace sammy {

class BinaryClauseCompressor {
  public:
    explicit BinaryClauseCompressor(SimplifyDatastructure* simplify_ds)
        : simplify_ds(simplify_ds), m_clauses(simplify_ds->extract_clause_db()),
          m_out_edges(2 * simplify_ds->original_num_vars()),
          m_visited(2 * simplify_ds->original_num_vars()),
          m_children(2 * simplify_ds->original_num_vars()),
          m_remove(2 * simplify_ds->original_num_vars()) {}

    void compute_binary_clause_graph() {
        for (const auto& cl : simplify_ds->clauses()) {
            if (cl.size() != 2)
                continue;
            Lit l1 = cl[0], l2 = cl[1];
            m_out_edges[lit::negate(l1)].push_back(l2);
            m_out_edges[lit::negate(l2)].push_back(l1);
        }
    }

    void transitively_reduce_dag() {
        std::vector<std::pair<Lit, Lit>> binaries;
        const Lit nl = 2 * simplify_ds->original_num_vars();
        for (Lit l = 0; l < nl; ++l) {
            if (simplify_ds->is_eliminated(lit::var(l))) {
                ++l;
                continue;
            }
            m_visited.clear();
            m_children.clear();
            m_remove.clear();
            const auto& succs = m_out_edges[l];
            std::for_each(succs.begin(), succs.end(),
                          [&](Lit succ) { m_children.insert(succ); });
            for (Lit succ : succs) {
                if (!m_visited.check_insert(succ)) {
                    m_remove.insert(succ);
                    continue;
                }
                m_dfs_stack.push_back(succ);
                while (!m_dfs_stack.empty()) {
                    Lit curr = m_dfs_stack.back();
                    m_dfs_stack.pop_back();
                    for (Lit next : m_out_edges[curr]) {
                        if (!m_visited.check_insert(next)) {
                            if (m_children.count(next)) {
                                m_remove.insert(next);
                            }
                        } else {
                            m_dfs_stack.push_back(next);
                        }
                    }
                }
            }
            for (Lit succ : succs) {
                if (!m_remove.count(succ)) {
                    Lit l1 = lit::negate(l), l2 = succ;
                    binaries.emplace_back((std::min)(l1, l2),
                                          (std::max)(l1, l2));
                }
            }
        }
        std::sort(binaries.begin(), binaries.end());
        binaries.erase(std::unique(binaries.begin(), binaries.end()),
                       binaries.end());
        simplify_ds->replace_all_binaries(binaries);
    }

    void compute_transitive_implication_graph(ThreadGroup<void>* tgroup) {
        SharedDBPropagator propagator_ctx{&m_clauses};
        tgroup->parallel_foreach_iterator(
            Lit(0), Lit(2 * simplify_ds->original_num_vars()), propagator_ctx,
            [&](SharedDBPropagator& propagator, Lit l) {
                if (simplify_ds->is_eliminated(lit::var(l)) ||
                    !propagator.is_open(l))
                {
                    return;
                }
                if (!propagator.push_level(l)) {
                    throw std::logic_error("BinaryClauseCompressor should only "
                                           "be used after simplification!");
                }
                auto& out_edges = m_out_edges[l];
                auto implied = detail::implied_literals(propagator);
                std::copy(implied.begin(), implied.end(),
                          std::back_inserter(out_edges));
                propagator.pop_level();
            });
    }

    void transitive_reduction() {
        std::vector<std::pair<Lit, Lit>> binaries;
        for (Lit l = 0, nl = 2 * simplify_ds->original_num_vars(); l != nl; ++l)
        {
            if (simplify_ds->is_eliminated(lit::var(l)))
                continue;
            m_visited.clear();
            for (Lit succ : m_out_edges[l]) {
                for (Lit succ2 : m_out_edges[succ]) {
                    m_visited.insert(succ2);
                }
            }
            for (Lit succ : m_out_edges[l]) {
                if (!m_visited.count(succ)) {
                    Lit l1 = lit::negate(l);
                    Lit l2 = succ;
                    binaries.emplace_back(std::min(l1, l2), std::max(l1, l2));
                }
            }
        }
        std::sort(binaries.begin(), binaries.end());
        binaries.erase(std::unique(binaries.begin(), binaries.end()),
                       binaries.end());
        simplify_ds->replace_all_binaries(binaries);
    }

  private:
    using Row = DynamicBitset;
    using Matrix = std::vector<Row>;
    SimplifyDatastructure* simplify_ds;
    ClauseDB m_clauses;
    std::vector<std::vector<Lit>> m_out_edges;
    std::vector<Lit> m_dfs_stack;
    StampSet<Lit, std::uint16_t> m_visited;
    StampSet<Lit, std::uint16_t> m_children;
    StampSet<Lit, std::uint16_t> m_remove;
};

} // namespace sammy

#endif
==> ./barrage_worker_with_core.h <==
#ifndef SAMMY_BARRAGE_WORKER_WITH_CORE_H_INCLUDED_
#define SAMMY_BARRAGE_WORKER_WITH_CORE_H_INCLUDED_

#include "barrage.h"
#include "thread_interrupt.h"

namespace sammy {

template <typename CoreType>
class PortfolioElementWithCore : public PortfolioElement {
  public:
    using CoreFactoryType = std::function<std::unique_ptr<CoreType>(
        PortfolioSolver*, PortfolioElementWithCore*)>;

    PortfolioElementWithCore(PortfolioSolver* solver,
                             CoreFactoryType core_factory,
                             const std::string& description)
        : PortfolioElement(solver), m_core_factory(std::move(core_factory)),
          m_description(description) {}

    EventRecorder* get_mutable_recorder() { return &m_recorder; }

    std::mutex& get_mutex() noexcept { return mutex; }

  protected:
    void main() override {
        if (should_terminate.load()) {
            return;
        }
        set_interrupt_flag_ptr(&should_terminate);
        try {
            p_construct_core();
        } catch (const InterruptError&) {
            return;
        }
        if (should_terminate.load()) {
            return;
        }
        m_core->main();
    }

    void interrupt_if_necessary(const InterruptionCheckInfo& info) override {
        if (m_core) {
            if (should_terminate.load()) {
                m_core->termination_flag_set();
            } else {
                m_core->interrupt_if_necessary(info);
            }
        }
        // otherwise, we are constructing the core and,
        // if termination is wanted, should_terminate is automatically set
        // which causes the thread to interrupt
        events = 0;
    }

    const EventRecorder* get_recorder() const override { return &m_recorder; }

    /**
     * If this element has any, synchronize the event recorder.
     */
    virtual void synchronize_recorder(const EventRecorder& other) override {
        m_recorder.synchronize_with(other);
    }

    /**
     * If this element has any, set the quiet flag of the event recorder.
     */
    virtual void set_recorder_quiet(bool quiet) override {
        m_recorder.set_print_events(!quiet);
    }

    /**
     * Get a description of this element (it it has any).
     */
    virtual std::string get_description() const override {
        return m_description;
    }

    std::unique_ptr<CoreType> m_core;
    CoreFactoryType m_core_factory;
    EventRecorder m_recorder;
    std::string m_description;

  private:
    void p_construct_core() {
        std::unique_ptr<CoreType> core = m_core_factory(this->solver, this);
        {
            std::unique_lock l{mutex};
            m_core = std::move(core);
        }
    }
};

} // namespace sammy

#endif
==> ./all.md <==
==> ./sammy.h <==
#ifndef SAMMY_H_INCLUDED_
#define SAMMY_H_INCLUDED_

#endif
==> ./simplification.h <==
#ifndef SAMMY_SIMPLIFICATION_H_INCLUDED_
#define SAMMY_SIMPLIFICATION_H_INCLUDED_

#include "bounded_variable_elimination.h"
#include "compress_binaries.h"
#include "detect_equalities.h"
#include "eliminate_subsumed.h"
#include "simplify_datastructure.h"
#include "vivify.h"

namespace sammy {

inline SimplifiedInstance run_simplifier(SimplifyDatastructure& simplifier,
                                         SimplificationStats& stats) {
    bool changed;
    Var nv = simplifier.original_num_vars();
    auto before = Clock::now();
    eliminate_subsumed(simplifier.clauses(), nv, &stats);
    assert(!simplifier.has_duplicate_binary_clause());
    while (detect_failed_and_equal_literals(simplifier, &stats))
        ;
    eliminate_subsumed(simplifier.clauses(), nv, &stats);
    assert(!simplifier.has_duplicate_binary_clause());
    do {
        changed = bounded_variable_elimination(simplifier, 20, &stats);
        if (changed)
            eliminate_subsumed(simplifier.clauses(), nv, &stats);
        assert(!simplifier.has_duplicate_binary_clause());
        changed |= vivify(simplifier, &stats);
        changed |= detect_failed_and_equal_literals(simplifier, &stats);
        if (changed)
            eliminate_subsumed(simplifier.clauses(), nv, &stats);
        assert(!simplifier.has_duplicate_binary_clause());
        stats.simplification_rounds += 1;
    } while (changed);
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
inline SimplifiedInstance remove_subsumed(const SimplifiedInstance& instance) {
    const auto n_all = instance.formula.num_vars();
    auto clause_list = to_clause_list<CVec>(instance.formula);
    eliminate_subsumed(clause_list, n_all);
    return SimplifiedInstance{instance.new_to_old, ClauseDB(n_all, clause_list),
                              instance.num_concrete};
}

} // namespace sammy

#endif
==> ./thread_interrupt.h <==
#ifndef SAMMY_THREAD_INTERRUPT_H_INCLUDED_
#define SAMMY_THREAD_INTERRUPT_H_INCLUDED_

#include <atomic>
#include <exception>
#include <stdexcept>

namespace sammy {

class InterruptError : public std::exception {
  public:
    const char* what() const noexcept override { return "Interrupted"; }
};

inline std::atomic<bool>*& get_interrupt_flag_ptr() noexcept {
    thread_local std::atomic<bool>* flag = nullptr;
    return flag;
}

inline void set_interrupt_flag_ptr(std::atomic<bool>* flag) noexcept {
    get_interrupt_flag_ptr() = flag;
}

inline bool peek_interrupt_flag() noexcept {
    std::atomic<bool>* iflag = get_interrupt_flag_ptr();
    if (!iflag)
        return false;
    return iflag->load();
}

inline bool get_and_clear_interrupt_flag() noexcept {
    std::atomic<bool>* iflag = get_interrupt_flag_ptr();
    if (!iflag)
        return false;
    bool fv = iflag->load();
    if (!fv)
        return false;
    fv = iflag->exchange(false);
    return fv;
}

inline void trigger_interrupt() {
    std::atomic<bool>* iflag = get_interrupt_flag_ptr();
    if (!iflag)
        throw std::logic_error("Interrupt flag not set!");
    iflag->store(true);
}

inline void throw_if_interrupted() {
    if (get_and_clear_interrupt_flag()) {
        throw InterruptError{};
    }
}

} // namespace sammy

#endif
==> ./gurobi_clique_solver_g2.h <==
#ifndef SAMMY_GUROBI_CLIQUE_SOLVER_G2_H_INCLUDED_
#define SAMMY_GUROBI_CLIQUE_SOLVER_G2_H_INCLUDED_

#include "fast_clique.h"
#include "gurobi.h"
#include "initial_coloring_heuristic.h"
#include "literals.h"
#include "lower_bound_mip_settings.h"
#include "output.h"
#include "rng.h"
#include "universe_subgraph.h"
#include "vertex_operations.h"

namespace sammy {

/**
 * Gurobi cut-and-price clique solver working on an extensible
 * subgraph of the conflicts in the interaction universe.
 */
class GurobiCliqueSolverG2 {
  public:
    template <typename BitsetType>
    GurobiCliqueSolverG2(UniverseSubgraph* subgraph, Var n_all,
                         EventRecorder* recorder,
                         const std::vector<BitsetType>& initial_sample,
                         const std::vector<Vertex>& initial_mut_ex_set,
                         const LowerBoundMIPConfig& config)
        : m_subgraph(subgraph), m_recorder(recorder), m_config(config),
          m_model(sammy::gurobi_environment(config.quiet_gurobi)),
          m_nall(n_all), m_indset_builder(subgraph->independent_set_builder()),
          m_clique_builder(subgraph->clique_builder()),
          m_subgraph_upper_bound(initial_sample.size()),
          m_best_sample(p_to_sample(initial_sample)),
          m_best_mes(initial_mut_ex_set) {
        m_model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
        m_model.set(GRB_IntParam_Threads,
                    int(m_subgraph->thread_group().num_threads() + 1));
        m_model.set(GRB_IntParam_Method, GRB_METHOD_DUAL);
        p_create_initial_vars();
        add_configurations(initial_sample.begin(), initial_sample.end());
    }

    /**
     * Make this solver abortable.
     * This means that the 'abort' method can
     * be used to attempt to abort long-running methods
     * before they are completed.
     * Such methods will behave as if a timeout occurred.
     */
    void make_abortable() {
        class GurobiAbortCallback : public GRBCallback {
          public:
            GurobiAbortCallback(GurobiCliqueSolverG2* that) : that(that) {}

            void callback() override {
                if (that->m_aborted.load()) {
                    this->abort();
                }
            }

          private:
            GurobiCliqueSolverG2* that;
        };
        if (!m_abortable) {
            m_abortable = true;
            m_callback = std::make_unique<GurobiAbortCallback>(this);
            m_model.setCallback(m_callback.get());
        }
    }

    /**
     * Attempt to stop long-running methods at the next opportunity.
     * Such methods will return a timeout result.
     */
    void abort() {
        if (!m_abortable)
            throw std::logic_error("Trying to abort non-abortable solver!");
        m_aborted.store(true);
    }

    /**
     * Get the thread group for this solver.
     */
    ThreadGroup<void>& thread_group() const noexcept {
        return m_subgraph->thread_group();
    }

    std::size_t num_lp_variables() const {
        return std::size_t(m_model.get(GRB_IntAttr_NumVars));
    }

    std::size_t num_lp_constraints() const {
        return std::size_t(m_model.get(GRB_IntAttr_NumConstrs));
    }

    std::size_t num_lp_nonzeros() const {
        return std::size_t(m_model.get(GRB_IntAttr_NumNZs));
    }

    std::size_t num_lp_solve_calls() const { return m_num_solve_calls; }

    std::size_t num_graph_vertices() const noexcept { return m_subgraph->n(); }

    /**
     * Add configurations as constraints.
     */
    template <typename ConfigurationIterator>
    inline void add_configurations(ConfigurationIterator cbeg,
                                   ConfigurationIterator cend);

    /**
     * Add a configuration as constraint.
     */
    template <typename ConfigurationType>
    inline void add_configuration(const ConfigurationType& cfg);

    /**
     * Add a set of vertices to the model and the subgraph.
     * Does not add vertices to independent sets, but adds
     * vertex variables to all independent set constraints
     * for independent sets the vertices are already added to.
     */
    inline void add_vertices(const std::vector<Vertex>& vertices);

    /**
     * Add the vertices from vertices using add_vertices that are not already in
     * the graph.
     */
    inline void add_new_vertices(const std::vector<Vertex>& vertices) {
        std::vector<Vertex> actual_vertices;
        std::copy_if(vertices.begin(), vertices.end(),
                     std::back_inserter(actual_vertices),
                     [&](Vertex v) { return !m_subgraph->has_vertex(v); });
        if (!actual_vertices.empty())
            add_vertices(actual_vertices);
    }

    /**
     * Solve the linear relaxation of the current model
     * and try to extract (essentially by greedy rounding) a
     * new clique from the solution. The solution may
     * still have violated non-edges, but any extracted clique is valid.
     */
    inline SolverState solve_relaxation(
        double time_allowed = std::numeric_limits<double>::infinity());

    /**
     * Solve the full relaxation, i.e., the linear relaxation including
     * all non-edge constraints on the current subgraph.
     */
    inline SolverState solve_full_relaxation(
        double time_allowed = std::numeric_limits<double>::infinity());

    /**
     * Solve the problem as MIP.
     * Temporarily transforms the problem into a MIP.
     */
    inline SolverState
    solve_as_mip(double time_allowed = std::numeric_limits<double>::infinity());

    /**
     * Find and prohibit violated nonedge constraints.
     * Returns true if a violated nonedge was found.
     */
    inline bool prohibit_violated_nonedges();

    /**
     * Add an independent set encoded by a propagator.
     * If the propagator encodes a complete configuration,
     * this adds a complete independent set; otherwise, it
     * adds an incomplete set that may be extended later.
     */
    inline void add_potentially_incomplete_independent_set(
        const SharedDBPropagator& propagator);

    /**
     * Add an independent set given as list of vertex indices.
     */
    inline void
    add_direct_independent_set(const std::vector<std::size_t>& indices);

    /**
     * Price a given range of vertices, and add those that could potentially
     * improve the optimal solution to the current relaxation.
     * Return true if it found any potentially interesting vertices, and false
     * if the current relaxation is optimal even if all priced vertices were
     * to be included in the current subgraph.
     */
    template <typename VertexIterator>
    inline std::size_t price_vertices(VertexIterator vbeg, VertexIterator vend);

    /**
     * Get the best mutually exclusive set found so far.
     */
    const std::vector<Vertex>& get_best_mes() const noexcept {
        return m_best_mes;
    }

    /**
     * Get the best solution, i.e., the best sample.
     */
    const std::vector<std::vector<bool>>& get_best_solution() const noexcept {
        return m_best_sample;
    }

    /**
     * Update the best (primal) solution.
     */
    void
    update_best_solution(const std::vector<std::vector<bool>>& new_solution) {
        if (new_solution.size() >= m_best_sample.size())
            return;
        m_best_sample = new_solution;
        if (m_subgraph_upper_bound > m_best_sample.size()) {
            m_subgraph_upper_bound = m_best_sample.size();
        }
    }

    double get_last_value() const noexcept { return m_last_objective_value; }

    std::size_t get_subgraph_ub() const noexcept {
        return m_subgraph_upper_bound;
    }

    inline void count_compressible_constraints(
        std::size_t min_iterations_without_relevance) {
        const auto ncmp = m_complete_independent_sets.size();
        const auto ninc = m_inc_independent_sets.size();
        std::size_t ccmp = 0, cinc = 0;
        for (std::size_t i = 0; i < ncmp; ++i) {
            if (m_complete_independent_set_zero_since[i] >=
                min_iterations_without_relevance)
            {
                ++ccmp;
            }
        }
        for (std::size_t i = 0; i < ninc; ++i) {
            if (m_inc_independent_set_zero_since[i] >=
                min_iterations_without_relevance)
            {
                ++cinc;
            }
        }
        double percent_reduction = 100.0 * (ccmp + cinc) / double(ncmp + ninc);
        m_recorder->store_event("COUNTED_COMPRESSIBLE",
                                {{"incomplete_reducible", cinc},
                                 {"complete_reducible", ccmp},
                                 {"percent_reducible", percent_reduction}},
                                "incomplete_reducible", "complete_reducible",
                                "percent_reducible");
    }

    /**
     * Try to greedily add new vertices to existing independent sets,
     * strengthening the corresponding constraints. Usually should be
     * tried before generating new independent sets.
     * Returns true if one of the strengthened constraints is now violated.
     */
    inline bool greedy_add_to_cutting_planes();

    /**
     * Try to greedily generate new, violated independent set constraints.
     * Returns true if a violated new constraint could be found.
     */
    inline bool greedy_generate_cutting_planes();

    /**
     * Incorporate a mutually exclusive set found externally.
     */
    inline void external_improved_mes(std::vector<Vertex> external_mes);

    /**
     * Add the (up to) num_constraints incomplete or complete
     * independent sets with the highest dual weights to the given
     * primal solver as initial partial configurations.
     */
    inline void
    export_highest_weights_to_primal(ColoringHeuristicSolver& primal,
                                     std::size_t num_constraints);

    /**
     * Get the support of the current fractional solution,
     * i.e., the vertices with non-zero value.
     */
    inline std::vector<Vertex> get_fractional_mes_support() const;

    /**
     * Use a vector of complete configurations to try to find violated cutting
     * planes.
     */
    inline std::size_t
    cuts_from_primal(const std::vector<SharedDBPropagator>& primal_configs,
                     std::size_t max_num_cuts);

    /**
     * Use a sample of complete configurations to try to find violated cutting
     * planes.
     */
    inline std::size_t
    cuts_from_sample(const std::vector<std::vector<bool>>& sample);

  private:
    // subgraph we are working on
    UniverseSubgraph* m_subgraph;
    // event recorder
    EventRecorder* m_recorder;
    // configuration values
    LowerBoundMIPConfig m_config;

    // gurobi model
    GRBModel m_model;

    // collection of as-of-yet incomplete independent sets
    // with buffers for their dual values and constraint handles
    std::unique_ptr<double[]> m_inc_independent_dual{nullptr};
    std::vector<GRBConstr> m_inc_independent_set_constraints;
    std::vector<SharedDBPropagator> m_inc_independent_sets;
    std::vector<std::uint32_t> m_inc_independent_set_zero_since;

    // collection of complete independent sets (i.e., configurations)
    // with buffers for their dual values and constraint handles
    std::unique_ptr<double[]> m_complete_independent_dual{nullptr};
    std::vector<GRBConstr> m_complete_independent_set_constraints;
    std::vector<DynamicBitset> m_complete_independent_sets;
    std::vector<std::uint32_t> m_complete_independent_set_zero_since;

    /**
     * Encoding of a directly-represented independent set.
     */
    struct DirectIndependentSet {
        // vertices in the independent set
        HashSet<std::size_t> vertex_indices;

        // literals are incompatible if any vertex
        // implies their negation
        DynamicBitset compatible_literals;
    };

    // collection of directly-represented independent sets
    // with buffers for their dual values and constraint handles
    std::unique_ptr<double[]> m_direct_independent_dual{nullptr};
    std::vector<GRBConstr> m_direct_independent_set_constraints;
    std::vector<DirectIndependentSet> m_direct_independent_sets;

    // number of features
    Var m_nall;

    // vector of all variables and buffer of primal values in last solution
    std::vector<GRBVar> m_vars;
    std::unique_ptr<double[]> m_last_solution_values{nullptr};

    using IndependentSetBuilder = UniverseSubgraph::IndependentSetBuilder;
    using CliqueBuilder = UniverseSubgraph::CliqueBuilder;
    IndependentSetBuilder m_indset_builder;
    CliqueBuilder m_clique_builder;

    // the currently extracted solution (vertex_index, value pairs), sorted by
    // decreasing value
    using OrderedSolution = std::vector<std::pair<std::size_t, double>>;
    using OSolIter = OrderedSolution::const_iterator;
    OrderedSolution m_ordered_solution;

    // the last objective value
    double m_last_objective_value;
    // the best integral value we can hope for on the current subgraph
    std::size_t m_subgraph_upper_bound;
    // if we have added new constraints by cutting, we want to solve
    // the relaxation to have dual values
    bool m_have_new_constraints = false;

    // the best sample we currently know
    std::vector<std::vector<bool>> m_best_sample;
    // the best mutually exclusive set we currently know
    std::vector<Vertex> m_best_mes;

    // buffer for a gurobi linear expression
    GRBLinExpr m_expr_buffer;

    // buffer for violated nonedges
    std::vector<std::pair<OSolIter, OSolIter>> m_violated_nonedges;

    // possibly, a callback if we're doing a MIP solve
    // or an abortable solve
    std::unique_ptr<GRBCallback> m_callback;
    std::unique_ptr<GRBCallback> m_old_callback;

    // mutex for when we need it during parallel operations
    std::mutex ds_lock;

    // if this solver is made abortable,
    // m_abortable is true. if the search
    // shall be aborted, m_aborted is set to true
    bool m_abortable{false};
    std::atomic<bool> m_aborted{false};

    // number of calls to optimize()
    std::size_t m_num_solve_calls{0};

    /**
     * Initialize the vector of variables.
     */
    inline void p_create_initial_vars();

    /**
     * Execute callable once for each literal.
     */
    template <typename Callable>
    inline void p_foreach_literal(Callable&& callable);

    const std::vector<std::vector<bool>>&
    p_to_sample(const std::vector<std::vector<bool>>& s) {
        return s;
    }

    std::vector<std::vector<bool>>
    p_to_sample(const std::vector<DynamicBitset>& s) {
        std::vector<std::vector<bool>> result;
        std::transform(s.begin(), s.end(), std::back_inserter(result),
                       [](const DynamicBitset& b) {
                           return static_cast<std::vector<bool>>(b);
                       });
        return result;
    }

    /**
     * Extract an ordered solution from the model, i.e.,
     * refresh/replace m_ordered_solution.
     */
    inline void
    p_extract_ordered_solution(bool possible_suboptimal_mip = false);

    /**
     * Turn a fractional bound on the objective value into an integer one,
     * leaving a small gap for numerical issues.
     */
    inline std::size_t
    p_make_objective_value_integral(double fractional) const noexcept;

    /**
     * Try to use rounding to turn the current fractional solution into a
     * clique.
     */
    inline void p_ordered_solution_to_clique();

    inline void p_begin_mip();
    inline void p_end_mip();

    /**
     * Identify violated nonedges.
     */
    inline bool p_identify_violated_nonedges();

    /**
     * Attempt to prohibit a violated non-edge by extending an
     * existing constraint to cover that non-edge.
     */
    inline bool p_prohibit_extend_existing(std::size_t v_index,
                                           DynamicBitset& covered_vne);

    /**
     * Update the lists of dual values.
     */
    inline void p_update_dual_values();

    /**
     * Add the given vertices to independent sets.
     * If the score of vertices increases to a safe level,
     * may delete some vertices from the list.
     */
    inline void
    p_add_vertices_to_independent_sets(std::vector<Vertex>& vertices,
                                       VertexMapTo<double>& vertex_scores);

    /**
     * Compute the total dual value of all existing
     * independent sets already containing v.
     */
    inline double p_price_in_existing(Vertex v) const noexcept;

    /**
     * Refresh an incomplete independent set constraint,
     * i.e., regenerate it based on its independent set
     * (to handle the case of indirectly added vertices).
     * Called in greedy constraint strengthening
     * when we detect an independent set to be violated.
     */
    inline void p_refresh_inc_constraint(std::size_t inc_index);

    /**
     * Collect vertices (either 'good' vertices, which correspond
     * to a dual constraint violated by at least 0.01, or 'possible'
     * vertices which are not explicitly in the current subgraph
     * but are only within 0.01 of corresponding to a violated dual constraint).
     */
    template <typename VertexIterator>
    inline bool
    p_collect_good_or_possible(VertexIterator vbeg, VertexIterator vend,
                               std::vector<Vertex>& good_vertices,
                               std::vector<Vertex>& possible_vertices,
                               VertexMapTo<double>& interesting_scores);

    /**
     * Categorize the given good vertices depending on their reduced costs,
     * and use that information to reduce the number of good vertices.
     */
    inline void p_categorize_good_vertices(
        std::vector<Vertex>& good, VertexMapTo<double>& interesting_scores,
        std::size_t num_categories, std::size_t goal_vertices);

    /**
     * Extract a complete independent set (constraint) and add it as initially
     * present configuration to the given primal heuristic solver.
     */
    inline void p_complete_to_primal(ColoringHeuristicSolver& primal,
                                     const DynamicBitset& complete_indset);

    /**
     * Extract an incomplete independent set (constraint) and add it as
     * initially present configuration to the given primal heuristic solver.
     */
    inline void p_incomplete_to_primal(ColoringHeuristicSolver& primal,
                                       const SharedDBPropagator& propagator);

    /**
     * Generate directly-represented independent sets.
     * Used if cheap cut procedures fail.
     */
    inline bool p_greedy_generate_direct_independent_sets();

    /**
     * Check if a given index can be added to a direct independent set
     * that it is not already part of.
     */
    inline bool p_direct_is_compatible(const DirectIndependentSet& direct,
                                       std::size_t index);

    /**
     * Strengthen direct independent set constraints
     * by adding vertices to them, based on the current
     * relaxed value of the vertices.
     */
    inline std::size_t p_greedy_add_to_direct();
};

namespace detail {

inline bool is_in_complete(const DynamicBitset& c, Vertex v) {
    using sammy::lit::negative;
    using sammy::lit::var;
    return c[var(v.first)] != negative(v.first) &&
           c[var(v.second)] != negative(v.second);
}

inline bool is_in_incomplete(const SharedDBPropagator& ic, Vertex v) {
    return ic.is_true(v.first) && ic.is_true(v.second);
}

} // namespace detail

void GurobiCliqueSolverG2::p_create_initial_vars() {
    const auto n = m_subgraph->n();
    m_vars.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        m_vars.emplace_back(m_model.addVar(0.0, 1.0, 1.0, GRB_CONTINUOUS));
    }
}

template <typename ConfigurationIterator>
inline void
GurobiCliqueSolverG2::add_configurations(ConfigurationIterator cbeg,
                                         ConfigurationIterator cend) {
    using ConfigType =
        typename std::iterator_traits<ConfigurationIterator>::value_type;
    std::for_each(cbeg, cend,
                  [this](const ConfigType& cfg) { add_configuration(cfg); });
}

template <typename ConfigurationType>
inline void
GurobiCliqueSolverG2::add_configuration(const ConfigurationType& cfg) {
    auto literal_value = [&cfg](Lit l) {
        return cfg[lit::var(l)] != lit::negative(l);
    };
    m_have_new_constraints = true;
    m_complete_independent_dual.reset();
    m_complete_independent_sets.emplace_back(cfg);
    m_complete_independent_set_zero_since.push_back(0);
    const auto n = m_subgraph->n();
    m_expr_buffer.clear();
    for (std::size_t vi = 0; vi < n; ++vi) {
        Vertex v = m_subgraph->vertex(vi);
        if (literal_value(v.first) && literal_value(v.second)) {
            m_expr_buffer += m_vars[vi];
        }
    }
    m_complete_independent_set_constraints.push_back(
        m_model.addConstr(m_expr_buffer <= 1));
}

template <typename Callable>
inline void GurobiCliqueSolverG2::p_foreach_literal(Callable&& callable) {
    for (Lit l = 0, nall_lit = 2 * m_nall; l < nall_lit; ++l) {
        std::invoke(std::forward<Callable>(callable), l);
    }
}

SolverState GurobiCliqueSolverG2::solve_relaxation(double time_allowed) {
    if (m_subgraph_upper_bound <= m_best_mes.size() && !m_have_new_constraints)
    {
        return SolverState::OPTIMUM_ON_SUBGRAPH;
    }
    if (time_allowed <= 0.0) {
        return SolverState::TIMEOUT_NO_IMPROVEMENT;
    }
    if (std::isfinite(time_allowed)) {
        m_model.set(GRB_DoubleParam_TimeLimit, time_allowed);
    }
    // m_recorder->store_event("BEGIN_SOLVE_RELAXATION");
    m_model.optimize();
    ++m_num_solve_calls;
    m_have_new_constraints = false;
    if (m_model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
        m_recorder->store_event("TIMEOUT_SOLVE_RELAXATION");
        return SolverState::TIMEOUT_NO_IMPROVEMENT;
    }
    std::size_t prev_clique_size = m_best_mes.size();
    p_extract_ordered_solution();
    // m_recorder->store_event("DONE_SOLVE_RELAXATION", {{"bound",
    // m_last_objective_value}}, "bound");
    p_ordered_solution_to_clique();
    std::size_t new_clique_size = m_best_mes.size();
    bool is_now_opt = (new_clique_size >= m_subgraph_upper_bound);
    return is_now_opt ? SolverState::OPTIMUM_ON_SUBGRAPH
                      : (new_clique_size > prev_clique_size
                             ? SolverState::IMPROVEMENT_FOUND
                             : SolverState::NO_IMPROVEMENT_FOUND);
}

void GurobiCliqueSolverG2::p_extract_ordered_solution(
    bool possible_suboptimal_mip) {
    m_ordered_solution.clear();
    m_last_solution_values.reset(
        m_model.get(GRB_DoubleAttr_X, m_vars.data(), int(m_vars.size())));
    const auto n = m_subgraph->n();
    for (std::size_t i = 0; i < n; ++i) {
        if (m_last_solution_values[i] > 1.0e-6) {
            m_ordered_solution.emplace_back(i, m_last_solution_values[i]);
        }
    }
    std::sort(
        m_ordered_solution.begin(), m_ordered_solution.end(),
        [](const auto& p1, const auto& p2) { return p1.second > p2.second; });
    if (possible_suboptimal_mip) {
        m_last_objective_value = m_model.get(GRB_DoubleAttr_ObjBound);
        m_subgraph_upper_bound =
            p_make_objective_value_integral(m_last_objective_value);
    } else {
        m_last_objective_value = m_model.get(GRB_DoubleAttr_ObjVal);
        m_subgraph_upper_bound =
            p_make_objective_value_integral(m_last_objective_value);
        p_update_dual_values();
    }
}

std::size_t
GurobiCliqueSolverG2::p_make_objective_value_integral(double v) const noexcept {
    auto i = std::make_signed_t<std::size_t>(std::floor(v + 0.01));
    if (i < 0)
        return 0;
    return std::size_t(i);
}

void GurobiCliqueSolverG2::p_ordered_solution_to_clique() {
    if (m_subgraph_upper_bound > m_best_mes.size()) {
        m_clique_builder.reset_vertices();
        m_clique_builder.greedily_extend_weighted(m_ordered_solution.begin(),
                                                  m_ordered_solution.end());
        m_clique_builder.randomly_make_maximal(sammy::rng());
        if (m_clique_builder.size() > m_best_mes.size()) {
            m_recorder->store_event(
                "IMPROVED_LB",
                {{"lb", m_clique_builder.size()}, {"method", "LP rounding"}},
                "lb", "method");
            m_best_mes = m_clique_builder.get_vertices();
        }
    }
}

SolverState GurobiCliqueSolverG2::solve_full_relaxation(double time_allowed) {
    if (m_subgraph_upper_bound == m_best_mes.size() && !m_have_new_constraints)
    {
        return SolverState::OPTIMUM_ON_SUBGRAPH;
    }
    if (time_allowed <= 0.0) {
        return SolverState::TIMEOUT_NO_IMPROVEMENT;
    }
    m_recorder->store_event(
        "BEGIN_SOLVE_FULL_RELAXATION",
        {{"incomplete_sets", m_inc_independent_sets.size()},
         {"complete_sets", m_complete_independent_sets.size()},
         {"n", m_subgraph->n()}},
        "incomplete_sets", "complete_sets", "n");
    auto before = Clock::now();
    bool any_improvement_found = false;
    for (;;) {
        double time_remaining =
            time_allowed - seconds_between(before, Clock::now());
        switch (solve_relaxation(time_remaining)) {
        case SolverState::IMPROVEMENT_FOUND:
            any_improvement_found = true;
            break;
        case SolverState::NO_IMPROVEMENT_FOUND:
            break;
        case SolverState::OPTIMUM_ON_SUBGRAPH:
            m_recorder->store_event("DONE_SOLVE_FULL_RELAXATION",
                                    {{"bound", m_subgraph_upper_bound},
                                     {"status", "optimum on subgraph"},
                                     {"best_mes", m_best_mes.size()}},
                                    "bound", "status", "best_mes");
            return SolverState::OPTIMUM_ON_SUBGRAPH;
        case SolverState::TIMEOUT_NO_IMPROVEMENT:
            m_recorder->store_event("TIMEOUT_SOLVE_FULL_RELAXATION");
            return any_improvement_found ? SolverState::TIMEOUT_IMPROVEMENT
                                         : SolverState::TIMEOUT_NO_IMPROVEMENT;
        case SolverState::TIMEOUT_IMPROVEMENT:
            m_recorder->store_event("TIMEOUT_SOLVE_FULL_RELAXATION");
            return SolverState::TIMEOUT_IMPROVEMENT;
        }
        if (!prohibit_violated_nonedges()) {
            m_recorder->store_event("DONE_SOLVE_FULL_RELAXATION",
                                    {{"bound", m_subgraph_upper_bound},
                                     {"best_mes", m_best_mes.size()}},
                                    "bound", "best_mes");
            return any_improvement_found ? SolverState::IMPROVEMENT_FOUND
                                         : SolverState::NO_IMPROVEMENT_FOUND;
        }
    }
}

/**
 * Identify violated non-edges.
 */
bool GurobiCliqueSolverG2::p_identify_violated_nonedges() {
    m_violated_nonedges.clear();
    for (OSolIter c = m_ordered_solution.begin(), e = m_ordered_solution.end();
         c != e; ++c)
    {
        double weight = c->second;
        std::size_t index = c->first;
        double thresh = 1.01 - weight;
        if (thresh > weight)
            break;
        const auto& row_c = m_subgraph->matrix_row(index);
        for (OSolIter d = std::next(c); d != e; ++d) {
            double weight_d = d->second;
            if (weight_d < thresh)
                break;
            std::size_t index_d = d->first;
            if (!row_c[index_d])
                m_violated_nonedges.emplace_back(c, d);
        }
    }
    return !m_violated_nonedges.empty();
}

bool GurobiCliqueSolverG2::p_prohibit_extend_existing(
    std::size_t v_index, DynamicBitset& covered_vne) {
    auto handle_other = [&](SharedDBPropagator& target_prop, std::size_t vio) {
        if (covered_vne[vio])
            return;
        auto [vit2, wit2] = m_violated_nonedges[vio];
        Vertex v = m_subgraph->vertex(vit2->first),
               w = m_subgraph->vertex(wit2->first);
        if (push_vertex_pair(target_prop, v, w)) {
            covered_vne[vio] = true;
        }
    };
    auto [vit, wit] = m_violated_nonedges[v_index];
    for (std::size_t ind_index = 0, n_inc = m_inc_independent_sets.size();
         ind_index < n_inc; ++ind_index)
    {
        std::size_t vi = vit->first, wi = wit->first;
        Vertex v = m_subgraph->vertex(vi), w = m_subgraph->vertex(wi);
        SharedDBPropagator& target_prop = m_inc_independent_sets[ind_index];
        if (target_prop.is_false(v.first) || target_prop.is_false(v.second) ||
            target_prop.is_false(w.first) || target_prop.is_false(w.second))
        {
            continue;
        }
        if (push_vertex_pair(target_prop, v, w)) {
            covered_vne[v_index] = true;
            auto nvio = m_violated_nonedges.size();
            for (std::size_t neindex = v_index + 1; neindex < nvio; ++neindex) {
                handle_other(target_prop, neindex);
            }
            p_refresh_inc_constraint(ind_index);
            return true;
        }
    }
    return false;
}

bool GurobiCliqueSolverG2::prohibit_violated_nonedges() {
    if (!p_identify_violated_nonedges()) {
        return false;
    }
    DynamicBitset covered_vne(m_violated_nonedges.size(), false);
    OSolIter last_first = m_ordered_solution.end();
    const auto nve = m_violated_nonedges.size();
    SharedDBPropagator& prop = m_subgraph->get_propagator();
    auto handle_other_ve = [&](std::size_t j) {
        if (!covered_vne[j]) {
            auto [xit, yit] = m_violated_nonedges[j];
            Vertex x = m_subgraph->vertex(xit->first);
            Vertex y = m_subgraph->vertex(yit->first);
            if (push_vertex_pair(prop, x, y)) {
                covered_vne[j].set();
            }
        }
    };
    for (std::size_t i = 0; i < nve; ++i) {
        auto [vit, wit] = m_violated_nonedges[i];
        if (vit != last_first && !covered_vne[i]) {
            if (p_prohibit_extend_existing(i, covered_vne))
                continue;
            last_first = vit;
            Vertex v = m_subgraph->vertex(vit->first);
            Vertex w = m_subgraph->vertex(wit->first);
            reset_and_push_noresolve(prop, v);
            if (push_vertex(prop, w) < 0) {
                // can happen if we've learned new clauses;
                // essentially, this is a 'hidden' edge in our graph
                covered_vne[i].set();
                m_subgraph->nonedge_to_edge(vit->first, wit->first);
                continue;
            }
            covered_vne[i].set();
            for (std::size_t j = i + 1; j < nve; ++j) {
                handle_other_ve(j);
            }
            for (std::size_t j = 0; j < i; ++j) {
                handle_other_ve(j);
            }
            add_potentially_incomplete_independent_set(prop);
        }
    }
    m_have_new_constraints = true;
    return true;
}

void GurobiCliqueSolverG2::add_direct_independent_set(
    const std::vector<std::size_t>& indices) {
    SharedDBPropagator& propagator = m_subgraph->get_propagator();
    DynamicBitset compatible_literals(2 * m_nall, true);
    for (std::size_t index : indices) {
        propagator.reset_or_throw();
        Vertex v = m_subgraph->vertex(index);
        if (push_vertex(propagator, v) < 0) {
            throw std::logic_error("Broken vertex in independent set!");
        }
        for (Lit l : propagator.get_trail()) {
            compatible_literals[lit::negate(l)].reset();
        }
    }
    m_direct_independent_dual.reset();
    m_direct_independent_sets.push_back(
        {HashSet<std::size_t>(indices.begin(), indices.end()),
         std::move(compatible_literals)});
    m_expr_buffer.clear();
    for (std::size_t index : indices) {
        m_expr_buffer += m_vars[index];
    }
    m_direct_independent_set_constraints.push_back(
        m_model.addConstr(m_expr_buffer <= 1));
}

void GurobiCliqueSolverG2::add_potentially_incomplete_independent_set(
    const SharedDBPropagator& propagator) {
    if (propagator.get_trail().size() == m_nall) {
        std::vector<bool> new_cfg(m_nall, false);
        for (Lit l : propagator.get_trail()) {
            if (!lit::negative(l)) {
                new_cfg[lit::var(l)] = true;
            }
        }
        add_configuration(new_cfg);
    } else {
        m_inc_independent_dual.reset();
        m_inc_independent_sets.push_back(propagator);
        m_inc_independent_set_zero_since.push_back(0);
        m_expr_buffer.clear();
        for (std::size_t vi = 0, n = m_subgraph->n(); vi != n; ++vi) {
            Vertex v = m_subgraph->vertex(vi);
            if (propagator.is_true(v.first) && propagator.is_true(v.second)) {
                m_expr_buffer += m_vars[vi];
            }
        }
        m_inc_independent_set_constraints.push_back(
            m_model.addConstr(m_expr_buffer <= 1));
    }
}

void GurobiCliqueSolverG2::add_vertices(const std::vector<Vertex>& vertices) {
    if (vertices.empty())
        return;
    m_subgraph_upper_bound = m_best_sample.size();
    const auto old_size = m_subgraph->n();
    m_subgraph->add_vertices(vertices);
    const auto new_size = m_subgraph->n();
    if (old_size + vertices.size() != new_size) {
        throw std::logic_error(
            "Duplicate or already existing vertices passed to add_vertices!");
    }
    for (Vertex v : vertices) {
        GRBColumn column;
        for (std::size_t i = 0, ni = m_complete_independent_sets.size();
             i != ni; ++i)
        {
            if (detail::is_in_complete(m_complete_independent_sets[i], v)) {
                column.addTerm(1.0, m_complete_independent_set_constraints[i]);
            }
        }
        for (std::size_t i = 0, ni = m_inc_independent_sets.size(); i != ni;
             ++i)
        {
            if (detail::is_in_incomplete(m_inc_independent_sets[i], v)) {
                column.addTerm(1.0, m_inc_independent_set_constraints[i]);
            }
        }
        m_vars.push_back(m_model.addVar(0.0, 1.0, 1.0, GRB_CONTINUOUS, column));
    }
}

void GurobiCliqueSolverG2::p_categorize_good_vertices(
    std::vector<Vertex>& good, VertexMapTo<double>& interesting_scores,
    std::size_t num_categories, std::size_t goal_vertices) {
    if (goal_vertices >= good.size())
        return;
    std::vector<std::size_t> category_count(num_categories, 0);
    const double delta_cat = 1.0 / num_categories;
    auto category_of = [&](double value) -> std::size_t {
        if (value <= 0)
            return 0;
        if (value >= 1)
            return num_categories - 1;
        return std::size_t(value / delta_cat);
    };
    for (Vertex v : good) {
        double val = interesting_scores.at(v);
        ++category_count[category_of(val)];
    }
    if (category_count.front() > goal_vertices) {
        goal_vertices = (std::max)(std::size_t(0.2 * category_count.front()),
                                   goal_vertices);
    }
    std::size_t take_all_from_below = 0;
    std::size_t total_count = 0;
    double last_ratio = 1.0;
    for (std::size_t c = 0; c < num_categories; ++c) {
        std::size_t this_cat = category_count[c];
        std::size_t count_after = total_count + this_cat;
        if (count_after <= goal_vertices) {
            ++take_all_from_below;
            total_count = count_after;
        } else {
            last_ratio = double(goal_vertices - total_count) / this_cat;
            break;
        }
    }
    auto& rng = sammy::rng();
    good.erase(std::remove_if(
                   good.begin(), good.end(),
                   [&](Vertex v) {
                       std::size_t cat = category_of(interesting_scores.at(v));
                       if (cat < take_all_from_below)
                           return false;
                       if (cat > take_all_from_below)
                           return true;
                       return bool(std::uniform_real_distribution<double>{
                                       0.0, 1.0}(rng) > last_ratio);
                   }),
               good.end());
}

template <typename VertexIterator>
std::size_t GurobiCliqueSolverG2::price_vertices(VertexIterator vbeg,
                                                 VertexIterator vend) {
    std::vector<Vertex> good_vertices, possible_vertices;
    VertexMapTo<double> interesting_scores;
    if (!p_collect_good_or_possible(vbeg, vend, good_vertices,
                                    possible_vertices, interesting_scores))
    {
        return 0;
    }
    std::size_t ret_result = good_vertices.size() + possible_vertices.size();
    if (!good_vertices.empty()) {
        std::size_t total_good = good_vertices.size();
        p_categorize_good_vertices(good_vertices, interesting_scores, 20,
                                   std::size_t(0.2 * m_subgraph->n()));
        std::size_t total_taken = good_vertices.size();
        p_add_vertices_to_independent_sets(good_vertices, interesting_scores);
        std::size_t added_to_independent_sets =
            total_taken - good_vertices.size();
        m_recorder->store_event(
            "ADDED_VERTICES",
            {{"added_to_graph", good_vertices.size()},
             {"added_to_indsets", added_to_independent_sets},
             {"considered", total_taken},
             {"total_good", total_good},
             {"type", "good"}},
            "added_to_graph", "added_to_indsets", "considered", "total_good");
        add_vertices(good_vertices);
    } else {
        std::size_t analyzed = possible_vertices.size();
        p_add_vertices_to_independent_sets(possible_vertices,
                                           interesting_scores);
        m_recorder->store_event("ADDED_VERTICES",
                                {{"new_count", possible_vertices.size()},
                                 {"total_count", analyzed},
                                 {"type", "possible"}},
                                "new_count", "total_count", "type");
        add_vertices(possible_vertices);
    }
    return ret_result;
}

void GurobiCliqueSolverG2::p_update_dual_values() {
    m_complete_independent_dual.reset(m_model.get(
        GRB_DoubleAttr_Pi, m_complete_independent_set_constraints.data(),
        int(m_complete_independent_set_constraints.size())));
    m_inc_independent_dual.reset(
        m_model.get(GRB_DoubleAttr_Pi, m_inc_independent_set_constraints.data(),
                    int(m_inc_independent_set_constraints.size())));
    m_direct_independent_dual.reset(m_model.get(
        GRB_DoubleAttr_Pi, m_direct_independent_set_constraints.data(),
        int(m_direct_independent_set_constraints.size())));
    for (std::size_t i = 0; i < m_complete_independent_sets.size(); ++i) {
        if (m_complete_independent_dual[i] > 1.0e-5) {
            m_complete_independent_set_zero_since[i] = 0;
        } else {
            ++m_complete_independent_set_zero_since[i];
        }
    }
    for (std::size_t i = 0; i < m_inc_independent_sets.size(); ++i) {
        if (m_inc_independent_dual[i] > 1.0e-5) {
            m_inc_independent_set_zero_since[i] = 0;
        } else {
            ++m_inc_independent_set_zero_since[i];
        }
    }
}

double GurobiCliqueSolverG2::p_price_in_existing(Vertex v) const noexcept {
    double score = 0.0;
    for (std::size_t i = 0, ni = m_complete_independent_sets.size(); i < ni;
         ++i)
    {
        if (detail::is_in_complete(m_complete_independent_sets[i], v)) {
            score += m_complete_independent_dual[i];
            if (score >= 1.02)
                return 1.02;
        }
    }
    for (std::size_t i = 0, ni = m_inc_independent_sets.size(); i < ni; ++i) {
        if (detail::is_in_incomplete(m_inc_independent_sets[i], v)) {
            score += m_inc_independent_dual[i];
            if (score >= 1.02)
                return 1.02;
        }
    }
    return score;
}

void GurobiCliqueSolverG2::p_add_vertices_to_independent_sets(
    std::vector<Vertex>& vertices, VertexMapTo<double>& vertex_scores) {
    for (Vertex v : vertices) {
        double& weight = vertex_scores.at(v);
        for (std::size_t i = 0, ni = m_inc_independent_sets.size(); i != ni;
             ++i)
        {
            SharedDBPropagator& prop = m_inc_independent_sets[i];
            if (!detail::is_in_incomplete(prop, v) && push_vertex(prop, v) >= 0)
            {
                weight += m_inc_independent_dual[i];
                if (weight >= 1.01)
                    break;
            }
        }
    }
    auto has_high_weight = [&vertex_scores](Vertex v) {
        return vertex_scores.at(v) >= 1.01;
    };
    vertices.erase(
        std::remove_if(vertices.begin(), vertices.end(), has_high_weight),
        vertices.end());
}

bool GurobiCliqueSolverG2::greedy_add_to_cutting_planes() {
    std::size_t num_refreshed = 0;
    std::size_t num_strengthened = 0;
    for (std::size_t i = 0, ni = m_inc_independent_sets.size(); i != ni; ++i) {
        double could_add = 0.0, current_weight = 0.0;
        SharedDBPropagator& prop = m_inc_independent_sets[i];
        for (const auto& e : m_ordered_solution) {
            std::size_t vi = e.first;
            double wi = e.second;
            Vertex v = m_subgraph->vertex(vi);
            if (detail::is_in_incomplete(prop, v)) {
                current_weight += wi;
            } else if (can_push(prop, v)) {
                could_add += wi;
            }
        }
        if (current_weight > 1.01) {
            // we have, possibly indirectly, added
            // present vertices to this independent set
            p_refresh_inc_constraint(i);
            ++num_strengthened;
            ++num_refreshed;
        } else if (current_weight + could_add > 1.01) {
            for (const auto& e : m_ordered_solution) {
                Vertex v = m_subgraph->vertex(e.first);
                if (detail::is_in_incomplete(prop, v)) {
                    continue;
                }
                if (push_vertex(prop, v) >= 0) {
                    current_weight += e.second;
                    m_model.chgCoeff(m_inc_independent_set_constraints[i],
                                     m_vars[e.first], 1.0);
                }
            }
            if (current_weight > 1.01) {
                ++num_strengthened;
            }
        }
    }
    if (num_strengthened == 0) {
        num_strengthened = p_greedy_add_to_direct();
    }
    if (num_strengthened > 0) {
        m_have_new_constraints = true;
        m_recorder->store_event("STRENGTHENED_CONSTRAINTS",
                                {{"num_strengthened", num_strengthened},
                                 {"num_refreshed", num_refreshed}},
                                "num_strengthened", "num_refreshed");
    }
    return num_strengthened > 0;
}

void GurobiCliqueSolverG2::p_refresh_inc_constraint(std::size_t inc_index) {
    m_expr_buffer.clear();
    const SharedDBPropagator& prop = m_inc_independent_sets[inc_index];
    for (std::size_t vi = 0, n = m_subgraph->n(); vi != n; ++vi) {
        Vertex v = m_subgraph->vertex(vi);
        if (detail::is_in_incomplete(prop, v)) {
            m_expr_buffer += m_vars[vi];
        }
    }
    GRBConstr constraint = m_model.addConstr(m_expr_buffer <= 1);
    assert(inc_index < m_inc_independent_set_constraints.size());
    m_model.remove(m_inc_independent_set_constraints[inc_index]);
    m_inc_independent_set_constraints[inc_index] = constraint;
}

std::size_t GurobiCliqueSolverG2::p_greedy_add_to_direct() {
    std::size_t num_strengthened = 0;
    auto& propagator = m_subgraph->get_propagator();
    for (std::size_t direct_index = 0, dn = m_direct_independent_sets.size();
         direct_index < dn; ++direct_index)
    {
        auto& direct = m_direct_independent_sets[direct_index];
        double current_weight = std::accumulate(
            m_ordered_solution.begin(), m_ordered_solution.end(), 0.0,
            [&](double acc, const auto& e) {
                if (direct.vertex_indices.count(e.first)) {
                    return acc + e.second;
                }
                return acc;
            });
        for (auto [index, value] : m_ordered_solution) {
            if (direct.vertex_indices.count(index))
                continue;
            Vertex v = m_subgraph->vertex(index);
            if (!direct.compatible_literals[v.first] ||
                !direct.compatible_literals[v.second])
            {
                continue;
            }
            if (!p_direct_is_compatible(direct, index)) {
                continue;
            }
            current_weight += value;
            direct.vertex_indices.insert(index);
            reset_and_push_noresolve(propagator, v);
            for (Lit l : propagator.get_trail()) {
                direct.compatible_literals[lit::negate(l)].reset();
            }
            m_model.chgCoeff(m_direct_independent_set_constraints[direct_index],
                             m_vars[index], 1.0);
        }
        if (current_weight > 1.01) {
            ++num_strengthened;
        }
    }
    return num_strengthened;
}

bool GurobiCliqueSolverG2::p_greedy_generate_direct_independent_sets() {
    if (m_ordered_solution.size() < 6) {
        return false;
    }
    auto ordered_solution_copy = m_ordered_solution;
    std::size_t initial_part_size = m_ordered_solution.size() / 3;
    auto initial_part_begin = ordered_solution_copy.begin();
    auto second_part_begin = initial_part_begin + initial_part_size;
    auto second_part_end = ordered_solution_copy.end();
    auto& rng = sammy::rng();
    double largest_weight = 0.0;
    std::vector<std::size_t> largest_set;
    m_indset_builder.reset_vertices();
    largest_weight = m_indset_builder.greedily_extend_weighted(
        initial_part_begin, second_part_end);
    largest_set = m_indset_builder.get_indices();
    for (std::size_t it = 0; it < 100; ++it) {
        std::shuffle(initial_part_begin, second_part_begin, rng);
        std::shuffle(second_part_begin, second_part_end, rng);
        m_indset_builder.reset_vertices();
        double weight = m_indset_builder.greedily_extend_weighted(
            initial_part_begin, second_part_end);
        weight += m_indset_builder.greedily_extend_weighted(second_part_begin,
                                                            second_part_end);
        if (weight > 1.01) {
            if (weight > largest_weight) {
                largest_weight = weight;
                largest_set = m_indset_builder.get_indices();
            }
        }
    }
    if (largest_weight > 1.01) {
        add_direct_independent_set(largest_set);
        return true;
    }
    // if nothing is found, return false
    return false;
}

bool GurobiCliqueSolverG2::greedy_generate_cutting_planes() {
    using IEntry = std::pair<std::vector<Vertex>, double>;
    SharedDBPropagator& propagator = m_subgraph->get_propagator();
    std::vector<IEntry> violated;
    std::vector<Vertex> current_pushed;
    double added_weights = 0.0;
    auto handle_entry = [&](std::pair<std::size_t, double> e) {
        Vertex v = m_subgraph->vertex(e.first);
        if (push_vertex(propagator, v) >= 0) {
            added_weights += e.second;
            current_pushed.push_back(v);
        }
    };
    for (OSolIter cbegin = m_ordered_solution.begin(),
                  cend = m_ordered_solution.end();
         cbegin != cend; ++cbegin)
    {
        added_weights = 0.0;
        current_pushed.clear();
        propagator.reset_or_throw();
        for (const auto& e : IteratorRange{cbegin, cend})
            handle_entry(e);
        for (const auto& e : IteratorRange{m_ordered_solution.cbegin(), cbegin})
            handle_entry(e);
        if (added_weights >= 1.01) {
            violated.emplace_back(current_pushed, added_weights);
        }
    }
    if (violated.empty()) {
        return p_greedy_generate_direct_independent_sets();
    }
    auto pos = violated.end();
    if (violated.size() > 10) {
        pos = violated.begin() + 10;
        std::nth_element(violated.begin(), pos, violated.end(),
                         [](const IEntry& e1, const IEntry& e2) {
                             return e1.second > e2.second;
                         });
    }
    double max_vio = 1.0;
    for (auto i = violated.begin(); i != pos; ++i) {
        if (i->second > max_vio)
            max_vio = i->second;
        propagator.reset_or_throw();
        for (Vertex v : i->first) {
            push_vertex(propagator, v);
        }
        add_potentially_incomplete_independent_set(propagator);
    }
    m_recorder->store_event(
        "GREEDY_CUTTING_PLANES",
        {{"count", pos - violated.begin()}, {"max_violation", max_vio - 1.0}},
        "count", "max_violation");
    m_have_new_constraints = true;
    return true;
}

template <typename VertexIterator>
bool GurobiCliqueSolverG2::p_collect_good_or_possible(
    VertexIterator vbeg, VertexIterator vend,
    std::vector<Vertex>& good_vertices, std::vector<Vertex>& possible_vertices,
    VertexMapTo<double>& interesting_scores) {
    assert(m_inc_independent_sets.empty() ||
           m_inc_independent_dual.get() != nullptr);
    assert(m_complete_independent_sets.empty() ||
           m_complete_independent_dual.get() != nullptr);
    struct Context {
        bool have_good = false;
        bool many_good = false;
    };
    const auto n = m_subgraph->n();
    thread_group().parallel_foreach(
        vbeg, vend, Context{}, [&](Context& ctx, Vertex v) {
            // if(ctx.many_good) return;
            if (m_subgraph->has_vertex(v))
                return;
            double ex_score = p_price_in_existing(v);
            if (ex_score >= 1.01)
                return;
            if (ex_score >= 0.99 && ctx.have_good)
                return;
            std::unique_lock l{ds_lock};
            ctx.have_good = !good_vertices.empty();
            ctx.many_good = (good_vertices.size() > 0.2 * n);
            // if(ctx.many_good) return;
            if (ex_score >= 0.99 && ctx.have_good)
                return;
            if (ex_score < 0.99) {
                ctx.have_good = true;
                if (!possible_vertices.empty()) {
                    possible_vertices.clear();
                    interesting_scores.clear();
                }
                good_vertices.push_back(v);
                interesting_scores[v] = ex_score;
            } else {
                possible_vertices.push_back(v);
                interesting_scores[v] = ex_score;
            }
        });
    return !good_vertices.empty() || !possible_vertices.empty();
}

void GurobiCliqueSolverG2::export_highest_weights_to_primal(
    ColoringHeuristicSolver& primal, std::size_t num_constraints) {
    struct ConstraintInfo {
        ConstraintInfo(std::size_t i, bool c, double d) noexcept
            : index(i), is_complete(c), dual_weight(d) {}

        std::size_t index;
        bool is_complete;
        double dual_weight;

        bool operator<(const ConstraintInfo& cinfo) const noexcept {
            return dual_weight > cinfo.dual_weight;
        }
    };

    std::vector<ConstraintInfo> cinfo;
    const auto nc = m_complete_independent_sets.size(),
               ni = m_inc_independent_sets.size();
    cinfo.reserve(nc + ni);
    for (std::size_t ic = 0; ic < nc; ++ic) {
        if (m_complete_independent_dual[ic] > 0.01) {
            cinfo.emplace_back(ic, true, m_complete_independent_dual[ic]);
        }
    }
    for (std::size_t ii = 0; ii < ni; ++ii) {
        if (m_inc_independent_dual[ii] > 0.01) {
            cinfo.emplace_back(ii, false, m_inc_independent_dual[ii]);
        }
    }
    auto end = cinfo.end();
    if (num_constraints < cinfo.size()) {
        end = cinfo.begin() + num_constraints;
        std::nth_element(cinfo.begin(), end, cinfo.end());
    }
    for (ConstraintInfo& ci : IteratorRange{cinfo.begin(), end}) {
        if (ci.is_complete) {
            p_complete_to_primal(primal, m_complete_independent_sets[ci.index]);
        } else {
            p_incomplete_to_primal(primal, m_inc_independent_sets[ci.index]);
        }
    }
}

void GurobiCliqueSolverG2::p_complete_to_primal(
    ColoringHeuristicSolver& primal, const DynamicBitset& complete_indset) {
    SharedDBPropagator propagator = m_subgraph->get_propagator();
    propagator.reset_or_throw();
    for (Var v = 0; v < complete_indset.size(); ++v) {
        Lit l =
            complete_indset[v] ? lit::positive_lit(v) : lit::negative_lit(v);
        if (propagator.is_false(l) ||
            (!propagator.is_true(l) && !propagator.push_level(l)))
        {
            throw std::logic_error(
                "Failed to push literal from complete configuration!");
        }
    }
    primal.add_color_class(std::move(propagator));
}

void GurobiCliqueSolverG2::p_incomplete_to_primal(
    ColoringHeuristicSolver& primal, const SharedDBPropagator& propagator) {
    primal.add_color_class(propagator);
}

std::vector<Vertex> GurobiCliqueSolverG2::get_fractional_mes_support() const {
    std::vector<Vertex> result;
    for (const auto& e : m_ordered_solution) {
        if (e.second > 0.01) {
            result.push_back(m_subgraph->vertex(e.first));
        } else
            break;
    }
    return result;
}

std::size_t GurobiCliqueSolverG2::cuts_from_primal(
    const std::vector<SharedDBPropagator>& primal_configs,
    std::size_t max_num_cuts) {
    using CVElement = std::pair<std::size_t, double>;
    std::vector<CVElement> cut_values;
    for (std::size_t i = 0, ncfg = primal_configs.size(); i < ncfg; ++i) {
        const SharedDBPropagator& prop = primal_configs[i];
        double total = 0.0;
        for (const auto& e : m_ordered_solution) {
            Vertex v = m_subgraph->vertex(e.first);
            if (prop.is_true(v.first) && prop.is_true(v.second)) {
                total += e.second;
            }
        }
        if (total > 1.01) {
            cut_values.emplace_back(i, total);
        }
    }
    if (cut_values.empty())
        return 0;
    auto end = cut_values.end();
    if (cut_values.size() > max_num_cuts) {
        end = cut_values.begin() + max_num_cuts;
        std::nth_element(cut_values.begin(), end, cut_values.end(),
                         [&](const CVElement& p1, const CVElement& p2) {
                             return p1.second > p2.second;
                         });
    }
    for (const CVElement& e : IteratorRange{cut_values.begin(), end}) {
        add_potentially_incomplete_independent_set(primal_configs[e.first]);
    }
    m_have_new_constraints = true;
    return std::size_t(end - cut_values.begin());
}

SolverState GurobiCliqueSolverG2::solve_as_mip(double time_allowed) {
    if (time_allowed <= 0)
        return SolverState::TIMEOUT_NO_IMPROVEMENT;
    if (std::isfinite(time_allowed)) {
        m_model.set(GRB_DoubleParam_TimeLimit, time_allowed);
    }
    p_begin_mip();
    m_model.optimize();
    m_num_solve_calls++;
    auto status = m_model.get(GRB_IntAttr_Status);
    auto solcount = m_model.get(GRB_IntAttr_SolCount);
    if (solcount > 0) {
        std::size_t old_mes_size = get_best_mes().size();
        p_extract_ordered_solution(true);
        p_ordered_solution_to_clique();
        p_end_mip();
        if (status == GRB_OPTIMAL)
            return SolverState::OPTIMUM_ON_SUBGRAPH;
        return get_best_mes().size() > old_mes_size
                   ? SolverState::TIMEOUT_IMPROVEMENT
                   : SolverState::TIMEOUT_NO_IMPROVEMENT;
    }
    p_end_mip();
    return SolverState::TIMEOUT_NO_IMPROVEMENT;
}

void GurobiCliqueSolverG2::p_begin_mip() {
    class Callback : public GRBCallback {
      public:
        explicit Callback(GurobiCliqueSolverG2* that) noexcept : that(that) {}

        void callback() override {
            if (that->m_abortable && that->m_aborted.load()) {
                GRBCallback::abort();
                return;
            }
            if (where != GRB_CB_MIPSOL) {
                return;
            }
            std::unique_ptr<double[]> sol_values{this->getSolution(
                that->m_vars.data(), int(that->m_vars.size()))};
            sol_vertices.clear();
            for (std::size_t i = 0, n = that->m_vars.size(); i < n; ++i) {
                if (sol_values[i] > 0.5)
                    sol_vertices.push_back(i);
            }
            for (std::size_t vi = 0, vn = sol_vertices.size(); vi < vn; ++vi) {
                std::size_t v = sol_vertices[vi];
                for (std::size_t vj = 0; vj < vi; ++vj) {
                    std::size_t w = sol_vertices[vj];
                    if (!that->m_subgraph->is_edge(v, w)) {
                        this->addLazy(that->m_vars[v] + that->m_vars[w] <= 1);
                    }
                }
            }
        }

      private:
        GurobiCliqueSolverG2* that;
        std::vector<std::size_t> sol_vertices;
    };
    for (GRBVar v : m_vars) {
        v.set(GRB_CharAttr_VType, GRB_BINARY);
    }
    m_old_callback.reset(m_callback.release());
    m_callback = std::make_unique<Callback>(this);
    m_model.set(GRB_IntParam_LazyConstraints, 1);
    m_model.setCallback(m_callback.get());
}

void GurobiCliqueSolverG2::p_end_mip() {
    for (GRBVar v : m_vars) {
        v.set(GRB_CharAttr_VType, GRB_CONTINUOUS);
    }
    m_model.set(GRB_IntParam_LazyConstraints, 0);
    m_callback.reset(m_old_callback.release());
    m_model.setCallback(m_callback.get());
}

void GurobiCliqueSolverG2::external_improved_mes(
    std::vector<Vertex> external_mes) {
    if (external_mes.size() <= m_best_mes.size())
        return;
    add_new_vertices(external_mes);
    m_best_mes = std::move(external_mes);
}

std::size_t GurobiCliqueSolverG2::cuts_from_sample(
    const std::vector<std::vector<bool>>& sample) {
    auto literal_value = [&](const std::vector<bool>& config, Lit l) -> bool {
        bool v = config[lit::var(l)];
        return lit::negative(l) ^ v;
    };
    auto vertex_value = [&](const std::vector<bool>& config, Vertex v) {
        return literal_value(config, v.first) &&
               literal_value(config, v.second);
    };
    std::size_t result = 0;
    for (std::size_t i = 0, ncfg = sample.size(); i < ncfg; ++i) {
        const auto& config = sample[i];
        double total = 0.0;
        for (const auto& e : m_ordered_solution) {
            Vertex v = m_subgraph->vertex(e.first);
            if (vertex_value(config, v))
                total += e.second;
        }
        if (total > 1.01) {
            add_configuration(config);
            ++result;
        }
    }
    return result;
}

bool GurobiCliqueSolverG2::p_direct_is_compatible(
    const DirectIndependentSet& direct, std::size_t index) {
    return std::all_of(direct.vertex_indices.begin(),
                       direct.vertex_indices.end(),
                       [&](std::size_t other_index) {
                           return !m_subgraph->is_edge(index, other_index);
                       });
}

} // namespace sammy

#endif
==> ./simplify_datastructure.h <==
#ifndef SAMMY_SIMPLIFY_DATASTRUCTURE_H_INCLUDED_
#define SAMMY_SIMPLIFY_DATASTRUCTURE_H_INCLUDED_

#include "algorithm_ex.h"
#include "clause_db.h"
#include "eliminate_subsumed.h"
#include "literals.h"
#include <cassert>

namespace sammy {

struct SimplifiedInstance {
    std::vector<Var> new_to_old;
    ClauseDB formula;
    Var num_concrete;
};

class SimplifyDatastructure {
  public:
    explicit SimplifyDatastructure(const ClauseDB& formula, Var num_concrete)
        : m_variable_eliminated(formula.num_vars(), false),
          m_variable_concrete(formula.num_vars(), false),
          m_var_presence(formula.num_vars()) {
        std::fill_n(m_variable_concrete.begin(), num_concrete, true);
        for (Lit u : formula.unary_literals()) {
            m_clauses.emplace_back(1, u);
        }
        for (auto [l1, l2] : formula.binary_clauses()) {
            m_clauses.emplace_back(std::initializer_list<Lit>{l1, l2});
        }
        for (CRef c = 1, n = formula.literal_db_size(); c < n;
             c = formula.next_clause(c))
        {
            auto lits = formula.lits_of(c);
            m_clauses.emplace_back(lits.begin(), lits.end());
        }
        assert(!has_tautology_or_double_literal());
    }

    /**
     * A slowish (for debugging only) routine that checks if the given set
     * of literals is a (not necessarily strict) superset of any clause in the
     * formula.
     * @param set A set (must have .count(element) method) of literals.
     * @return true iff a clause of the formula subsumes the given set of
     * literals.
     */
    template <typename SetType> bool has_clause_subsuming(const SetType& set) {
        auto is_contained = [&](Lit l) -> bool { return set.count(l); };
        auto subsumes = [&](const SCVec& v) -> bool {
            return std::all_of(v.begin(), v.end(), is_contained);
        };
        return std::any_of(m_clauses.begin(), m_clauses.end(), subsumes);
    }

    /**
     * A debug method for checking that all clauses are non-tautological
     * and no clause contains several copies of the same literal.
     */
    bool has_tautology_or_double_literal() {
        for (const SCVec& cl : m_clauses) {
            m_var_presence.clear();
            for (Lit l : cl) {
                Var v = lit::var(l);
                if (!m_var_presence.check_insert(v)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * A debug method for finding duplicate binary clauses.
     */
    bool has_duplicate_binary_clause() {
        EdgeSet eset;
        for (const SCVec& cl : m_clauses) {
            if (cl.size() != 2)
                continue;
            std::pair<Lit, Lit> edge{std::min(cl[0], cl[1]),
                                     std::max(cl[0], cl[1])};
            if (!eset.insert(edge).second) {
                return true;
            }
        }
        return false;
    }

    bool is_eliminated(Var old_var) const noexcept {
        return m_variable_eliminated[old_var];
    }

    bool is_concrete(Var old_var) const noexcept {
        return m_variable_concrete[old_var];
    }

    void mark_eliminated(Var v) noexcept { m_variable_eliminated[v] = true; }

    void push_reconstruction_clause(const SCVec& c) {
        m_reconstruction_stack.emplace_back(c);
    }

    void push_reconstruction_clause(SCVec&& c) {
        m_reconstruction_stack.emplace_back(std::move(c));
    }

    const std::vector<SCVec>& get_reconstruction_stack() const noexcept {
        return m_reconstruction_stack;
    }

    void remove_eliminated_clauses() {
        auto clause_eliminated = [&](const SCVec& c) {
            return std::find_if(c.begin(), c.end(), [&](Lit l) {
                       return m_variable_eliminated[lit::var(l)];
                   }) != c.end();
        };
        m_clauses.erase(std::remove_if(m_clauses.begin(), m_clauses.end(),
                                       clause_eliminated),
                        m_clauses.end());
        eliminate_subsumed(m_clauses, original_num_vars());
    }

    bool apply_fixes_and_equalities(const std::vector<Lit>& old_to_new) {
        bool result = false;
        for (Var v = 0, nv = old_to_new.size() / 2; v < nv; ++v) {
            if (m_variable_eliminated[v])
                continue;
            Lit p = lit::positive_lit(v);
            Lit m = old_to_new[p];
            if (m == simplify::fixed_positive() ||
                m == simplify::fixed_negative())
            {
                p_apply_fix(v, m);
                result = true;
            } else if (m != p) {
                p_apply_remap(v, m);
                result = true;
            }
        }
        if (result) {
            p_fix_and_eq_rewrite_clauses(old_to_new);
        }
        return result;
    }

    ClauseDB extract_clause_db() {
        ClauseDB result{Var(m_variable_eliminated.size())};
        for (const SCVec& cl : m_clauses) {
            result.add_clause(cl.data(), cl.data() + cl.size());
        }
        return result;
    }

    void add_clauses(const ClauseDB& clause_db,
                     ClauseCounts old_clause_counts) {
        auto new_u = clause_db.unary_literals(
            old_clause_counts.unary_clause_end, clause_db.num_unaries());
        auto new_b = clause_db.binary_clauses(
            old_clause_counts.binary_clause_end, clause_db.num_binaries());
        for (Lit u : new_u) {
            m_clauses.emplace_back(1, u);
        }
        for (auto [l1, l2] : new_b) {
            m_clauses.emplace_back(std::initializer_list<Lit>{l1, l2});
        }
        for (CRef c = old_clause_counts.long_clause_end,
                  n = clause_db.literal_db_size();
             c < n; c = clause_db.next_clause(c))
        {
            auto lits = clause_db.lits_of(c);
            m_clauses.emplace_back(lits.begin(), lits.end());
        }
    }

    Var original_num_vars() const noexcept {
        return Var(m_variable_eliminated.size());
    }

    void sort_clauses() {
        auto sort_clause = [](SCVec& clause) {
            std::sort(clause.begin(), clause.end());
        };
        auto cmp_clauses = [](const SCVec& c1, const SCVec& c2) -> bool {
            return c1.size() < c2.size() ||
                   (c1.size() == c2.size() &&
                    std::lexicographical_compare(c1.begin(), c1.end(),
                                                 c2.begin(), c2.end()));
        };
        std::for_each(m_clauses.begin(), m_clauses.end(), sort_clause);
        std::sort(m_clauses.begin(), m_clauses.end(), cmp_clauses);
    }

    SimplifiedInstance compress() {
        Var old_nv = m_variable_eliminated.size();
        std::vector<Var> new_to_old;
        Var new_count = 0;
        std::vector<Lit> old_to_new(2 * old_nv, simplify::eliminated());
        for (Var v = 0; v < old_nv; ++v) {
            if (m_variable_eliminated[v])
                continue;
            if (m_variable_concrete[v]) {
                new_to_old.push_back(v);
                old_to_new[lit::positive_lit(v)] = lit::positive_lit(new_count);
                old_to_new[lit::negative_lit(v)] = lit::negative_lit(new_count);
                ++new_count;
            }
        }
        Var new_num_concrete = new_count;
        for (Var v = 0; v < old_nv; ++v) {
            if (m_variable_eliminated[v])
                continue;
            if (!m_variable_concrete[v]) {
                new_to_old.push_back(v);
                old_to_new[lit::positive_lit(v)] = lit::positive_lit(new_count);
                old_to_new[lit::negative_lit(v)] = lit::negative_lit(new_count);
                ++new_count;
            }
        }
        return SimplifiedInstance{std::move(new_to_old),
                                  p_compress_clauses(new_count, old_to_new),
                                  new_num_concrete};
    }

    std::vector<bool>
    reconstruct_solution(const SimplifiedInstance& simp,
                         const std::vector<bool>& simp_sol) const {
        Var new_num_var = simp_sol.size();
        Var old_num_var = m_variable_eliminated.size();
        std::vector<bool> result(old_num_var, false);
        for (Var v = 0; v < new_num_var; ++v) {
            result[simp.new_to_old[v]] = simp_sol[v];
        }
        auto is_satisfied = [&](Lit l) {
            return lit::negative(l) ^ result[lit::var(l)];
        };
        IteratorRange stack_order{m_reconstruction_stack.rbegin(),
                                  m_reconstruction_stack.rend()};
        for (const SCVec& rc : stack_order) {
            if (std::find_if(rc.begin(), rc.end(), is_satisfied) == rc.end()) {
                result[lit::var(rc[0])].flip();
            }
        }
        return result;
    }

    std::vector<Vertex>
    reconstruct_lb(const SimplifiedInstance& simp,
                   const std::vector<Vertex>& simp_lb) const {
        std::vector<Vertex> result;
        result.reserve(simp_lb.size());
        auto transform_literal = [&](Lit l) {
            Var ov = simp.new_to_old[lit::var(l)];
            return lit::negative(l) ? lit::negative_lit(ov)
                                    : lit::positive_lit(ov);
        };
        auto transform_vertex = [&](Vertex simp_v) {
            return Vertex{transform_literal(simp_v.first),
                          transform_literal(simp_v.second)};
        };
        std::transform(simp_lb.begin(), simp_lb.end(),
                       std::back_inserter(result), transform_vertex);
        return result;
    }

    std::vector<std::vector<bool>>
    reconstruct_sample(const SimplifiedInstance& simplified_inst,
                       const std::vector<std::vector<bool>>& sample) const {
        std::vector<std::vector<bool>> result;
        result.reserve(sample.size());
        auto reconstruct = [&](const std::vector<bool>& cfg) {
            return reconstruct_solution(simplified_inst, cfg);
        };
        std::transform(sample.begin(), sample.end(), std::back_inserter(result),
                       reconstruct);
        return result;
    }

    std::vector<bool> reduce_solution(const SimplifiedInstance& simp,
                                      const std::vector<bool>& full_sol) const {
        Var nnew = simp.formula.num_vars();
        std::vector<bool> result(nnew, false);
        for (Var vnew = 0; vnew < nnew; ++vnew) {
            result[vnew] = full_sol[simp.new_to_old[vnew]];
        }
        return result;
    }

    std::vector<SCVec>& clauses() noexcept { return m_clauses; }

    const std::vector<SCVec>& clauses() const noexcept { return m_clauses; }

    void
    replace_all_binaries(const std::vector<std::pair<Lit, Lit>>& binaries) {
        auto is_long = [](const SCVec& cl) { return cl.size() > 2; };
        auto pair_to_scvec = [](const std::pair<Lit, Lit>& b) {
            return SCVec{std::initializer_list<Lit>{b.first, b.second}};
        };

        std::vector<SCVec> result;
        auto cnt_long =
            std::count_if(m_clauses.begin(), m_clauses.end(), is_long);
        result.reserve(binaries.size() + cnt_long);
        std::transform(binaries.begin(), binaries.end(),
                       std::back_inserter(result), pair_to_scvec);
        std::copy_if(m_clauses.begin(), m_clauses.end(),
                     std::back_inserter(result), is_long);
        m_clauses.swap(result);
    }

  private:
    ClauseDB p_compress_clauses(Var new_count,
                                const std::vector<Lit>& old_to_new) {
        ClauseDB result{new_count};
        CVec buffer;
        for (const SCVec& c : m_clauses) {
            buffer.clear();
            std::transform(c.begin(), c.end(), std::back_inserter(buffer),
                           [&](Lit l) { return old_to_new[l]; });
            assert(std::find(buffer.begin(), buffer.end(),
                             simplify::eliminated()) == buffer.end());
            result.add_clause(buffer.data(), buffer.data() + buffer.size());
        }
        return result;
    }

    void p_fix_and_eq_rewrite_clause(SCVec& c,
                                     const std::vector<Lit>& old_to_new) {
        auto transform_literal = [&](Lit l) { return old_to_new[l]; };
        auto not_fixed_negative = [&](Lit l) {
            return old_to_new[l] != simplify::fixed_negative();
        };
        auto fixed_positive = [&](Lit l) {
            return old_to_new[l] == simplify::fixed_positive();
        };
        if (std::any_of(c.begin(), c.end(), fixed_positive)) {
            c.clear();
            return;
        }
        c.erase(copy_transformed_if(c.begin(), c.end(), c.begin(),
                                    transform_literal, not_fixed_negative),
                c.end());
        if (c.empty())
            throw UNSATError();
        std::sort(c.begin(), c.end());
        c.erase(std::unique(c.begin(), c.end()), c.end());
        auto taut_pair = find_pair_if(c.begin(), c.end(), [](Lit l1, Lit l2) {
            return l1 == lit::negate(l2);
        });
        if (taut_pair != c.end()) {
            c.clear();
        }
    }

    void p_fix_and_eq_rewrite_clauses(const std::vector<Lit>& old_to_new) {
        std::for_each(m_clauses.begin(), m_clauses.end(), [&](SCVec& c) {
            p_fix_and_eq_rewrite_clause(c, old_to_new);
        });
        m_clauses.erase(
            std::remove_if(m_clauses.begin(), m_clauses.end(),
                           [](const SCVec& v) { return v.empty(); }),
            m_clauses.end());
    }

    void p_apply_fix(Var v, Lit l) {
        m_variable_eliminated[v] = true;
        Lit slit =
            (lit::negative(l) ? lit::negative_lit(v) : lit::positive_lit(v));
        m_reconstruction_stack.emplace_back(1, slit);
    }

    void p_apply_remap(Var v, Lit m) {
        if (m_variable_concrete[v]) {
            m_variable_concrete[lit::var(m)] = true;
        }
        Lit p = lit::positive_lit(v);
        Lit n = lit::negative_lit(v);
        Lit mneg = lit::negate(m);
        m_variable_eliminated[v] = true;
        m_reconstruction_stack.emplace_back(std::initializer_list<Lit>{n, m});
        m_reconstruction_stack.emplace_back(
            std::initializer_list<Lit>{p, mneg});
    }

    std::vector<SCVec> m_clauses;
    std::vector<SCVec> m_reconstruction_stack;
    std::vector<bool> m_variable_eliminated;
    std::vector<bool> m_variable_concrete;
    StampSet<Var, std::uint16_t> m_var_presence;
};

} // namespace sammy

#endif
==> ./barrage_lns_subproblem.h <==
#ifndef SAMMY_BARRAGE_LNS_SUBPROBLEM_H_INCLUDED_
#define SAMMY_BARRAGE_LNS_SUBPROBLEM_H_INCLUDED_

#include "dynamic_bitset.h"
#include "literals.h"
#include <memory>

namespace sammy {

/**
 * Structure for general representation of a subproblem
 * as encountered by our LNS.
 */
struct LNSSubproblem {
    std::vector<Vertex> uncovered_universe;
    std::vector<Vertex> mutually_exclusive_set;
    std::vector<DynamicBitset> removed_configurations;
    Lit num_concrete;
};

} // namespace sammy

#endif
==> ./barrage_worker_cnp.h <==
#ifndef SAMMY_BARRAGE_WORKER_CNP_H_INCLUDED_
#define SAMMY_BARRAGE_WORKER_CNP_H_INCLUDED_

#include "barrage.h"
#include "barrage_worker_with_core.h"
#include "gurobi_clique_solver_g2.h"
#include "universe_subgraph.h"

namespace sammy {

static constexpr std::size_t MAX_CUTS_FROM_PRIMAL = 5;
static constexpr double MIN_VERTICES_PER_PRICING = 10.0;
static constexpr double MIN_RELATIVE_PER_PRICING = 0.01;
static constexpr std::size_t MAX_CUT_ROUNDS_PER_PRICING = 40;
static constexpr double PRIMAL_TO_DUAL_GOAL_RATIO = 0.1;
static constexpr std::size_t UNIVERSE_SAMPLE_SIZE = 100'000;
static constexpr std::size_t SUBOPTIMAL_PRICE_SAMPLE_SIZE = 10'000;

/**
 * Cut-and-price portfolio core.
 * Runs a cut-and-price clique solver on a dynamic
 * subgraph of the full universe.
 *
 * Note: Core constructors are called and run in the worker thread.
 * They may raise InterruptError to signal that they detected the termination
 * flag being set. Does not have a need for a separate interruption/termination
 * flag, so does not have to update the already automatically set interruption
 * flag.
 */
class CutAndPricePortfolioCore {
  public:
    static std::unique_ptr<CutAndPricePortfolioCore>
    factory(PortfolioSolver* solver,
            PortfolioElementWithCore<CutAndPricePortfolioCore>* our_element) {
        return std::make_unique<CutAndPricePortfolioCore>(
            solver, our_element->get_mutable_recorder());
    }

    CutAndPricePortfolioCore(PortfolioSolver* solver,
                             EventRecorder* local_recorder)
        : m_solver(solver), m_clauses(solver->get_clauses()),
          m_single_thread(0), m_local_recorder(*local_recorder),
          m_subgraph(&m_clauses, &m_single_thread,
                     &solver->get_infeasibility_map(),
                     p_initial_vertex_set(solver)),
          m_base_solver(
              p_extend(m_subgraph), m_clauses.num_vars(), &m_local_recorder,
              solver->get_best_solution().assignments_as<DynamicBitset>(),
              solver->get_best_mes(), p_get_config()) {
        m_clauses.mark_frozen();
        m_base_solver.make_abortable();
    }

    void main() {
        try {
            // our interrupt flag is the element's termination flag
            while (!get_and_clear_interrupt_flag()) {
                if (p_check_global_optimality()) {
                    m_mes_optimal = true;
                    break;
                }
                if (!p_main_loop_should_continue()) {
                    break;
                }
            }
        } catch (const InterruptError&) {
            // timeout exception generated by some code
        }
        if (m_mes_optimal) {
            m_solver->report_global(
                "MES_OPTIMALITY_PROVED",
                {{"size", m_base_solver.get_best_mes().size()}}, "size");
        }
        m_local_recorder.store_event("CNP_WORKER_EXITING");
    }

    /**
     * Aside from handling the termination flag,
     * we do not need to check for other interrupt reasons.
     */
    void interrupt_if_necessary(const InterruptionCheckInfo& /*info*/) {}

    /**
     * Called by the PortfolioElementWithCore
     * if it finds the termination flag to be set.
     */
    void termination_flag_set() { m_base_solver.abort(); }

  private:
    bool p_main_loop_should_continue() {
        switch (m_base_solver.solve_full_relaxation()) {
        case SolverState::TIMEOUT_IMPROVEMENT:
        case SolverState::TIMEOUT_NO_IMPROVEMENT:
        default:
            return false;

        case SolverState::OPTIMUM_ON_SUBGRAPH:
            m_solver->report_mes(m_base_solver.get_best_mes(), "Cut & Price");
            m_solver->add_clique(m_base_solver.get_best_mes());
            if (p_subgraph_optimal()) {
                m_mes_optimal = true;
                return false;
            }
            return true;

        case SolverState::IMPROVEMENT_FOUND:
            m_solver->report_mes(m_base_solver.get_best_mes(), "Cut & Price");
            m_solver->add_clique(m_base_solver.get_best_mes());
            /* fall-through */
        case SolverState::NO_IMPROVEMENT_FOUND:
            if (!p_subgraph_non_optimal()) {
                m_solver->report_global("CNP_SEPARATION_FAILED", {});
                return false;
            }
            return true;
        }
    }

    bool p_check_global_optimality() {
        if (m_solver->get_best_mes_size() > m_base_solver.get_best_mes().size())
        {
            m_base_solver.external_improved_mes(m_solver->get_best_mes());
        }
        if (m_solver->get_best_solution_size() <
            m_base_solver.get_best_solution().size())
        {
            auto assignments = m_solver->get_best_solution()
                                   .assignments_as<std::vector<bool>>();
            m_base_solver.update_best_solution(assignments);
            m_base_solver.cuts_from_sample(assignments);
        }
        return m_base_solver.get_best_mes().size() >=
               m_base_solver.get_best_solution().size();
    }

    bool p_report_pricing(std::size_t added_vertices, bool sg_optimal = true) {
        if (added_vertices) {
            m_local_recorder.store_event("PRICING_ADDED_VERTICES",
                                         {{"num_added", added_vertices}},
                                         "num_added");
        } else {
            if (sg_optimal) {
                m_local_recorder.store_event("PRICING_FOUND_OPTIMALITY");
            } else {
                m_local_recorder.store_event("PRICING_FAILED");
            }
        }
        return added_vertices == 0;
    }

    bool p_subgraph_optimal() {
        m_cutting_planes_since_pricing = 0;
        double goal_vertices_added = MIN_RELATIVE_PER_PRICING * m_subgraph.n();
        goal_vertices_added =
            (std::max)(goal_vertices_added, MIN_VERTICES_PER_PRICING);
        std::size_t goal_added = std::size_t(std::ceil(goal_vertices_added));
        std::size_t added = 0;
        p_price_vertices(added, goal_added, m_solver->get_best_spawners()) ||
            p_price_vertices(added, goal_added, m_solver->get_all_spawners()) ||
            p_price_vertices(added, goal_added, m_solver->get_coloring_order());
        if (added > 0)
            return p_report_pricing(added);
        const auto& icache = m_solver->implied_cache();
        if (icache.have_reduced_universe()) {
            const auto& reduced = icache.get_reduced_universe();
            std::size_t usize = reduced.size();
            if (usize > 4 * UNIVERSE_SAMPLE_SIZE) {
                auto sample = sample_from_range(reduced, UNIVERSE_SAMPLE_SIZE,
                                                sammy::rng());
                p_price_vertices(added, goal_added, sample);
                if (added > 0)
                    return p_report_pricing(added);
            }
            p_price_vertices(added, goal_added, reduced);
        } else {
            std::size_t usize = m_solver->get_universe_size();
            auto& infmap = m_solver->get_infeasibility_map();
            if (usize > 4 * UNIVERSE_SAMPLE_SIZE) {
                auto sample =
                    infmap.sample_vertices(UNIVERSE_SAMPLE_SIZE, usize);
                p_price_vertices(added, goal_added, sample);
                if (added > 0)
                    return p_report_pricing(added);
            }
            auto all_vertices = infmap.collect_vertices(usize);
            p_price_vertices(added, goal_added, all_vertices);
        }
        return p_report_pricing(added);
    }

    bool p_subgraph_non_optimal() {
        if (++m_cutting_planes_since_pricing >= MAX_CUT_ROUNDS_PER_PRICING) {
            m_cutting_planes_since_pricing = 0;
            if (p_price_suboptimal())
                return true;
            return p_cut();
        }
        if (p_cut())
            return true;
        m_cutting_planes_since_pricing = 0;
        return p_price_suboptimal(true);
    }

    bool p_price_vertices(std::size_t& added, std::size_t goal_added,
                          const std::vector<Vertex>& vertices) {
        std::size_t pos_vertices =
            m_base_solver.price_vertices(vertices.begin(), vertices.end());
        added += pos_vertices;
        return added >= goal_added;
    }

    bool p_price_suboptimal(bool last_chance = false) {
        double goal_vertices_added = MIN_RELATIVE_PER_PRICING * m_subgraph.n();
        goal_vertices_added =
            (std::max)(goal_vertices_added, MIN_VERTICES_PER_PRICING);
        std::size_t goal_added = std::size_t(std::ceil(goal_vertices_added));
        std::size_t added = 0;
        p_price_vertices(added, goal_added, m_solver->get_best_spawners()) ||
            p_price_vertices(added, goal_added, m_solver->get_all_spawners()) ||
            p_price_vertices(added, goal_added, m_solver->get_coloring_order());
        if (added > 0)
            return p_report_pricing(added, false);
        const auto& icache = m_solver->implied_cache();
        if (icache.have_reduced_universe()) {
            const auto& reduced = icache.get_reduced_universe();
            std::size_t usize = reduced.size();
            if (usize > 4 * SUBOPTIMAL_PRICE_SAMPLE_SIZE) {
                auto sample = sample_from_range(
                    reduced, SUBOPTIMAL_PRICE_SAMPLE_SIZE, sammy::rng());
                p_price_vertices(added, goal_added, sample);
                if (added > 0)
                    return p_report_pricing(added, false);
            }
            if (!last_chance)
                return p_report_pricing(added, false);
            p_price_vertices(added, goal_added, reduced);
        } else {
            std::size_t usize = m_solver->get_universe_size();
            auto& infmap = m_solver->get_infeasibility_map();
            if (usize > 4 * SUBOPTIMAL_PRICE_SAMPLE_SIZE) {
                auto sample =
                    infmap.sample_vertices(SUBOPTIMAL_PRICE_SAMPLE_SIZE, usize);
                p_price_vertices(added, goal_added, sample);
                if (added > 0)
                    return p_report_pricing(added, false);
            }
            if (!last_chance)
                return p_report_pricing(added, false);
            auto all_vertices = infmap.collect_vertices(usize);
            p_price_vertices(added, goal_added, all_vertices);
        }
        return p_report_pricing(added, false);
    }

    bool p_cut() {
        if (m_base_solver.greedy_add_to_cutting_planes() ||
            m_base_solver.greedy_generate_cutting_planes())
        {
            return true;
        }
        return false;
    }

    static std::vector<Vertex> p_initial_vertex_set(PortfolioSolver* solver) {
        const auto& all_spawners = solver->get_all_spawners();
        const auto& best_spawners = solver->get_best_spawners();
        std::vector<Vertex> vertex_set =
            (best_spawners.size() < 2000 && all_spawners.size() < 20000)
                ? all_spawners
                : best_spawners;
        const auto& mes = solver->get_best_mes();
        vertex_set.insert(vertex_set.end(), mes.begin(), mes.end());
        std::sort(vertex_set.begin(), vertex_set.end());
        vertex_set.erase(std::unique(vertex_set.begin(), vertex_set.end()),
                         vertex_set.end());
        return vertex_set;
    }

    static UniverseSubgraph* p_extend(UniverseSubgraph& subgraph) {
        subgraph.extend_matrix_by_propagation();
        return &subgraph;
    }

    static LowerBoundMIPConfig p_get_config() {
        LowerBoundMIPConfig config;
        config.quiet_gurobi = true;
        return config;
    }

    PortfolioSolver* m_solver;
    ClauseDB& m_clauses;
    ThreadGroup<void> m_single_thread;
    EventRecorder& m_local_recorder;
    UniverseSubgraph m_subgraph;
    GurobiCliqueSolverG2 m_base_solver;
    bool m_mes_optimal = false;
    std::size_t m_cutting_planes_since_pricing = 0;
};

} // namespace sammy

#endif
==> ./lower_bound_mip_settings.h <==
#ifndef SAMMY_LOWER_BOUND_SOLVER_STATE_H_INCLUDED_
#define SAMMY_LOWER_BOUND_SOLVER_STATE_H_INCLUDED_

#include <cstddef>
#include <limits>

namespace sammy {

enum class SolverState {
    OPTIMUM_ON_SUBGRAPH,
    NO_IMPROVEMENT_FOUND,
    IMPROVEMENT_FOUND,
    TIMEOUT_NO_IMPROVEMENT,
    TIMEOUT_IMPROVEMENT
};

struct LowerBoundMIPConfig {
    bool solve_exactly_full_mip = false;
    bool solve_exactly_cut_and_price = false;
    bool quiet_gurobi = false;
    std::size_t stop_cuts_iteration_window = 5;
    std::size_t separation_rounds_before_pricing = 5;
    double stop_cuts_below_relative_improvement = 0.03;
    double total_mip_timeout = std::numeric_limits<double>::infinity();
};

struct SATColoringConfig {
    bool initial_sat_coloring = false;
    double initial_sat_coloring_timeout = 10.0;
};

} // namespace sammy

#endif
==> ./sat_lns.h <==
#ifndef SAMMY_SAT_LNS_H_INCLUDED_
#define SAMMY_SAT_LNS_H_INCLUDED_

#include "barrage.h"
#include "barrage_lns_subproblem.h"
#include "clause_db.h"
#include "shared_db_propagator.h"
#include <typeinfo>
#include <variant>

namespace sammy {

/**
 * Create a non-incremental SAT solver
 * that attempts to improve the current
 * solution by at least one configuration.
 */
template <typename BasicSatSolver> class FixedMESSATImprovementSolver {
  public:
    using SatSolver = BasicSatSolver;
    using Lit = typename SatSolver::Lit;
    using SLit = sammy::Lit;
    using LitOrVal = std::variant<Lit, bool>;

    static std::string name() {
        return std::string("SAT<") + BasicSatSolver::name() + ">";
    }

    /**
     * Estimate the total size of the formula (in bytes)
     * that would be created for a given subproblem.
     */
    static std::size_t estimate_formula_size(const ClauseDB* clause_db,
                                             const LNSSubproblem& subproblem) {
        // we need some estimate of the fixed memory per variable
        static constexpr std::size_t bytes_per_var = 16;
        // we need some estimate of the fixed memory per binary clause
        static constexpr std::size_t bytes_per_binary = 8;
        // we need some estimate of the fixed memory per long clause
        static constexpr std::size_t bytes_per_long = 4 + 2 * 8;
        // how many bytes does each longer-clause entry need?
        static constexpr std::size_t bytes_per_long_entry = 4;
        std::stringstream msg;

        std::size_t num_configurations =
            subproblem.removed_configurations.size() - 1;
        std::size_t num_vars = subproblem.removed_configurations[0].size();
        std::size_t num_uncovered = subproblem.uncovered_universe.size();
        std::size_t num_concrete = subproblem.num_concrete;
        std::size_t ncl = clause_db->num_clauses() - clause_db->num_unaries();
        std::size_t nbin = clause_db->num_binaries();
        std::size_t nlong = ncl - nbin;
        std::size_t nlong_entries = clause_db->literal_db_size() - nlong;

        // size for 'normal' variables (not coverage variables)
        std::size_t num_normal_vars = num_vars * num_configurations;
        std::size_t normal_var_size = num_normal_vars * bytes_per_var;
        // size for 'normal' clauses
        std::size_t num_normal_long = nlong * num_configurations;
        std::size_t num_normal_bin = nbin * num_configurations;
        std::size_t normal_clause_size =
            num_normal_bin * bytes_per_binary +
            num_normal_long * bytes_per_long +
            nlong_entries * num_configurations * bytes_per_long_entry;
        // size for 'vertical' clauses
        std::size_t vertical_clause_count = 2 * num_concrete;
        std::size_t vertical_clause_length = num_configurations;
        std::size_t vertical_clause_size =
            (bytes_per_long * vertical_clause_count) +
            (bytes_per_long_entry * vertical_clause_length *
             vertical_clause_count);
        // size for 'coverage' variables
        std::size_t coverage_var_size =
            num_uncovered * num_configurations * bytes_per_var;
        std::size_t coverage_var_clause_size =
            num_uncovered * num_configurations *
            (2 * bytes_per_binary + bytes_per_long + 3 * bytes_per_long_entry);
        std::size_t coverage_clause_size =
            num_uncovered *
            (bytes_per_long + bytes_per_long_entry * num_configurations);
        return normal_var_size + normal_clause_size + vertical_clause_size +
               coverage_var_size + coverage_var_clause_size +
               coverage_clause_size;
    }

    /**
     * Create a new solver for the given subproblem.
     */
    FixedMESSATImprovementSolver(PortfolioSolver* portfolio,
                                 LNSSubproblem&& subproblem,
                                 SharedDBPropagator prop,
                                 EventRecorder* recorder, std::size_t worker_id)
        : m_solver(), m_subproblem(std::move(subproblem)), m_recorder(recorder),
          m_worker_id(worker_id), m_propagator(std::move(prop)),
          m_num_configs_allowed(m_subproblem.removed_configurations.size() - 1),
          m_config_vars() {
        try {
            m_estimated_formula_size =
                estimate_formula_size(&m_propagator.db(), m_subproblem);
            OutputObject event_data{
                {"num_uncovered", m_subproblem.uncovered_universe.size()},
                {"num_removed", m_subproblem.removed_configurations.size()},
                {"mes_size", m_subproblem.mutually_exclusive_set.size()},
                {"solver_name", m_solver.name()},
                {"num_configs_allowed", m_num_configs_allowed}};
            p_store_event("SAT_LNS_CONSTRUCTION_BEGIN", event_data,
                          "num_uncovered", "num_removed", "mes_size",
                          "solver_name", "num_configs_allowed");
            if (m_num_configs_allowed <
                m_subproblem.mutually_exclusive_set.size())
            {
                m_infeasible_by_construction = true;
                return;
            }
            // make variables for configurations and ensure their consistency
            p_make_config_vars();
            assert(m_config_vars.size() == m_num_configs_allowed);
            assert(std::all_of(
                m_config_vars.begin(), m_config_vars.end(), [&](const auto& v) {
                    return v.size() == m_propagator.db().num_vars();
                }));
            // make clauses for all individual literals in uncovered
            // interactions
            p_make_vertical_clauses();
            // ensure that every uncovered universe element is covered
            p_make_coverage_clauses();
            p_store_event("SAT_LNS_CONSTRUCTION_DONE", std::move(event_data),
                          "num_uncovered", "num_removed", "mes_size",
                          "solver_name", "num_configs_allowed");
        } catch (InterruptError&) {
            subproblem = std::move(m_subproblem);
            throw;
        }
    }

    /**
     * Get the number of allowed configurations.
     */
    std::size_t allowed_configurations() const noexcept {
        return m_subproblem.removed_configurations.size() - 1;
    }

    /**
     * Get the number of vertices in the mutually exclusive set.
     */
    std::size_t mes_size() const noexcept {
        return m_subproblem.mutually_exclusive_set.size();
    }

    /**
     * Get the number of interactions in the uncovered universe.
     */
    std::size_t num_uncovered() const noexcept {
        return m_subproblem.uncovered_universe.size();
    }

    /**
     * Get the vertices in the mutually exclusive set.
     */
    const std::vector<Vertex>& mes_vertices() const noexcept {
        return m_subproblem.mutually_exclusive_set;
    }

    /**
     * Move out the subproblem, including the
     * (potentially huge) list of interactions.
     */
    LNSSubproblem move_out_subproblem() noexcept {
        return std::move(m_subproblem);
    }

    /**
     * Get the solution of the solver, after solve()
     * returned true.
     */
    const std::vector<DynamicBitset>& get_solution() const {
        return m_current_solution;
    }

    /**
     * Abort the solve.
     */
    void abort() { m_solver.terminate(); }

    /**
     * Solve the problem; return
     * nullopt if interrupted, false if
     * an improvement is impossible,
     * and true if an improvement was found.
     */
    std::optional<bool> solve() {
        if (m_infeasible_by_construction) {
            p_store_event("SAT_LNS_INFEASIBLE_BY_CONSTRUCTION");
            return false;
        }
        if (m_subproblem.uncovered_universe.empty()) {
            p_store_event("SAT_LNS_EMPTY_UNIVERSE");
            return true;
        }
        OutputObject event_data{
            {"num_uncovered", num_uncovered()},
            {"num_removed", m_subproblem.removed_configurations.size()},
            {"mes_size", mes_size()},
            {"solver_name", m_solver.name()}};
        if (get_and_clear_interrupt_flag()) {
            return std::nullopt;
        }
        p_store_event("SAT_LNS_BEGIN_SAT_SOLVE", event_data, "num_uncovered",
                      "num_removed", "mes_size", "solver_name");
        std::optional<bool> res = m_solver.solve();
        const char* event_name = !res ? "SAT_LNS_ABORT_SAT_SOLVE"
                                      : (*res ? "SAT_LNS_IMPROVEMENT_FOUND"
                                              : "SAT_LNS_SOLUTION_WAS_OPTIMAL");
        p_store_event(event_name, event_data, "num_uncovered", "num_removed",
                      "mes_size", "solver_name");
        if (!res) {
            return std::nullopt;
        }
        if (!*res) {
            return false;
        }
        p_extract_solution();
        return true;
    }

  private:
    // basic SAT-solver backend
    SatSolver m_solver;
    // subproblem information
    LNSSubproblem m_subproblem;
    // event recorder
    EventRecorder* m_recorder;
    // worker id for events
    std::size_t m_worker_id;
    // propagator for the SAT model
    SharedDBPropagator m_propagator;
    // flag to set when, during construction, we find infeasibility
    bool m_infeasible_by_construction{false};
    // how many configurations are needed?
    std::size_t m_num_configs_allowed;
    // (outer) index k: configuration index [0, allowed_configurations());
    // inner index: feature index [0, propagator().db().num_vars())
    std::vector<std::vector<LitOrVal>> m_config_vars;
    // buffer for solver literals/clauses
    std::vector<Lit> m_buffer;
    // buffer for the current solution
    std::vector<DynamicBitset> m_current_solution;
    // stats
    std::size_t m_added_binaries{0};
    std::size_t m_added_long{0};
    std::size_t m_added_literals{0};
    std::size_t m_estimated_formula_size{0};

    // ----------------------------- IMPLEMENTATION --------------------------
    /**
     * Store an event in the recorder.
     */
    template <typename... Args>
    void p_store_event(const char* name, OutputObject data, Args&&... keys) {
        if (!m_recorder)
            return;
        data["worker_id"] = m_worker_id;
        m_recorder->store_event(name, std::move(data),
                                std::forward<Args>(keys)...);
    }

    /**
     * Store an event in the recorder.
     */
    void p_store_event(const char* name) {
        p_store_event(name, {{"worker_id", m_worker_id}});
    }

    void p_make_config_vars() {
        for (std::size_t config_index = 0; config_index < m_num_configs_allowed;
             ++config_index)
        {
            m_propagator.reset_or_throw();
            if (config_index < m_subproblem.mutually_exclusive_set.size()) {
                if (push_vertex(
                        m_propagator,
                        m_subproblem.mutually_exclusive_set[config_index]) < 0)
                {
                    throw std::logic_error("Invalid interaction in MES!");
                }
            }
            p_make_config_from_propagator();
        }
        m_propagator.reset_or_throw();
        assert(m_num_configs_allowed ==
               m_subproblem.removed_configurations.size() - 1);
        assert(m_num_configs_allowed == m_config_vars.size());
    }

    void p_make_config_from_propagator() {
        const auto nall = m_propagator.db().num_vars();
        std::vector<LitOrVal> config;
        config.reserve(nall);
        for (Var var = 0; var < nall; ++var) {
            SLit pos = lit::positive_lit(var);
            if (m_propagator.is_true(pos)) {
                config.emplace_back(std::in_place_type<bool>, true);
            } else if (m_propagator.is_false(pos)) {
                config.emplace_back(std::in_place_type<bool>, false);
            } else {
                config.emplace_back(std::in_place_type<Lit>,
                                    m_solver.new_var());
            }
        }
        assert(config.size() == nall);
        m_config_vars.emplace_back(std::move(config));
        assert(m_config_vars.back().size() == nall);
        p_config_insert_binaries();
        p_config_insert_long_clauses();
    }

    void p_config_insert_binaries() {
        const auto& db = m_propagator.db();
        const auto& config = m_config_vars.back();
        assert(db.num_vars() == config.size());
        for (auto clause : db.binary_clauses()) {
            LitOrVal lv1 = p_config_lit_for(config, clause.first);
            LitOrVal lv2 = p_config_lit_for(config, clause.second);
            std::visit(
                overloaded{[&](bool b1, bool b2) {
                               // clause is either violated (indicating a logic
                               // error) or already satisfied by propagation
                               if (!b1 && !b2) {
                                   throw std::logic_error(
                                       "Infeasible binary clause!");
                               }
                           },
                           [&](bool b1, Lit l2) {
                               // if b1 is true, clause is satisfied;
                               // if b1 is false, propagation should have made
                               // l2 true.
                               if (!b1) {
                                   throw std::logic_error(
                                       "Propagation should have made l2 true!");
                               }
                           },
                           [&](Lit l1, bool b2) {
                               // if b2 is true, clause is satisfied;
                               // if b2 is false, propagation should have made
                               // l1 true.
                               if (!b2) {
                                   throw std::logic_error(
                                       "Propagation should have made l1 true!");
                               }
                           },
                           [&](Lit l1, Lit l2) {
                               m_solver.add_short_clause(l1, l2);
                               ++m_added_binaries;
                           }},
                lv1, lv2);
        }
    }

    void p_config_insert_long_clauses() {
        const auto& db = m_propagator.db();
        const auto& config = m_config_vars.back();
        for (CRef clause = 1, dbs = db.literal_db_size(); clause < dbs;
             clause = db.next_clause(clause))
        {
            m_buffer.clear();
            bool got_true = false;
            for (SLit l : db.lits_of(clause)) {
                LitOrVal v = p_config_lit_for(config, l);
                if (std::holds_alternative<bool>(v)) {
                    if (*std::get_if<bool>(&v)) {
                        got_true = true;
                        break;
                    }
                    continue;
                }
                m_buffer.push_back(*std::get_if<Lit>(&v));
            }
            if (got_true)
                continue;
            if (m_buffer.empty()) {
                throw std::logic_error("Infeasible interaction in MES!");
            }
            if (m_buffer.size() == 1) {
                throw std::logic_error("Unit clause in MES - propagation "
                                       "should have assigned a value!");
            }
            m_solver.add_clause(m_buffer.begin(), m_buffer.end());
            m_added_long += 1;
            m_added_literals += m_buffer.size();
        }
    }

    void p_make_vertical_clauses() {
        if (m_infeasible_by_construction)
            return;
        HashSet<SLit> occurring;
        for (Vertex v : m_subproblem.uncovered_universe) {
            occurring.insert(v.first);
            occurring.insert(v.second);
        }
        for (SLit l : occurring) {
            m_buffer.clear();
            bool got_true = false;
            for (const auto& config : m_config_vars) {
                assert(config.size() == m_propagator.db().num_vars());
                LitOrVal lv = p_config_lit_for(config, l);
                if (std::holds_alternative<bool>(lv)) {
                    if (*std::get_if<bool>(&lv)) {
                        got_true = true;
                        break;
                    }
                    continue;
                }
                m_buffer.push_back(*std::get_if<Lit>(&lv));
            }
            if (got_true)
                continue;
            if (m_buffer.empty()) {
                m_infeasible_by_construction = true;
                break;
            }
            m_solver.add_clause(m_buffer.begin(), m_buffer.end());
            m_added_long += 1;
            m_added_literals += m_buffer.size();
        }
    }

    LitOrVal p_config_lit_for(const std::vector<LitOrVal>& config,
                              SLit literal) const {
        auto var = lit::var(literal);
        bool is_pos = !lit::negative(literal);
        assert(var < config.size());
        auto& entry = config[var];
        return std::visit(
            overloaded{
                [&](bool b) -> LitOrVal {
                    return LitOrVal{std::in_place_type<bool>, is_pos ? b : !b};
                },
                [&](Lit l) -> LitOrVal {
                    return LitOrVal{std::in_place_type<Lit>, is_pos ? l : -l};
                }},
            entry);
    }

    /**
     * Extract the solution from the solver
     * and store it in m_current_solution.
     */
    void p_extract_solution() {
        const auto nclasses = m_config_vars.size();
        const auto nvars = m_propagator.db().num_vars();
        auto model_map = m_solver.get_model();
        std::vector<DynamicBitset> result(nclasses);
        for (std::size_t class_index = 0; class_index < nclasses; ++class_index)
        {
            DynamicBitset config(nvars, false);
            for (Var v = 0; v < nvars; ++v) {
                config[v] = p_get_class_value(model_map, class_index, v);
            }
            result[class_index] = std::move(config);
        }
        m_current_solution = std::move(result);
    }

    /**
     * Get the value of a literal or boolean
     * constant, using a model map.
     */
    template <typename ModelMap>
    bool p_get_value(const ModelMap& model, LitOrVal v) const noexcept {
        return std::visit(overloaded{[&](bool b) { return b; },
                                     [&](Lit l) { return model[l]; }},
                          v);
    }

    /**
     * Get the value of a variable in a class.
     */
    template <typename ModelMap>
    bool p_get_class_value(const ModelMap& model, std::size_t class_index,
                           Var var) const noexcept {
        assert(class_index < m_config_vars.size());
        assert(var < m_config_vars[class_index].size());
        return p_get_value(model, m_config_vars[class_index][var]);
    }

    /**
     * Add variables and clauses to ensure that
     * every uncovered universe element is covered.
     */
    void p_make_coverage_clauses() {
        std::size_t count = 0;
        for (Vertex v : m_subproblem.uncovered_universe) {
            if (++count == 32768) {
                count = 0;
                throw_if_interrupted();
            }
            bool got_true = false;
            m_buffer.clear();
            // first scan partially-assigned configurations;
            // this avoids variable creation for vertices that
            // are already covered in a configuration
            for (const auto& config : m_config_vars) {
                LitOrVal l1 = p_config_lit_for(config, v.first);
                LitOrVal l2 = p_config_lit_for(config, v.second);
                if (std::holds_alternative<bool>(l1) ||
                    std::holds_alternative<bool>(l2))
                {
                    bool l1true = false;
                    bool l2true = false;
                    if (std::holds_alternative<bool>(l1)) {
                        if (*std::get_if<bool>(&l1)) {
                            l1true = true;
                        } else {
                            // fixed false in this configuration
                            continue;
                        }
                    }
                    if (std::holds_alternative<bool>(l2)) {
                        if (*std::get_if<bool>(&l2)) {
                            l2true = true;
                        } else {
                            // fixed false in this configuration
                            continue;
                        }
                    }
                    if (l1true && l2true) {
                        got_true = true;
                        break;
                    }
                    if (l1true) {
                        m_buffer.push_back(*std::get_if<Lit>(&l2));
                    } else {
                        m_buffer.push_back(*std::get_if<Lit>(&l1));
                    }
                }
            }
            if (got_true)
                continue;
            // m_buffer contains the beginnings of a clause ensuring
            // that some configuration covers v; we now complete it
            // using extra variables for the other classes where
            // v is covered but needs two literals to be set properly
            for (const auto& config : m_config_vars) {
                LitOrVal l1v = p_config_lit_for(config, v.first);
                LitOrVal l2v = p_config_lit_for(config, v.second);
                if (std::holds_alternative<bool>(l1v) ||
                    std::holds_alternative<bool>(l2v))
                {
                    continue;
                }
                Lit l1 = *std::get_if<Lit>(&l1v);
                Lit l2 = *std::get_if<Lit>(&l2v);
                Lit coverage_var = m_solver.new_var();
                m_solver.add_short_clause(coverage_var, -l1, -l2);
                m_solver.add_short_clause(-coverage_var, l1);
                m_solver.add_short_clause(-coverage_var, l2);
                m_added_long += 1;
                m_added_literals += 3;
                m_added_binaries += 2;
                m_buffer.push_back(coverage_var);
            }
            if (m_buffer.empty()) {
                m_infeasible_by_construction = true;
                break;
            }
            m_solver.add_clause(m_buffer.begin(), m_buffer.end());
            m_added_long += 1;
            m_added_literals += m_buffer.size();
        }
    }
};

/**
 * Non-incremental use of a SAT solver for improvement-by-one search.
 */
template <typename BaseSolver> class ImprovementSolver {
  public:
    using Lit = typename BaseSolver::Lit;
    using SLit = sammy::Lit;

    ImprovementSolver(SharedDBPropagator* propagator,
                      std::vector<Vertex> vertices,
                      std::vector<std::size_t> clique_indices, std::size_t k)
        : m_base_solver(), m_propagator(propagator),
          m_all_vertices(std::move(vertices)),
          m_clique_indices(std::move(clique_indices)), m_class_vars(k),
          m_vertex_vars(m_all_vertices.size()) {
        for (std::size_t i = 0; i < k; ++i) {
            if (i <= m_clique_indices.size())
                m_propagator->reset_or_throw();
            if (i < m_clique_indices.size()) {
                p_make_class_with_clique_vertex(
                    i, m_all_vertices[m_clique_indices[i]]);
            } else {
                p_make_unconstrained_class(i);
            }
            if (!p_clauses_for_class(i)) {
                m_construction_infeasible = true;
                break;
            }
        }
        if (!m_construction_infeasible) {
            p_clauses_for_lits();
        }
        if (!m_construction_infeasible) {
            for (std::size_t i = 0, n = m_all_vertices.size(); i < n; ++i) {
                if (!p_vertex_clause(i)) {
                    m_construction_infeasible = true;
                    break;
                }
            }
        }
    }

    std::optional<bool>
    solve(double time_limit = std::numeric_limits<double>::infinity()) {
        return m_base_solver.solve(time_limit);
    }

    std::vector<std::vector<bool>> get_solution() const {
        const auto nclasses = m_class_vars.size();
        const auto nvars = m_propagator->db().num_vars();
        auto model_map = m_base_solver.get_model();
        std::vector<std::vector<bool>> result(nclasses);
        for (std::size_t class_index = 0; class_index < nclasses; ++class_index)
        {
            std::vector<bool> config(nvars, false);
            for (Var v = 0; v < nvars; ++v) {
                config[v] = p_get_class_value(model_map, class_index, v);
            }
            result[class_index] = std::move(config);
        }
        return result;
    }

  private:
    template <typename ModelMap>
    bool p_get_value(const ModelMap& model, std::variant<Lit, bool> v) const {
        return std::visit(overloaded{[&](bool& b) { return b; },
                                     [&](Lit& l) { return model[l]; }},
                          v);
    }

    template <typename ModelMap>
    bool p_get_class_value(const ModelMap& model, std::size_t class_index,
                           Var var) const {
        return p_get_value(model, m_class_vars[class_index][var]);
    }

    bool p_vertex_clause(std::size_t vertex_index) {
        auto& vvars = m_vertex_vars[vertex_index];
        m_clause_buffer.clear();
        for (auto& var : vvars) {
            if (!std::visit(overloaded{[&](bool& b) -> bool { return !b; },
                                       [&](Lit& l) -> bool {
                                           m_clause_buffer.push_back(l);
                                           return true;
                                       }},
                            var))
            {
                return true;
            }
        }
        if (m_clause_buffer.empty())
            return false;
        m_base_solver.add_clause(m_clause_buffer.begin(),
                                 m_clause_buffer.end());
        return true;
    }

    std::variant<Lit, bool> p_class_var_for(std::size_t class_index,
                                            SLit literal) {
        Var v = lit::var(literal);
        bool neg = lit::negative(literal);
        std::variant<Lit, bool> cv = m_class_vars[class_index][v];
        std::visit(overloaded{[&](bool& b) {
                                  if (neg)
                                      b = !b;
                              },
                              [&](Lit& l) {
                                  if (neg)
                                      l = -l;
                              }},
                   cv);
        return cv;
    }

    void p_class_vars_from_prop(std::size_t i) {
        const auto& clauses = m_propagator->db();
        auto& cvars = m_class_vars[i];
        for (Var v = 0, nv = clauses.num_vars(); v < nv; ++v) {
            if (m_propagator->is_open(lit::positive_lit(v))) {
                cvars.emplace_back(std::in_place_index<0>,
                                   m_base_solver.new_var());
            } else {
                cvars.emplace_back(std::in_place_index<1>,
                                   m_propagator->is_true(lit::positive_lit(v)));
            }
        }
    }

    void p_vertex_vars_from_prop(std::size_t class_index) {
        for (std::size_t i = 0, n = m_all_vertices.size(); i < n; ++i) {
            Vertex v = m_all_vertices[i];
            auto v1 = p_class_var_for(class_index, v.first);
            auto v2 = p_class_var_for(class_index, v.second);
            auto& vvars = m_vertex_vars[i];
            auto one_fixed = [&](Lit& l1, bool& b2) {
                if (!b2) {
                    vvars.emplace_back(std::in_place_index<1>, false);
                } else {
                    vvars.emplace_back(std::in_place_index<0>, l1);
                }
            };
            std::visit(
                overloaded{[&](bool& b1, bool& b2) {
                               vvars.emplace_back(std::in_place_index<1>,
                                                  b1 && b2);
                           },
                           [&](Lit& l1, Lit& l2) {
                               Lit vvar = m_base_solver.new_var();
                               vvars.emplace_back(std::in_place_index<0>, vvar);
                               m_base_solver.add_short_clause(-vvar, l1);
                               m_base_solver.add_short_clause(-vvar, l2);
                               m_base_solver.add_short_clause(-l1, -l2, vvar);
                           },
                           [&](Lit& l1, bool& b2) { one_fixed(l1, b2); },
                           [&](bool& b1, Lit& l2) { one_fixed(l2, b1); }},
                v1, v2);
        }
    }

    void p_make_unconstrained_class(std::size_t i) {
        p_class_vars_from_prop(i);
        p_vertex_vars_from_prop(i);
    }

    void p_make_class_with_clique_vertex(std::size_t class_index,
                                         Vertex clique_vertex) {
        if (m_propagator->is_open(clique_vertex.first))
            m_propagator->push_level(clique_vertex.first);
        if (m_propagator->is_open(clique_vertex.second))
            m_propagator->push_level(clique_vertex.second);
        p_class_vars_from_prop(class_index);
        p_vertex_vars_from_prop(class_index);
    }

    bool p_clauses_for_class(std::size_t i) {
        const auto& db = m_propagator->db();
        for (auto [l1, l2] : db.binary_clauses()) {
            SLit ls[2] = {l1, l2};
            if (!p_to_clause(i, +ls, ls + 2))
                return false;
        }
        for (CRef c = 1, n = db.literal_db_size(); c < n; c = db.next_clause(c))
        {
            auto lits = db.lits_of(c);
            if (!p_to_clause(i, lits.begin(), lits.end()))
                return false;
        }
        return true;
    }

    template <typename Iterator>
    bool p_to_clause(std::size_t class_index, Iterator slit_begin,
                     Iterator slit_end) {
        bool found_true = false;
        m_lit_buffer.clear();
        for (SLit sl : IteratorRange{slit_begin, slit_end}) {
            if (m_propagator->is_open(sl)) {
                m_lit_buffer.push_back(sl);
            } else if (m_propagator->is_true(sl)) {
                found_true = true;
                break;
            }
        }
        if (found_true)
            return true;
        if (m_lit_buffer.empty())
            return false;
        m_clause_buffer.clear();
        std::transform(m_lit_buffer.begin(), m_lit_buffer.end(),
                       std::back_inserter(m_clause_buffer), [&](SLit sl) {
                           return std::get<Lit>(
                               p_class_var_for(class_index, sl));
                       });
        m_base_solver.add_clause(m_clause_buffer.begin(),
                                 m_clause_buffer.end());
        return true;
    }

    /**
     * For each literal that occurs in any
     * vertex that we have to cover: add a clause
     * that ensures that at least one configuration
     * has that literal set to true, unless that
     * is already trivially satisfied.
     */
    void p_clauses_for_lits() {
        {
            HashSet<SLit> occurring;
            for (Vertex v : m_all_vertices) {
                occurring.insert(v.first);
                occurring.insert(v.second);
            }
            m_lit_buffer.reserve(occurring.size());
            m_lit_buffer.assign(occurring.begin(), occurring.end());
        }
        std::sort(m_lit_buffer.begin(), m_lit_buffer.end());
        for (SLit l : m_lit_buffer) {
            m_clause_buffer.clear();
            bool got_true = false;
            for (std::size_t config_index = 0,
                             num_configs = m_class_vars.size();
                 config_index < num_configs; ++config_index)
            {
                auto lv = p_class_var_for(config_index, l);
                std::visit(
                    overloaded{[&](bool b) {
                                   if (b) {
                                       got_true = true;
                                   }
                               },
                               [&](Lit l) { m_clause_buffer.push_back(l); }},
                    lv);
                if (got_true)
                    break;
            }
            if (got_true)
                continue;
            if (m_clause_buffer.empty()) {
                m_construction_infeasible = true;
                break;
            }
            m_base_solver.add_clause(m_clause_buffer.begin(),
                                     m_clause_buffer.end());
        }
        m_lit_buffer.clear();
    }

    BaseSolver m_base_solver;
    SharedDBPropagator* m_propagator;
    std::vector<Vertex> m_all_vertices;
    std::vector<std::size_t> m_clique_indices;
    // m_class_vars[c][v] -> variable (or fixed value) for class c and basic
    // variable v
    std::vector<std::vector<std::variant<Lit, bool>>> m_class_vars;
    // m_vertex_vars[vi][c] -> variable (or fixed value) for class c and vertex
    // with index vi
    std::vector<std::vector<std::variant<Lit, bool>>> m_vertex_vars;
    bool m_construction_infeasible = false;
    std::vector<SLit> m_lit_buffer;
    std::vector<Lit> m_clause_buffer;
};

} // namespace sammy

#endif
==> ./simplification_stats.h <==
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
==> ./vivify.h <==
#ifndef SAMMY_VIVIFY_H_INCLUDED_
#define SAMMY_VIVIFY_H_INCLUDED_

#include "shared_db_propagator.h"
#include "simplification_stats.h"
#include "simplify_datastructure.h"
#include <cassert>

namespace sammy {

class ClauseVivifier {
  public:
    explicit ClauseVivifier(SimplifyDatastructure* simplify_ds)
        : simplify_ds(simplify_ds), clauses((simplify_ds->sort_clauses(),
                                             simplify_ds->extract_clause_db())),
          propagator(&clauses), first_long_index(p_find_first_long()),
          stats(&stats_buffer) {}

    void set_stats(SimplificationStats* stats) { this->stats = stats; }

    bool run() {
        if (!propagator.get_trail().empty())
            return true;
        CRef c = 1, n = clauses.literal_db_size();
        CRef coffset = 0;
        bool changed = false;
        abort_run = false;
        for (; c < n; ++coffset, c = clauses.next_clause(c)) {
            auto size = clauses.clause_length(c);
            if (size > 3) {
                break;
            }
            changed |= p_vivify_ternary(clauses.lits_of(c), coffset);
            if (abort_run)
                return changed;
        }
        for (; c < n; ++coffset, c = clauses.next_clause(c)) {
            changed |= p_vivify_longer(clauses.lits_of(c), coffset);
            if (abort_run)
                return changed;
        }
        return changed;
    }

  private:
    CRef p_find_first_long() {
        auto& sc = simplify_ds->clauses();
        SCVec dummy(3, NIL);
        auto first_long_iter = std::lower_bound(
            sc.begin(), sc.end(), dummy, [](const SCVec& v1, const SCVec& v2) {
                return v1.size() < v2.size();
            });
        return CRef(first_long_iter - sc.begin());
    }

    bool p_vivify_longer(ClauseDB::Lits lits_, CRef coffset) {
        const Lit* lits = lits_.begin();
        auto [first_diff_pushed, first_diff_lits] =
            std::mismatch(pushed.begin(), pushed.end(), lits);
        if (first_diff_pushed == pushed.end()) {
            // the latter clause is subsumed by the first;
            // let subsumption detection handle this
            return false;
        }
        p_vivify_longer_pop_mismatched(first_diff_pushed);
        if (abort_run) {
            return true;
        }
        auto& scclause = simplify_ds->clauses()[first_long_index + coffset];
        for (auto pos = first_diff_lits; pos != lits_.end() - 1; ++pos) {
            Lit lcur = *pos;
            if (!propagator.is_open(lcur)) {
                if (propagator.is_true(lcur)) {
                    // (-l1 & ... & -lk -> lcur) => shorten clause
                    scclause.assign(lits_.begin(), pos + 1);
                } else {
                    // (-l1 & ... & -lk -> -lcur) === (l1 v ... v lk v -lcur)
                    // => remove lcur by resolution
                    scclause.assign(lits_.begin(), pos);
                    scclause.insert(scclause.end(), pos + 1, lits_.end());
                    clauses.add_clause(scclause.begin(), scclause.end());
                }
                stats->clauses_strengthened_by_vivification += 1;
                return true;
            }
            if (!propagator.push_level(lit::negate(lcur))) {
                // (-l1 & ... & -lk & -lcur) -> conflict => shorten clause
                scclause.assign(lits_.begin(), pos + 1);
                clauses.add_clause(scclause.begin(), scclause.end());
                stats->clauses_strengthened_by_vivification += 1;
                propagator.pop_level();
                return true;
            }
            pushed.push_back(lcur);
        }
        return false;
    }

    void p_vivify_longer_pop_mismatched(
        std::vector<Lit>::const_iterator first_diff_pushed) {
        // pop the additional levels from the propagator
        for (auto i = first_diff_pushed; i != pushed.end(); ++i) {
            propagator.pop_level();
        }
        pushed.erase(first_diff_pushed, pushed.end());
        if (pushed.empty()) {
            propagator.incorporate_or_throw();
            if (!propagator.get_trail().empty()) {
                abort_run = true;
            }
        }
    }

    bool p_vivify_ternary(ClauseDB::Lits lits_, CRef coffset) {
        assert(propagator.get_current_level() == 0);
        const Lit* lits = lits_.begin();
        return p_vivify_ternary_pair_conflicts(lits[0], lits[1], lits[2],
                                               coffset) ||
               p_vivify_ternary_pair_conflicts(lits[0], lits[2], lits[1],
                                               coffset) ||
               p_vivify_ternary_pair_conflicts(lits[1], lits[2], lits[0],
                                               coffset);
    }

    bool p_vivify_ternary_pair_conflicts(Lit l1, Lit l2, Lit l3, CRef coffset) {
        if (!propagator.push_level(lit::negate(l1))) {
            stats->failed_literals += 1;
            propagator.pop_level();
            simplify_ds->clauses().emplace_back(1, l1);
            abort_run = true;
            return true;
        }
        Lit replace_clause[2] = {l1, l2};
        bool replace = false;
        if (!propagator.is_open(l2)) {
            // -l1 -> l2 or -l1 -> -l2
            // -l1 -> l2 === l1 v l2 => subsumes clause.
            // -l1 -> -l2 === l1 v -l2 => resolution with clause on l2 gives l1
            // v l3
            replace = true;
            if (propagator.is_false(l2)) {
                replace_clause[1] = l3;
            }
        }
        if (!replace) {
            if (!propagator.push_level(lit::negate(l2))) {
                // -l1 & -l2 -> conflict
                // learn l1 v l2 => subsumes clause
                replace = true;
            }
            propagator.pop_level();
        }
        if (replace) {
            clauses.add_clause(+replace_clause, replace_clause + 2);
            simplify_ds->clauses()[first_long_index + coffset].assign(
                +replace_clause, replace_clause + 2);
            stats->clauses_strengthened_by_vivification += 1;
        }
        propagator.pop_level();
        if (replace) {
            propagator.incorporate_or_throw();
            if (!propagator.get_trail().empty()) {
                abort_run = true;
            }
        }
        return replace;
    }

    SimplifyDatastructure* simplify_ds;
    ClauseDB clauses;
    SharedDBPropagator propagator;
    std::vector<Lit> pushed;
    CRef first_long_index;
    SimplificationStats stats_buffer;
    SimplificationStats* stats;
    bool abort_run = false;
};

inline bool vivify(SimplifyDatastructure& simplifier,
                   SimplificationStats* stats = nullptr) {
    ClauseVivifier vivifier(&simplifier);
    if (stats)
        vivifier.set_stats(stats);
    return vivifier.run();
}

} // namespace sammy

#endif
==> ./initial_coloring_heuristic.h <==
#ifndef SAMMY_INITIAL_COLORING_HEURISTIC_H_INCLUDED_
#define SAMMY_INITIAL_COLORING_HEURISTIC_H_INCLUDED_

#include "algorithm_ex.h"
#include "class_completer.h"
#include "cuda_iteration.h"
#include "error.h"
#include "greedysat.h"
#include "literals.h"
#include "pair_infeasibility_map.h"
#include "parallel_bit_filter.h"
#include "partial_solution.h"
#include <boost/iterator/transform_iterator.hpp>
#include <unordered_set>

namespace sammy {

#ifdef SAMMY_CUDA_SUPPORTED
namespace detail {

class CUDABitFilter {
  public:
    CUDABitFilter(const std::vector<DynamicBitset>& literals_in_class,
                  const std::vector<std::vector<Index>>& classes_with_literal,
                  const PairInfeasibilityMap* inf_map)
        : m_inf_map(inf_map), m_u32_per_bitset((*inf_map)[0].blocks().size() *
                                               sizeof(DynamicBitset::Block) /
                                               sizeof(std::uint32_t)),
          m_device_bit_data((p_fill_prepare_buffer(literals_in_class),
                             m_host_prepare_buffer)),
          m_device_classes_with_literal(
              (p_fill_offset_buffer(classes_with_literal),
               m_host_prepare_buffer)),
          m_device_classes_with_literal_offsets(m_host_prepare_offsets),
          m_device_output_buffer(GOAL_ROWS_PER_CALL() * m_u32_per_bitset) {}

    template <typename Callable /*(Lit lmin, Lit lmax)*/>
    void iterate_all_uncovered(Callable&& callable) {
        Lit n_concrete = m_inf_map->get_n_concrete();
        Lit nclit = 2 * n_concrete;
        Lit lmin = 0;
        std::vector<std::uint32_t> host_output_buffer;
        DynamicBitset row_buffer(nclit, true);
        while (lmin < nclit - 2) {
            Lit num_rows = nclit - lmin - 2;
            if (num_rows > GOAL_ROWS_PER_CALL()) {
                num_rows = GOAL_ROWS_PER_CALL();
            }
            cuda_call_bit_filter_kernel(
                m_device_bit_data.get(), m_u32_per_bitset,
                m_device_classes_with_literal.get(), nclit,
                m_device_classes_with_literal_offsets.get(),
                m_device_output_buffer.get(), lmin, num_rows);
            m_device_output_buffer.copy_to_host(host_output_buffer);
            for (Lit l = 0; l < num_rows; ++l) {
                Lit lm = lmin + l;
                row_buffer.set();
                row_buffer ^= (*m_inf_map)[lm];
                row_buffer.binary_subtract(
                    &host_output_buffer[l * m_u32_per_bitset]);
                for (Lit lmax : row_buffer.ones_from(lm + 1)) {
                    callable(lm, lmax);
                }
            }
            lmin += num_rows;
        }
    }

    std::size_t count_uncovered() {
        Lit n_concrete = m_inf_map->get_n_concrete();
        Lit nclit = 2 * n_concrete;
        Lit lmin = 0;
        std::vector<std::uint32_t> host_output_buffer;
        DynamicBitset row_buffer(nclit, true);
        std::size_t result = 0;
        while (lmin < nclit - 2) {
            Lit num_rows = nclit - lmin - 2;
            if (num_rows > GOAL_ROWS_PER_CALL()) {
                num_rows = GOAL_ROWS_PER_CALL();
            }
            cuda_call_bit_filter_kernel(
                m_device_bit_data.get(), m_u32_per_bitset,
                m_device_classes_with_literal.get(), nclit,
                m_device_classes_with_literal_offsets.get(),
                m_device_output_buffer.get(), lmin, num_rows);
            m_device_output_buffer.copy_to_host(host_output_buffer);
            for (Lit l = 0; l < num_rows; ++l) {
                Lit lm = lmin + l;
                row_buffer.set();
                row_buffer ^= (*m_inf_map)[lm];
                row_buffer.binary_subtract(
                    &host_output_buffer[l * m_u32_per_bitset]);
                result += row_buffer.count_from(lm + 1);
            }
            lmin += num_rows;
        }
        return result;
    }

  private:
    void
    p_fill_prepare_buffer(const std::vector<DynamicBitset>& literals_in_class) {
        if (literals_in_class.empty()) {
            throw std::invalid_argument("literals_in_class must not be empty");
        }
        m_host_prepare_buffer.reserve(literals_in_class.size() *
                                      m_u32_per_bitset);
        for (const DynamicBitset& bs : literals_in_class) {
            for (DynamicBitset::Block b : bs.blocks()) {
                to_prepare_buffer(b, m_host_prepare_buffer);
            }
        }
    }

    void p_fill_offset_buffer(
        const std::vector<std::vector<Index>>& classes_with_literal) {
        m_host_prepare_buffer.clear();
        m_host_prepare_offsets.reserve(classes_with_literal.size() + 1);
        m_host_prepare_offsets.push_back(0);
        for (const std::vector<Index>& cls : classes_with_literal) {
            std::copy(cls.begin(), cls.end(),
                      std::back_inserter(m_host_prepare_buffer));
            m_host_prepare_offsets.push_back(m_host_prepare_buffer.size());
        }
    }

    // pair infeasibility map
    const PairInfeasibilityMap* m_inf_map;
    // words per bitset in m_host_prepare_buffer and on the device
    std::size_t m_u32_per_bitset;
    // buffer on the host to prepare the bitsets for copy,
    // or the classes with literal list
    std::vector<std::uint32_t> m_host_prepare_buffer;
    // host buffer for the offsets in classes with literal
    std::vector<std::uint32_t> m_host_prepare_offsets;

    // pointer to array of concatenated bitsets
    CUDADevicePointer<const std::uint32_t> m_device_bit_data;
    // pointer to array of classes with literal
    CUDADevicePointer<const std::uint32_t> m_device_classes_with_literal;
    // pointer to offsets within m_device_classes_with_literal to the start of
    // each literal
    CUDADevicePointer<const std::uint32_t>
        m_device_classes_with_literal_offsets;
    // pointer to output matrix rows
    CUDADevicePointer<std::uint32_t> m_device_output_buffer;
};

} // namespace detail

template <typename Callable>
static inline void cuda_iterate_all_uncovered(
    const std::vector<DynamicBitset>& literals_in_class,
    const std::vector<std::vector<Index>>& classes_with_literal,
    const PairInfeasibilityMap* inf_map, Callable&& callable) {
    if (literals_in_class.empty()) {
        Lit nclit = 2 * inf_map->get_n_concrete();
        DynamicBitset row_buffer(nclit, true);
        for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
            row_buffer.set();
            row_buffer ^= (*inf_map)[lmin];
            for (Lit lmax : row_buffer.ones_from(lmin + 1)) {
                std::forward<Callable>(callable)(lmin, lmax);
            }
        }
        return;
    }

    detail::CUDABitFilter bit_filter(literals_in_class, classes_with_literal,
                                     inf_map);
    bit_filter.iterate_all_uncovered(std::forward<Callable>(callable));
}

static inline std::size_t cuda_count_uncovered(
    const std::vector<DynamicBitset>& literals_in_class,
    const std::vector<std::vector<Index>>& classes_with_literal,
    const PairInfeasibilityMap* inf_map) {
    if (literals_in_class.empty()) {
        Lit nclit = 2 * inf_map->get_n_concrete();
        DynamicBitset row_buffer(nclit, true);
        std::size_t result = 0;
        for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
            row_buffer.set();
            row_buffer ^= (*inf_map)[lmin];
            result += row_buffer.count_from(lmin + 1);
        }
        return result;
    }

    detail::CUDABitFilter bit_filter(literals_in_class, classes_with_literal,
                                     inf_map);
    return bit_filter.count_uncovered();
}
#endif

/**
 * A class containing a set of propagators,
 * managing the incremental construction of
 * several configurations simultaneously.
 */
template <typename EventListener> class ColorClasses {
    /**
     * @brief A handler for the various events that
     *        can occur during our operations.
     */
    EventListener* event_handler;

    /**
     * @brief The database of all clauses,
     *        shared by all propagators.
     */
    ClauseDB* all_clauses;

    /**
     * @brief A propagator that is kept 'empty',
     *        i.e., with no (permanent) decisions;
     *        it gets all clauses and serves as a source to
     *        copy-init new color classes from
     *        (much faster than re-initialization).
     */
    SharedDBPropagator empty_class;

    /**
     * @brief Each color class is represented by the
     *        literals in the trail of a propagator.
     */
    std::vector<SharedDBPropagator> color_classes;

    /**
     * @brief For each literal, a list of class indices
     *        in which that literal is forced to true.
     */
    std::vector<std::vector<Index>> classes_with_literal;

    /**
     * @brief The list of vertices that spawned the color classes,
     *        i.e., that were first given that color.
     */
    std::vector<Vertex> class_spawners;

    /**
     * @brief An internal buffer of the phases
     *        of variables for p_try_completion.
     */
    std::vector<bool> m_phases;

    /**
     * @brief An internal buffer of the decisions
     *        for p_try_completion.
     */
    std::vector<Lit> m_decisions;

    /**
     * @brief Add a new color class.
     *
     * @param prop The propagator for the new class.
     * @param spawner The vertex (feature literal pair) that spawned the class.
     */
    void p_add_class(SharedDBPropagator&& prop, Vertex spawner) {
        class_spawners.emplace_back((std::min)(spawner.first, spawner.second),
                                    (std::max)(spawner.first, spawner.second));
        color_classes.emplace_back(std::move(prop));
        Index idx = num_classes() - 1;
        for (Lit l : color_classes.back().get_trail()) {
            classes_with_literal[l].push_back(idx);
        }
        event_handler->new_color_class(idx);
    }

    bool p_spawn_class_one_is_false(Vertex spawner) const noexcept {
        if (empty_class.is_false(spawner.first) ||
            empty_class.is_false(spawner.second))
        {
            if (empty_class.is_false(spawner.first)) {
                event_handler->feature_literal_infeasible(spawner.first);
            }
            if (empty_class.is_false(spawner.second)) {
                event_handler->feature_literal_infeasible(spawner.second);
            }
            return true;
        }
        return false;
    }

    void p_init_phases(SharedDBPropagator& propagator) {
        const Lit nv = all_clauses->num_vars();
        m_phases.assign(nv, false);
        for (Lit v = 0; v < nv; ++v) {
            if (propagator.is_true(lit::positive_lit(v))) {
                m_phases[v] = true;
            }
        }
    }

    /**
     * @brief Get a spawner from a propagator by finding the first concrete
     * literals on different levels.
     */
    std::optional<Vertex> p_spawner_from_levels(const SharedDBPropagator& prop,
                                                Lit n_concrete) {
        const std::uint32_t nlvl = prop.get_current_level();
        if (nlvl < 2)
            return std::nullopt;
        const Lit nclit = 2 * n_concrete;
        auto&& is_concrete = [&](Lit l) { return l < nclit; };
        for (std::uint32_t level = 1; level <= nlvl; ++level) {
            auto lbeg = prop.level_begin(level);
            auto lend = prop.level_end(level);
            auto lf = std::find_if(lbeg, lend, is_concrete);
            if (lf != lend) {
                for (++level; level <= nlvl; ++level) {
                    auto lbeg2 = prop.level_begin(level);
                    auto lend2 = prop.level_end(level);
                    auto lf2 = std::find_if(lbeg2, lend2, is_concrete);
                    if (lf2 != lend2) {
                        Lit l1 = *lf, l2 = *lf2;
                        return std::optional<Vertex>{
                            std::in_place, std::min(l1, l2), std::max(l1, l2)};
                    }
                }
                break;
            }
        }
        return std::nullopt;
    }

    /**
     * @brief Called if conflict resolution took us below
     *        the level introduced by the last assumption (decision
     *        that we wanted to keep) in p_try_completion.
     *
     * @param propagator
     * @return true If we could fix the situation, i.e., reintroduce a decision,
     *              or a decision turned into an implied literal.
     * @return false If the set of assumptions is infeasible.
     */
    bool p_try_completion_jumped_decision(SharedDBPropagator& propagator) {
        for (Lit l : m_decisions) {
            // if any of our assumptions became forced to false, return false.
            if (propagator.is_false(l))
                return false;
            // if any of our assumptions became open, try closing it.
            if (propagator.is_open(l)) {
                // if that yields a conflict, resolve and return false.
                if (!propagator.push_level(l)) {
                    if (!propagator.resolve_conflicts())
                        throw UNSATError();
                    return false;
                }
            }
        }
        // remove all decisions that are decisions no more
        m_decisions.erase(
            std::remove_if(m_decisions.begin(), m_decisions.end(),
                           [&](Lit d) { return !propagator.is_decision(d); }),
            m_decisions.end());
        return true;
    }

    bool p_try_completion_push(SharedDBPropagator& propagator, Lit v) {
        Lit d = m_phases[v] ? lit::positive_lit(v) : lit::negative_lit(v);
        if (!propagator.push_level(d)) {
            m_phases[v].flip();
            if (!propagator.resolve_conflicts()) {
                throw UNSATError();
            }
            if (std::size_t(propagator.get_current_level()) <
                m_decisions.size())
            {
                if (!p_try_completion_jumped_decision(propagator))
                    return false;
            }
            if (propagator.is_open(d)) {
                return p_try_completion_push(propagator, v);
            }
        }
        return true;
    }

    /**
     * @brief Try to complete the given partial assignment
     *        (actually solving SAT with a very simple variable selection
     * strategy).
     *
     * @param propagator
     * @return true If the solve was successful (the propagator is returned to a
     * state with the same assumptions as it had initially).
     * @return false If the solve was unsuccessful (the propagator is in some
     * unspecified state); the given assumptions are infeasible, and clauses
     * have been learned to reflect that.
     */
    bool p_try_completion(SharedDBPropagator& propagator) {
        p_init_phases(propagator);
        const Lit nv = all_clauses->num_vars();
        m_decisions = propagator.get_decisions();
        bool any_open = true;
        while (any_open) {
            any_open = false;
            for (Lit v = 0; v < nv; ++v) {
                Lit p = lit::positive_lit(v);
                if (propagator.is_open(p)) {
                    any_open = true;
                    if (!p_try_completion_push(propagator, v))
                        return false;
                }
            }
        }
        while (std::size_t(propagator.get_current_level()) > m_decisions.size())
        {
            propagator.pop_level();
        }
        return true;
    }

    struct ConflictResolutionHandler {
        void assignment_forced(Lit l) const {
            that->classes_with_literal[l].push_back(cls_index);
            that->event_handler->literal_added_to_class(cls_index, l);
        }

        void assignment_undone(Lit l) const {
            auto& list = that->classes_with_literal[l];
            list.erase(std::remove(list.begin(), list.end(), cls_index),
                       list.end());
            that->event_handler->literal_removed_from_class(cls_index, l);
        }

        ColorClasses* that;
        Index cls_index;
    };

  public:
    std::vector<SharedDBPropagator>
    remove_classes(const std::vector<std::size_t>& sorted_classes) {
        if (sorted_classes.empty())
            return {};
        std::vector<SharedDBPropagator> result;
        result.reserve(sorted_classes.size());
        for (std::size_t idx : sorted_classes) {
            result.emplace_back(std::move(color_classes[idx]));
        }
        class_spawners.erase(
            remove_indices(class_spawners.begin(), class_spawners.end(),
                           sorted_classes.begin(), sorted_classes.end()),
            class_spawners.end());
        const auto n = color_classes.size();
        std::vector<Index> old_to_new(n, std::numeric_limits<Index>::max());
        auto iter = sorted_classes.begin();
        Index new_out = 0;
        for (Index i = 0; i < n; ++i) {
            if (iter != sorted_classes.end() && *iter == i) {
                ++iter;
            } else {
                old_to_new[i] = new_out++;
            }
        }
        color_classes.erase(
            remove_indices(color_classes.begin(), color_classes.end(),
                           sorted_classes.begin(), sorted_classes.end()),
            color_classes.end());
        for (auto& cwl : classes_with_literal) {
            std::transform(cwl.begin(), cwl.end(), cwl.begin(),
                           [&](Index c) { return old_to_new[c]; });
            cwl.erase(std::remove(cwl.begin(), cwl.end(),
                                  std::numeric_limits<Index>::max()),
                      cwl.end());
        }
        return result;
    }

    const std::vector<Vertex>& spawners() const noexcept {
        return class_spawners;
    }

    void class_changed(Index cls_index, const Bitset& old_true,
                       const Bitset& new_true, Bitset& buffer) {
        buffer = old_true;
        buffer -= new_true;
        for (Lit l : buffer.ones()) {
            auto& list = classes_with_literal[l];
            list.erase(std::remove(list.begin(), list.end(), cls_index),
                       list.end());
            event_handler->literal_removed_from_class(cls_index, l);
        }
        buffer = new_true;
        buffer -= old_true;
        for (Lit l : buffer.ones()) {
            classes_with_literal[l].push_back(cls_index);
            event_handler->literal_added_to_class(cls_index, l);
        }
    }

    const std::vector<std::vector<Index>>&
    get_classes_with_literal() const noexcept {
        return classes_with_literal;
    }

    const std::vector<Index>& with_literal(Lit l) const noexcept {
        return classes_with_literal[l];
    }

    IteratorRange<std::vector<Lit>::const_iterator>
    pretend_push(Lit l, Index cls_index) {
        auto& cls = color_classes[cls_index];
        if (!cls.push_level(l)) {
            return {cls.get_trail().end(), cls.get_trail().end()};
        } else {
            return {cls.current_level_begin(), cls.get_trail().end()};
        }
    }

    void pretend_pop(Index cls_index) {
        auto& cls = color_classes[cls_index];
        cls.pop_level();
    }

    explicit ColorClasses(EventListener* event_handler, ClauseDB* all_clauses)
        : event_handler(event_handler), all_clauses(all_clauses),
          empty_class(all_clauses),
          classes_with_literal(2 * all_clauses->num_vars()) {}

    bool push_vertex(Index cindex, Lit lmin, Lit lmax) {
        SharedDBPropagator& prop = color_classes[cindex];
        bool min_true = prop.is_true(lmin);
        bool max_true = prop.is_true(lmax);
        bool pushed_min = false, pushed_max = false;
        if (min_true & max_true) {
            return true;
        }
        if (prop.is_false(lmin) || prop.is_false(lmax))
            return false;
        if (!min_true) {
            if (!prop.push_level(lmin)) {
                ConflictResolutionHandler handler{this, cindex};
                if (!prop.resolve_conflicts(handler))
                    throw UNSATError();
                return false;
            }
            if (prop.is_false(lmax)) {
                prop.pop_level();
                return false;
            }
            max_true = prop.is_true(lmax);
            pushed_min = true;
        }
        if (!max_true) {
            if (!prop.push_level(lmax)) {
                if (pushed_min) {
                    prop.pop_level();
                    prop.pop_level();
                } else {
                    ConflictResolutionHandler handler{this, cindex};
                    if (!prop.resolve_conflicts(handler))
                        throw UNSATError();
                }
                return false;
            }
            pushed_max = true;
        }
        auto new_level = prop.get_current_level() - pushed_min - pushed_max + 1;
        auto new_lits_begin = prop.level_begin(new_level);
        auto new_lits_end = prop.get_trail().end();
        for (; new_lits_begin != new_lits_end; ++new_lits_begin) {
            Lit new_lit = *new_lits_begin;
            classes_with_literal[new_lit].push_back(cindex);
            event_handler->literal_added_to_class(cindex, new_lit);
        }
        return true;
    }

    void add_class(SharedDBPropagator&& propagator, Lit n_concrete) {
        Vertex spawner;
        auto lspawner = p_spawner_from_levels(propagator, n_concrete);
        if (lspawner)
            spawner = *lspawner;
        else {
            const Lit nclit = 2 * n_concrete;
            const auto& t = propagator.get_trail();
            for (auto beg = t.begin(), end = t.end(); beg != end; ++beg) {
                Lit l1 = *beg;
                if (l1 < nclit) {
                    for (++beg; beg != end; ++beg) {
                        Lit l2 = *beg;
                        if (l2 < nclit) {
                            spawner = std::make_pair(l1, l2);
                            beg = end;
                            break;
                        }
                    }
                    --beg;
                }
            }
            if (spawner.first > spawner.second)
                std::swap(spawner.first, spawner.second);
        }
        p_add_class(std::move(propagator), spawner);
    }

    void add_class(SharedDBPropagator&& propagator, Lit, Vertex spawner) {
        p_add_class(std::move(propagator), spawner);
    }

    SharedDBPropagator& operator[](Index idx) noexcept {
        return color_classes[idx];
    }

    const SharedDBPropagator& operator[](Index idx) const noexcept {
        return color_classes[idx];
    }

    const std::vector<SharedDBPropagator>& all() const noexcept {
        return color_classes;
    }

    ClauseDB& clauses() noexcept { return *all_clauses; }
    const ClauseDB& clauses() const noexcept { return *all_clauses; }

    Index num_classes() const noexcept { return color_classes.size(); }

    /**
     * @brief Create a new propagator at level 0
     *        that includes all clauses learned so far.
     *        More efficient than creating a propagator
     *        via its constructor (copy-constructs instead).
     * @return SharedDBPropagator
     */
    SharedDBPropagator make_empty_class() {
        assert(empty_class.get_current_level() == 0);
        empty_class.incorporate_new_clauses_at_level_0();
        return empty_class;
    }

    /**
     * @brief Spawn a new color class from the given pair.
     *
     * @param spawner
     * @return true If a new color class was successfully spawned.
     *              This guarantees that the given pair is actually feasible
     *              in the strong sense, i.e., it can be extended to a complete
     * assignment.
     * @return false The given feature literal pair is infeasible.
     *               This is also reported to the event handler.
     */
    bool spawn_class(Vertex spawner) {
        assert(empty_class.get_current_level() == 0);
        empty_class.incorporate_new_clauses_at_level_0();
        if (p_spawn_class_one_is_false(spawner))
            return false;
        SharedDBPropagator new_class(empty_class);
        bool ft = new_class.is_true(spawner.first);
        bool st = new_class.is_true(spawner.second);
        if (ft | st) {
            if (ft & st) {
                event_handler->pair_definitely_feasible(spawner.first,
                                                        spawner.second);
                p_add_class(std::move(new_class), spawner);
                return true;
            }
            if (!new_class.push_level(ft ? spawner.second : spawner.first)) {
                event_handler->feature_literal_pair_infeasible(spawner.first,
                                                               spawner.second);
                new_class.resolve_or_throw();
                return false;
            }
        } else {
            if (!new_class.push_level(spawner.first)) {
                event_handler->feature_literal_infeasible(spawner.first);
                new_class.resolve_or_throw();
                return false;
            }
            if (new_class.is_false(spawner.second)) {
                event_handler->feature_literal_pair_infeasible(spawner.first,
                                                               spawner.second);
                return false;
            }
            if (!new_class.is_true(spawner.second) &&
                !new_class.push_level(spawner.second))
            {
                event_handler->feature_literal_pair_infeasible(spawner.first,
                                                               spawner.second);
                new_class.resolve_or_throw();
                return false;
            }
        }
        if (!event_handler->is_pair_known_feasible(spawner.first,
                                                   spawner.second))
        {
            if (!p_try_completion(new_class)) {
                event_handler->feature_literal_pair_infeasible(spawner.first,
                                                               spawner.second);
                return false;
            }
            event_handler->pair_definitely_feasible(spawner.first,
                                                    spawner.second);
        }
        p_add_class(std::move(new_class), spawner);
        return true;
    }

    void initialize_feasibilities() {
        empty_class.incorporate_or_throw();
        const auto nl = 2 * all_clauses->num_vars();
        for (Lit l = 0; l < nl; ++l) {
            if (empty_class.is_true(l))
                continue;
            if (empty_class.is_false(l)) {
                event_handler->feature_literal_infeasible(l);
                continue;
            }
            if (!empty_class.push_level(l)) {
                event_handler->feature_literal_infeasible(l);
                empty_class.resolve_or_throw();
                continue;
            }
            const auto& trail = empty_class.get_trail();
            for (std::size_t ti = 1, s = trail.size(); ti < s; ++ti) {
                Lit f = lit::negate(trail[ti]);
                event_handler->feature_literal_pair_infeasible(l, f);
            }
            empty_class.pop_level();
        }
        empty_class.incorporate_or_throw();
    }

    std::size_t total_memory_usage() const noexcept {
        std::size_t total_prop_usage = 0;
        for (const SharedDBPropagator& c : color_classes) {
            total_prop_usage += c.total_memory_usage();
        }
        std::size_t total_with_literal_usage = 0;
        for (const auto& v : classes_with_literal) {
            total_with_literal_usage += v.capacity() * sizeof(Index);
        }
        return empty_class.total_memory_usage() + total_prop_usage +
               (color_classes.capacity() - color_classes.size()) *
                   sizeof(SharedDBPropagator) +
               classes_with_literal.capacity() * sizeof(std::vector<Index>) +
               total_with_literal_usage +
               class_spawners.capacity() * sizeof(Vertex) +
               m_phases.capacity() / CHAR_BIT +
               m_decisions.capacity() * sizeof(Lit) + sizeof(ColorClasses);
    }

    /* void print_memory_stats() const {
        std::size_t total_prop_usage = 0;
        std::size_t max_idx = 0, i = 0, max_usage = 0;
        for (const SharedDBPropagator& c : color_classes) {
            std::size_t usage = c.total_memory_usage();
            total_prop_usage += usage;
            if (usage > max_usage) {
                max_usage = usage;
                max_idx = i;
            }
            ++i;
        }
        std::size_t total_with_literal_usage =
            classes_with_literal.capacity() * sizeof(std::vector<Index>);
        for (const auto& v : classes_with_literal) {
            total_with_literal_usage += v.capacity() * sizeof(Index);
        }
        std::cout << "    Propagators " << mibibytes(total_prop_usage) << " ("
                  << (mibibytes(total_prop_usage) / color_classes.size())
                  << " MiB per propagator)\n";
        std::cout << "    Largest propagator: \n";
        color_classes[max_idx].print_memory_stats();
        std::cout << "    Empty class "
                  << mibibytes(empty_class.total_memory_usage()) << " MiB\n";
        std::cout << "    CWL lists   " << mibibytes(total_with_literal_usage)
                  << " MiB\n";
    } */

    void reset_coloring() {
        empty_class.incorporate_new_clauses_at_level_0();
        color_classes.clear();
        for (auto& cc : classes_with_literal) {
            cc.clear();
        }
        class_spawners.clear();
    }

    template <typename FullAssignment>
    void incorporate_colors(const std::vector<FullAssignment>& colors,
                            const std::vector<Vertex>& spawners) {
        for (std::size_t i = 0, s = colors.size(); i < s; ++i) {
            SharedDBPropagator propagator{empty_class};
            propagator.incorporate_assignment(colors[i]);
            p_add_class(std::move(propagator), spawners[i]);
        }
    }
};

struct QueueVertexEntry {
    Lit lmin, lmax;
    std::uint32_t class_count;
    Bitset available_classes;

    QueueVertexEntry(Lit lmin, Lit lmax, std::uint32_t class_count,
                     Bitset available_classes)
        : lmin(lmin), lmax(lmax), class_count(class_count),
          available_classes(std::move(available_classes)) {}

    bool class_unavailable(Index cls) noexcept {
        if (available_classes[cls]) {
            available_classes[cls] = false;
            --class_count;
            return true;
        }
        return false;
    }

    bool class_available(Index cls) noexcept {
        if (!available_classes[cls]) {
            available_classes[cls] = true;
            ++class_count;
            return true;
        }
        return false;
    }

    std::size_t total_memory_usage() const noexcept {
        return sizeof(QueueVertexEntry) + available_classes.bytes_used();
    }
};

struct EntryComesBefore {
    bool operator()(const QueueVertexEntry& e1,
                    const QueueVertexEntry& e2) const noexcept {
        return e1.class_count < e2.class_count;
    }
};

template <typename IndexMap, typename Compare = EntryComesBefore>
class IndexedVertexHeap {
  public:
    QueueVertexEntry& operator[](std::size_t s) noexcept {
        assert(s < m_entries.size());
        return m_entries[s];
    }

    const QueueVertexEntry& operator[](std::size_t s) const noexcept {
        assert(s < m_entries.size());
        return m_entries[s];
    }

    Index index_of(Lit l1, Lit l2) const noexcept {
        return ordered_index_of((std::min)(l1, l2), (std::max)(l1, l2));
    }

    Index ordered_index_of(Lit lmin, Lit lmax) const noexcept {
        return m_indices.index_of(lmin, lmax);
    }

    std::vector<QueueVertexEntry>& entries() noexcept { return m_entries; }

    const std::vector<QueueVertexEntry>& entries() const noexcept {
        return m_entries;
    }

    QueueVertexEntry& top() noexcept { return m_entries.front(); }

    const QueueVertexEntry& top() const noexcept { return m_entries.front(); }

    bool empty() const noexcept { return m_entries.empty(); }

    void pop() {
        assert(!empty());
        p_swap(0, m_entries.size() - 1);
        p_pop_back();
        lowered_priority(0); // works for empty queue as well
    }

    std::size_t size() const noexcept { return m_entries.size(); }

    template <typename... Args> void push(Args&&... args) {
        Index idx = m_entries.size();
        m_entries.emplace_back(std::forward<Args>(args)...);
        m_indices.add(m_entries.back().lmin, m_entries.back().lmax, idx);
        increased_priority(idx);
    }

    void lowered_priority(std::size_t idx) noexcept {
        QueueVertexEntry* cur = &m_entries[idx];
        std::size_t cpos = idx, hsize = m_entries.size();
        while ((cpos = (cpos << 1) + 1) < hsize) {
            QueueVertexEntry* c1 = &m_entries[cpos - 1];
            QueueVertexEntry* c2 = &m_entries[cpos];
            if (m_comp(*c1, *c2)) {
                if (m_comp(*c1, *cur)) {
                    p_swap(cpos - 1, idx);
                    idx = cpos - 1;
                    cur = c1;
                } else {
                    return;
                }
            } else {
                if (m_comp(*c2, *cur)) {
                    p_swap(cpos, idx);
                    idx = cpos;
                    cur = c2;
                } else {
                    return;
                }
            }
        }
        if (cpos == hsize) {
            QueueVertexEntry* c1 = &m_entries[cpos - 1];
            if (m_comp(*c1, *cur)) {
                p_swap(idx, cpos - 1);
            }
        }
    }

    void increased_priority(std::size_t idx) noexcept {
        QueueVertexEntry* cur = &m_entries[idx];
        for (;;) {
            std::size_t par = idx >> 1;
            QueueVertexEntry* p = &m_entries[par];
            if (m_comp(*cur, *p)) {
                p_swap(par, idx);
                idx = par;
                cur = p;
            } else {
                return;
            }
        }
    }

    void clear() {
        m_entries.clear();
        m_indices.clear();
    }

    template <typename... CompareArgs>
    IndexedVertexHeap(std::size_t n_concrete, CompareArgs&&... c)
        : m_entries(), m_indices(n_concrete),
          m_comp(std::forward<CompareArgs>(c)...) {}

    std::size_t total_memory_usage() const noexcept {
        std::size_t total_queue_usage = 0;
        for (const auto& e : m_entries)
            total_queue_usage += e.total_memory_usage();
        return m_indices.total_memory_usage() + total_queue_usage +
               (m_entries.capacity() - m_entries.size()) *
                   sizeof(QueueVertexEntry);
    }

  private:
    std::vector<QueueVertexEntry> m_entries;
    IndexMap m_indices;
    Compare m_comp;

    void p_pop_back() noexcept {
        const auto& b = m_entries.back();
        m_indices.remove(b.lmin, b.lmax);
        m_entries.pop_back();
    }

    void p_swap(std::size_t i1, std::size_t i2) noexcept {
        auto& e1 = m_entries[i1];
        auto& e2 = m_entries[i2];
        m_indices.swap(e1.lmin, e1.lmax, e2.lmin, e2.lmax);
        std::swap(e1, e2);
    }
};

template <typename T> class LiteralPairMatrix {
  public:
    explicit LiteralPairMatrix(std::size_t n_concrete, const T& init)
        : row_length(n_concrete * 2),
          buffer(std::make_unique<T[]>(row_length * (n_concrete - 1))) {
        set_all(init);
    }

    void set_all(const T& value) noexcept {
        std::size_t n_concrete = row_length / 2;
        std::fill_n(buffer.get(), row_length * (n_concrete - 1), value);
    }

    std::size_t offset(Lit lmin, Lit lmax) const noexcept {
        std::size_t pairs_below = lmin >> 1;
        std::size_t row_offset =
            2 * pairs_below * (row_length - pairs_below - 1);
        std::size_t even_odd_offset = -(lmin & 1) & (row_length - lmin);
        return row_offset + even_odd_offset + lmax - lmin - 2;
    }

    T& operator()(Lit lmin, Lit lmax) noexcept {
        return buffer[offset(lmin, lmax)];
    }

    const T& operator()(Lit lmin, Lit lmax) const noexcept {
        return buffer[offset(lmin, lmax)];
    }

    template <typename Callable /*(Lit lmin, Lit lmax, T&)*/>
    void iterate_all(Callable&& callable) {
        T* current = buffer.get();
        for (Lit vmin = 0, nc = row_length / 2; vmin < nc; ++vmin) {
            Lit lmin = 2 * vmin;
            for (Lit vmax = vmin + 1; vmax < nc; ++vmax, current += 2) {
                Lit lmax = 2 * vmax;
                callable(lmin, lmax, current[0]);
                callable(lmin, lmax + 1, current[1]);
            }
            ++lmin;
            for (Lit vmax = vmin + 1; vmax < nc; ++vmax, current += 2) {
                Lit lmax = 2 * vmax;
                callable(lmin, lmax, current[0]);
                callable(lmin, lmax + 1, current[1]);
            }
        }
    }

    /**
     * @brief Iterates the 'rows' of the matrix.
     *        Each row corresponds to the entries
     *        for some literal row_lit and has two entries
     *        for each variable strictly greater than the
     *        one of row_lit.
     *
     * @tparam RowCallbacks
     * @param row_out
     */
    template <typename RowCallback /*(Lit row_lit, Lit first_other, T* row)*/>
    void iterate_rows(RowCallback&& row_out) {
        T* current = buffer.get();
        for (Lit vmin = 0, nc = row_length / 2; vmin < nc; ++vmin) {
            Lit lmin = 2 * vmin;
            std::forward<RowCallback>(row_out)(lmin, lmin + 2, current);
            current += 2 * (nc - vmin - 1);
            std::forward<RowCallback>(row_out)(lmin + 1, lmin + 2, current);
            current += 2 * (nc - vmin - 1);
        }
    }

    template <typename Callable /*(Lit lmin, Lit lmax, const T&)*/>
    void iterate_all(Callable&& callable) const noexcept {
        const T* current = buffer.get();
        for (Lit vmin = 0, nc = row_length / 2; vmin < nc; ++vmin) {
            Lit lmin = 2 * vmin;
            for (Lit vmax = vmin + 1; vmax < nc; ++vmax, current += 2) {
                Lit lmax = 2 * vmax;
                callable(lmin, lmax, current[0]);
                callable(lmin, lmax + 1, current[1]);
            }
            ++lmin;
            for (Lit vmax = vmin + 1; vmax < nc; ++vmax, current += 2) {
                Lit lmax = 2 * vmax;
                callable(lmin, lmax, current[0]);
                callable(lmin, lmax + 1, current[1]);
            }
        }
    }

    std::size_t total_memory_usage() const noexcept {
        const std::size_t n_concrete = row_length / 2;
        return sizeof(LiteralPairMatrix) +
               (row_length * (n_concrete - 1) * sizeof(T));
    }

  private:
    std::size_t row_length;
    std::unique_ptr<T[]> buffer;
};

class MatrixIndexMap {
  public:
    explicit MatrixIndexMap(std::size_t n_concrete)
        : m_matrix(n_concrete, NIL) {}

    void swap(Lit lmin1, Lit lmax1, Lit lmin2, Lit lmax2) {
        std::swap(m_matrix(lmin1, lmax1), m_matrix(lmin2, lmax2));
    }

    void remove(Lit lmin, Lit lmax) { m_matrix(lmin, lmax) = NIL; }

    void add(Lit lmin, Lit lmax, Index index) { m_matrix(lmin, lmax) = index; }

    Index index_of(Lit lmin, Lit lmax) const noexcept {
        return m_matrix(lmin, lmax);
    }

    std::size_t total_memory_usage() const noexcept {
        return m_matrix.total_memory_usage();
    }

    void clear() noexcept { m_matrix.set_all(NIL); }

  private:
    LiteralPairMatrix<Index> m_matrix;
};

class ColoringHeuristicSolver {
  public:
    ColoringHeuristicSolver(ClauseDB* all_clauses, Lit n_concrete,
                            ThreadGroup<void>* thread_pool)
        : all_clauses(all_clauses), n_concrete(n_concrete),
          color_classes(this, all_clauses),
          internal_inf_map(std::in_place, n_concrete),
          inf_map(&*internal_inf_map),
          vertex_queue(n_concrete, EntryComesBefore{}),
          class_completer(n_concrete, all_clauses->num_vars(), this),
          explicit_partners_of(2 * n_concrete, std::vector<Lit>{}),
          m_old_true(2 * all_clauses->num_vars(), false),
          m_new_true(2 * all_clauses->num_vars(), false),
          m_buffer(2 * all_clauses->num_vars(), false),
          m_bitset_buffer(thread_pool) {}

    ColoringHeuristicSolver(ClauseDB* all_clauses, Lit n_concrete,
                            ThreadGroup<void>* thread_pool,
                            PairInfeasibilityMap* infmap)
        : all_clauses(all_clauses), n_concrete(n_concrete),
          color_classes(this, all_clauses), internal_inf_map(std::nullopt),
          inf_map(infmap), vertex_queue(n_concrete, EntryComesBefore{}),
          class_completer(n_concrete, all_clauses->num_vars(), this),
          explicit_partners_of(2 * n_concrete, std::vector<Lit>{}),
          m_old_true(2 * all_clauses->num_vars(), false),
          m_new_true(2 * all_clauses->num_vars(), false),
          m_buffer(2 * all_clauses->num_vars(), false),
          m_bitset_buffer(thread_pool) {}

    void initialize_feasibilities() {
        color_classes.initialize_feasibilities();
    }

    template <typename RNG>
    Vertex random_from_class(Index index, RNG& rng) const {
        const auto& t = color_classes[index].get_trail();
        const auto nclit = 2 * n_concrete;
        Lit l1, l2;
        std::uniform_int_distribution<std::size_t> idx_gen(0, t.size() - 1);
        do {
            l1 = t[idx_gen(rng)];
        } while (l1 >= nclit);
        do {
            l2 = t[idx_gen(rng)];
        } while (l2 >= nclit || l1 == l2);
        return {std::min(l1, l2), std::max(l1, l2)};
    }

    void extract_feasibilities() {
        if (extracted_feasibilities)
            return;

#ifdef SAMMY_CUDA_SUPPORTED
        if (n_concrete > 2048 && color_classes.all().size() > 512 &&
            should_use_cuda())
        {
            inf_map->cuda_incorporate_complete_classes(
                literals_in_class, m_bitset_buffer.thread_group());
            extracted_feasibilities = true;
            return;
        }
#endif

        inf_map->incorporate_complete_classes(literals_in_class,
                                              m_bitset_buffer.thread_group());
        extracted_feasibilities = true;
    }

    std::vector<Vertex> extract_coloring_order() {
        if (!extracted_feasibilities)
            extract_feasibilities();
        coloring_order.erase(
            std::remove_if(
                coloring_order.begin(), coloring_order.end(),
                [this](Vertex v) { return (*inf_map)(v.first, v.second); }),
            coloring_order.end());
        EdgeSet elements;
        coloring_order.erase(
            std::remove_if(
                coloring_order.begin(), coloring_order.end(),
                [&](Vertex v) { return !elements.insert(v).second; }),
            coloring_order.end());
        return coloring_order;
    }

    void feature_literal_infeasible(Lit l) {
        if (l < 2 * n_concrete) {
            inf_map->literal_infeasible(l);
        }
    }

    void feature_literal_pair_infeasible(Lit l1, Lit l2) {
        if ((std::max)(l1, l2) < 2 * n_concrete) {
            inf_map->literal_pair_infeasible(l1, l2);
        }
    }

    bool is_covered(Lit lmin, Lit lmax) const noexcept {
        if ((*inf_map)(lmin, lmax))
            return true;
        for (Index cls_index : color_classes.with_literal(lmin)) {
            if (literals_in_class[cls_index][lmax])
                return true;
        }
        return false;
    }

    Index find_covering_class(Lit lmin, Lit lmax) const {
        if ((*inf_map)(lmin, lmax))
            return NIL;
        for (Index cls_index : color_classes.with_literal(lmin)) {
            if (literals_in_class[cls_index][lmax])
                return cls_index;
        }
        throw std::out_of_range(
            "Vertex (" + std::to_string(lmin) + ", " + std::to_string(lmax) +
            ") is not covered and not marked as infeasible!");
    }

    void find_color_class_for_first_prefer_half(Lit lmin, Lit lmax) {
        // check if already colored or infeasible
        if (is_covered(lmin, lmax))
            return;
        // need to copy (concurrent modification!)
        candidates = color_classes.with_literal(lmin);
        // check classes containing lmin
        if (p_find_first_suitable_candidate_one_new(lmin, lmax, lmax))
            return;
        // check classes containing lmax
        candidates = color_classes.with_literal(lmax);
        if (p_find_first_suitable_candidate_one_new(lmin, lmax, lmin))
            return;
        // check all classes
        if (p_find_first_suitable_candidate_from_all(lmin, lmax))
            return;
        // if all else fails, create new class
        if (color_classes.spawn_class(Vertex(lmin, lmax))) {
            coloring_order.emplace_back(lmin, lmax);
        }
    }

    template <typename Callable /*(Lit lmin, Lit lmax)*/>
    void iterate_all_uncovered(Callable&& callable) const {
        Lit nclit = 2 * n_concrete;

#ifdef SAMMY_CUDA_SUPPORTED
        if (nclit > 2048 && color_classes.all().size() > 512 &&
            should_use_cuda())
        {
            try {
                cuda_iterate_all_uncovered(
                    literals_in_class, color_classes.get_classes_with_literal(),
                    inf_map, callable);
                return;
            } catch (const CUDAError& err) {
                std::cerr << "Not using CUDA because of an error: "
                          << err.what() << "\n";
                had_cuda_error(err);
            }
        }
#endif

        Bitset row_uncolored(nclit, true);
        for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
            row_uncolored.set();
            row_uncolored ^= (*inf_map)[lmin];
            const auto& with_lit = color_classes.with_literal(lmin);
            sammy::bitwise_filter(
                m_bitset_buffer, row_uncolored,
                p_make_bitset_transform_iterator(with_lit.begin()),
                p_make_bitset_transform_iterator(with_lit.end()));
            for (Lit lmax : row_uncolored.ones_from(lmin + 1)) {
                (std::forward<Callable>(callable))(lmin, lmax);
            }
        }
    }

    template <typename Callable /*(Lit lmin, Lit lmax)*/>
    void iterate_all_uniquely_covered(Callable&& callable) {
        Lit nclit = 2 * n_concrete;
        Bitset multi_covered(nclit, false);
        Bitset prev_covered(nclit, false);
        Bitset tmp(nclit, false);
        for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
            multi_covered.reset();
            prev_covered.reset();
            const auto& with_lit = color_classes.with_literal(lmin);
            for (Index cci : with_lit) {
                const auto& ccbits = literals_in_class[cci];
                tmp = prev_covered;
                tmp &= ccbits;
                multi_covered |= tmp;
                prev_covered |= ccbits;
            }
            tmp = multi_covered;
            tmp.flip();
            tmp -= (*inf_map)[lmin];
            for (Lit lmax : tmp.ones_from(lmin + 1)) {
                (std::forward<Callable>(callable))(lmin, lmax);
            }
        }
    }

    template <typename FullAssignment>
    void replace_colors(const std::vector<FullAssignment>& colors,
                        const std::vector<Vertex>& spawners) {
        reset_coloring();
        color_classes.incorporate_colors(colors, spawners);
    }

    PartialSolution get_partial_solution() const {
        return PartialSolution{all_clauses->num_vars(), inf_map,
                               all_classes().begin(), all_classes().end()};
    }

    LiteralPairMatrix<std::uint16_t> saturated_coverage_count_matrix() const {
        std::vector<std::uint16_t> coverage_buf(2 * all_clauses->num_vars(), 0);
        LiteralPairMatrix<std::uint16_t> result(n_concrete, 0);
        std::uint16_t* cv_buf_beg = &coverage_buf[2];
        std::uint16_t* cv_buf_end = coverage_buf.data() + (2 * n_concrete);
        result.iterate_rows([&](Lit row_lit, Lit /*next*/, std::uint16_t* row) {
            if (color_classes.all().size() <=
                std::numeric_limits<std::uint16_t>::max())
            {
                for (Index cci : color_classes.with_literal(row_lit)) {
                    const SharedDBPropagator& cc = color_classes[cci];
                    for (Lit l : cc.get_trail())
                        ++coverage_buf[l];
                }
            } else {
                for (Index cci : color_classes.with_literal(row_lit)) {
                    const SharedDBPropagator& cc = color_classes[cci];
                    for (Lit l : cc.get_trail())
                        p_saturated_add(coverage_buf[l], 1);
                }
            }
            for (std::uint16_t* cv_cur = cv_buf_beg; cv_cur != cv_buf_end;
                 ++cv_cur, ++row)
            {
                p_saturated_add(*row, *cv_cur);
                *cv_cur = 0;
            }
            if (row_lit & 1)
                cv_buf_beg += 2;
        });
        return result;
    }

    void find_color_class_for(Lit lmin, Lit lmax) {
        find_color_class_for_first_prefer_half(lmin, lmax);
    }

    bool complete_class(Index index) {
        SharedDBPropagator& prop = color_classes[index];
        if (prop.get_trail().size() == all_clauses->num_vars())
            return true;
        p_prop_to_buffer(prop, m_old_true);
        bool compl_result = class_completer.complete_class(prop);
        p_prop_to_buffer(prop, m_new_true);
        color_classes.class_changed(index, m_old_true, m_new_true, m_buffer);
        return compl_result;
    }

    bool complete_classes() {
        Index nc = color_classes.all().size();
        for (Index i = 0; i < nc; ++i) {
            if (!complete_class(i)) {
                if (!quiet) {
                    std::cout << "Could not complete class " << i
                              << " without rolling back a decision!\n";
                }
                return false;
            }
        }
        return true;
    }

    void add_greedy_false_class() {
        SharedDBPropagator prop = color_classes.make_empty_class();
        GreedySAT gsat(&prop, n_concrete, PreferFalse{});
        if (gsat.solve()) {
            add_color_class(std::move(prop));
        }
    }

    void color_nonlazy(bool use_initialization = true) {
        if (use_initialization) {
            add_greedy_false_class();
            auto initial = default_initial_vertices();
            p_bulk_add_vertices(initial.cbegin(), initial.cend());
            color_vertices_in_queue();
        }
        for (;;) {
            add_all_uncolored_to_queue();
            color_vertices_in_queue();
            if (everything_covered()) {
                if (complete_classes())
                    return;
            }
        }
    }

    void color_nonlazy(const std::vector<Vertex>& initial_vertices) {
        for (Vertex v : initial_vertices) {
            add_vertex_to_queue(v.first, v.second);
        }
        color_vertices_in_queue();
        for (;;) {
            add_all_uncolored_to_queue();
            color_vertices_in_queue();
            if (everything_covered()) {
                if (complete_classes())
                    return;
            }
        }
    }

    void color_lazy(const std::vector<Vertex>& initial_vertices_ = {}) {
        bool had_initial_vertices = true;
        if (initial_vertices_.empty() && coloring_order.empty()) {
            add_greedy_false_class();
            had_initial_vertices = false;
        }
        const std::vector<Vertex>& initial_vertices =
            !had_initial_vertices ? default_initial_vertices()
                                  : initial_vertices_;
        for (Vertex v : initial_vertices) {
            add_vertex_to_queue(v.first, v.second);
        }
        color_vertices_in_queue();
        std::size_t lazy_col_explicit = coloring_order.size();
        std::size_t lazy_col_next = p_lazy_color_next(lazy_col_explicit);
        if (had_initial_vertices &&
            lazy_col_explicit != color_classes.all().size())
        {
            for (Vertex v : initial_vertices) {
                if (!is_covered(v.first, v.second)) {
                    // may happen due to conflict resolution
                    find_color_class_for(v.first, v.second);
                }
            }
        }
        if (!quiet) {
            std::cout << "Colored initial " << lazy_col_explicit
                      << " vertices (" << color_classes.all().size()
                      << " classes)\n";
        }
        // print_memory_stats();
        for (;;) {
            std::size_t n_unc = num_uncovered();
            if (!quiet) {
                std::cout << "Uncovered: " << n_unc << "...\n";
            }
            if (n_unc == 0) {
                if (complete_classes()) {
                    return;
                }
                continue;
            }
            if (n_unc <= 2 * lazy_col_next) {
                add_all_uncolored_to_queue();
                lazy_col_explicit += n_unc;
            } else {
                lazy_col_explicit += lazy_sample_to_queue(lazy_col_next, n_unc);
            }
            // print_memory_stats();
            color_vertices_in_queue();
            lazy_col_next = p_lazy_color_next(lazy_col_explicit);
            if (!quiet) {
                std::cout << "Colored " << lazy_col_explicit
                          << " explicit vertices ("
                          << color_classes.all().size() << " classes) ...\n";
            }
        }
    }

    std::size_t lazy_sample_to_queue(std::size_t approx_vertices,
                                     std::size_t n_unc) {
        double sample_rate = double(approx_vertices) / n_unc;
        auto& rng = sammy::rng();
        std::geometric_distribution<std::size_t> skip_sample_dist(sample_rate);
        std::size_t skip = skip_sample_dist(rng);
        std::size_t count = 0;
        std::vector<Vertex> vbuf;
        iterate_all_uncovered([&](Lit lmin, Lit lmax) {
            if (!skip) {
                skip = skip_sample_dist(rng);
                vbuf.emplace_back(lmin, lmax);
                ++count;
            } else {
                --skip;
            }
        });
        p_bulk_add_vertices(vbuf.begin(), vbuf.end());
        return count;
    }

    std::vector<Vertex> default_initial_vertices() {
        std::vector<Vertex> result;
        result.reserve(n_concrete);
        std::vector<Lit> v2;
        v2.reserve(n_concrete);
        for (Lit l = 0; l < n_concrete; ++l) {
            v2.push_back(lit::positive_lit(l));
        }
        std::shuffle(v2.begin(), v2.end(), sammy::rng());
        SharedDBPropagator prop = color_classes.make_empty_class();
        for (Lit v1 = 0; v1 < n_concrete; ++v1) {
            Lit l1 = lit::positive_lit(v1);
            if (!prop.is_open(l1))
                continue;
            if (!prop.push_level(l1)) {
                prop.resolve_or_throw();
                continue;
            }
            auto start_from = v2.begin() + v1;
            for (; start_from != v2.end(); ++start_from) {
                Lit o = *start_from;
                if (!prop.is_open(o))
                    continue;
                bool pres = prop.push_level(o);
                prop.pop_level();
                if (!pres)
                    continue;
                break;
            }
            if (start_from == v2.end()) {
                for (start_from = v2.begin(); start_from != v2.end();
                     ++start_from)
                {
                    Lit o = *start_from;
                    if (!prop.is_open(o))
                        continue;
                    bool pres = prop.push_level(o);
                    prop.pop_level();
                    if (!pres)
                        continue;
                    break;
                }
            }
            if (start_from != v2.end()) {
                result.emplace_back(std::min(l1, *start_from),
                                    std::max(l1, *start_from));
            }
            prop.pop_level();
        }
        return result;
    }

    const std::vector<Vertex>& class_spawners() const noexcept {
        return color_classes.spawners();
    }

    void new_color_class(Index index) {
        for (auto& e : vertex_queue.entries()) {
            e.class_count += 1;
            e.available_classes.push_back(true);
        }
        literals_in_class.push_back(Bitset(2 * n_concrete, false));
        const auto& cc = color_classes[index];
        for (Lit lpos : cc.get_trail()) {
            literal_added_to_class(index, lpos);
        }
    }

    void literal_added_to_class(Index cindex, Lit l) {
        if (l >= 2 * n_concrete)
            return;
        literals_in_class[cindex][l] = true;
        l = lit::negate(l);
        auto& partners = explicit_partners_of[l];
        auto pend =
            std::remove_if(partners.begin(), partners.end(), [&](Lit o) {
                Index idx = vertex_queue.index_of(l, o);
                if (idx == NIL)
                    return true;
                if (vertex_queue[idx].class_unavailable(cindex)) {
                    vertex_queue.increased_priority(idx);
                }
                return false;
            });
        partners.erase(pend, partners.end());
    }

    void literal_removed_from_class(Index cindex, Lit l) {
        if (l >= 2 * n_concrete)
            return;
        literals_in_class[cindex][l] = false;
        SharedDBPropagator& cc = color_classes[cindex];
        l = lit::negate(l);
        auto& partners = explicit_partners_of[l];
        auto pend =
            std::remove_if(partners.begin(), partners.end(), [&](Lit o) {
                Index idx = vertex_queue.index_of(l, o);
                if (idx == NIL)
                    return true;
                if (!cc.is_false(o) &&
                    vertex_queue[idx].class_available(cindex)) {
                    vertex_queue.lowered_priority(idx);
                }
                return false;
            });
        partners.erase(pend, partners.end());
    }

    bool is_pair_known_feasible(Lit lmin, Lit lmax) const noexcept {
        return inf_map->is_definitely_feasible(lmin, lmax);
    }

    void pair_definitely_feasible(Lit lmin, Lit lmax) noexcept {
        inf_map->set_definitely_feasible(lmin, lmax);
    }

    void color_vertices_in_queue() {
        while (!vertex_queue.empty()) {
            QueueVertexEntry ve = std::move(vertex_queue.top());
            vertex_queue.pop();
            if (ve.class_count == 0) {
                if (color_classes.spawn_class(Vertex(ve.lmin, ve.lmax))) {
                    coloring_order.emplace_back(ve.lmin, ve.lmax);
                }
            } else {
                find_color_class_for(ve.lmin, ve.lmax);
            }
        }
    }

    void initialize_with_start_vertices(const std::vector<Vertex>& vertices) {
        for (Vertex v : vertices) {
            find_color_class_for(v.first, v.second);
        }
    }

    void add_vertex_to_queue(Lit lmin, Lit lmax) {
        const auto& cc_all = color_classes.all();
        auto s = cc_all.size();
        std::uint32_t count = s;
        Bitset bset(s, true);
        for (Index i = 0, s = cc_all.size(); i < s; ++i) {
            const auto& cf = cc_all[i];
            if (cf.is_false(lmin) || cf.is_false(lmax)) {
                bset[i] = false;
                --count;
            }
        }
        vertex_queue.push(lmin, lmax, count, std::move(bset));
    }

    template <typename Iterator>
    void add_vertices_to_queue(Iterator begin, Iterator end) {
        p_bulk_add_vertices(begin, end);
    }

    void add_all_uncolored_to_queue() {
        std::vector<Vertex> vbuf;
        iterate_all_uncovered(
            [&](Lit lmin, Lit lmax) { vbuf.emplace_back(lmin, lmax); });
        p_bulk_add_vertices(vbuf.begin(), vbuf.end());
    }

    void add_color_class(const SharedDBPropagator& cc) {
        SharedDBPropagator prop(cc);
        color_classes.add_class(std::move(prop), n_concrete);
    }

    void add_color_class(const SharedDBPropagator& cc, Vertex spawner) {
        SharedDBPropagator prop(cc);
        color_classes.add_class(std::move(prop), n_concrete, spawner);
    }

    void add_color_class(SharedDBPropagator&& cc, Vertex spawner) {
        color_classes.add_class(std::move(cc), n_concrete, spawner);
    }

    void add_color_class(SharedDBPropagator&& cc) {
        color_classes.add_class(std::move(cc), n_concrete);
    }

    const std::vector<SharedDBPropagator>& all_classes() const noexcept {
        return color_classes.all();
    }

    bool everything_covered() const {
        bool res = true;
        iterate_all_uncovered([&](Lit, Lit) { res = false; });
        return res;
    }

    std::size_t num_uncovered() const {
        Lit nclit = 2 * n_concrete;

#ifdef SAMMY_CUDA_SUPPORTED
        if (nclit > 2048 && color_classes.all().size() > 512 &&
            should_use_cuda())
        {
            try {
                auto result = cuda_count_uncovered(
                    literals_in_class, color_classes.get_classes_with_literal(),
                    inf_map);
                return result;
            } catch (const CUDAError& err) {
                std::cerr << "Not using CUDA because of an error: "
                          << err.what() << "\n";
                had_cuda_error(err);
            }
        }
#endif

        Bitset row_uncolored(nclit, true);
        Bitset prv_rows(nclit, true);
        std::size_t nu = 0;
        for (Lit lmin = 0; lmin < nclit; ++lmin) {
            prv_rows[lmin] = false;
            row_uncolored = prv_rows;
            row_uncolored -= (*inf_map)[lmin];
            const auto& with_lit = color_classes.with_literal(lmin);
            sammy::bitwise_filter(
                m_bitset_buffer, row_uncolored,
                p_make_bitset_transform_iterator(with_lit.begin()),
                p_make_bitset_transform_iterator(with_lit.end()));
            nu += row_uncolored.count();
        }
        return nu;
    }

    std::size_t total_memory_usage() const noexcept {
        return color_classes.total_memory_usage() +
               (internal_inf_map ? internal_inf_map->total_memory_usage() : 0) +
               vertex_queue.total_memory_usage() +
               class_completer.total_memory_usage() +
               p_mem_usage_explicit_partners_of() +
               candidates.capacity() * sizeof(Index) +
               p_mem_usage_literals_in_class(literals_in_class) +
               m_old_true.bytes_used() + m_new_true.bytes_used() +
               m_buffer.bytes_used() + sizeof(ColoringHeuristicSolver);
    }

    /**
     * @brief Reset the coloring (dropping all color classes and related
     * information). Retains quite a bit of internal memory (vector capacities
     * etc).
     */
    void reset_coloring() {
        color_classes.reset_coloring();
        vertex_queue.clear();
        for (auto& e : explicit_partners_of)
            e.clear();
        candidates.clear();
        literals_in_class.clear();
        coloring_order.clear();
    }

    /**
     * @brief Remove color classes given by the indices.
     * Deletes the coloring order information.
     * Returns the color classes that were removed.
     */
    std::vector<SharedDBPropagator>
    remove_color_classes(std::vector<std::size_t> classes) {
        if (!vertex_queue.empty()) {
            throw std::logic_error(
                "remove_color_classes called with non-empty vertex queue!");
        }
        std::for_each(explicit_partners_of.begin(), explicit_partners_of.end(),
                      [](CVec& v) { v.clear(); });
        std::sort(classes.begin(), classes.end());
        candidates.clear();
        coloring_order.clear();
        auto result = color_classes.remove_classes(classes);
        literals_in_class.erase(remove_indices(literals_in_class.begin(),
                                               literals_in_class.end(),
                                               classes.begin(), classes.end()),
                                literals_in_class.end());
        return result;
    }

    /**
     * @brief Get the resulting sample as list of lists of internal literals.
     *
     * @return std::vector<std::vector<Lit>>
     */
    std::vector<std::vector<Lit>> internal_solution(bool only_concrete) const {
        std::vector<std::vector<Lit>> result;
        for (const auto& cc : color_classes.all()) {
            result.emplace_back(cc.get_trail());
        }
        if (only_concrete) {
            const Lit nclit = 2 * n_concrete;
            for (auto& cc : result) {
                cc.erase(std::remove_if(cc.begin(), cc.end(),
                                        [&](Lit l) { return l >= nclit; }),
                         cc.end());
            }
        }
        return result;
    }

    /**
     * @brief Get the resulting sample as list of lists of external literals.
     * @return std::vector<ExternalClause>
     */
    std::vector<std::vector<ExternalLit>>
    external_solution(bool only_concrete) const {
        std::vector<std::vector<ExternalLit>> result;
        const Lit nclit = 2 * n_concrete;
        for (const auto& cc : color_classes.all()) {
            result.emplace_back();
            auto& r = result.back();
            for (const Lit l : cc.get_trail()) {
                if (!only_concrete || l < nclit) {
                    r.push_back(lit::externalize(l));
                }
            }
        }
        return result;
    }

    void set_quiet(bool q) noexcept { quiet = q; }

  private:
    struct ClassIndexToBitset {
        const Bitset& operator()(Index cls_index) const noexcept {
            return that->literals_in_class[cls_index];
        }

        const ColoringHeuristicSolver* that;
    };

    template <typename IndexIterator>
    boost::transform_iterator<ClassIndexToBitset, IndexIterator, const Bitset&>
    p_make_bitset_transform_iterator(IndexIterator iter) const {
        ClassIndexToBitset trans{this};
        return {iter, trans};
    }

    static void p_saturated_add(std::uint16_t& x, std::uint16_t v) {
        if (x > std::numeric_limits<std::uint16_t>::max() - v) {
            x = std::numeric_limits<std::uint16_t>::max();
        } else {
            x += v;
        }
    }

    std::size_t
    p_mem_usage_literals_in_class(const std::vector<Bitset>& v) const noexcept {
        return std::transform_reduce(
            v.begin(), v.end(), std::size_t(0), std::plus<>{},
            [](const Bitset& b) { return b.bytes_used(); });
    }

    std::size_t p_mem_usage_explicit_partners_of() const noexcept {
        std::size_t element_sum = 0;
        for (const auto& v : explicit_partners_of) {
            element_sum += v.capacity() * sizeof(Lit);
        }
        return explicit_partners_of.capacity() * sizeof(std::vector<Lit>) +
               element_sum;
    }

    std::size_t p_lazy_color_next(std::size_t already_colored) const {
        if (already_colored < 1'000'000) {
            if (already_colored < 256)
                return 256;
            return already_colored;
        }
        return 1'000'000;
    }

    /**
     * Much more efficient than adding all vertices individually;
     * that indeed becomes the bottleneck in color_lazy
     * if one does not use the bulk method.
     */
    template <typename Iterator>
    void p_bulk_add_vertices(Iterator begin, Iterator end) {
        std::vector<Bitset> cl_allowing_lit =
            p_compute_classes_allowing_literal();
        for (; begin != end; ++begin) {
            Vertex v = *begin;
            auto bs = make_larger_bitset(cl_allowing_lit[v.first]);
            bs &= cl_allowing_lit[v.second];
            vertex_queue.push(v.first, v.second, bs.count(), std::move(bs));
        }
    }

    std::vector<Bitset> p_compute_classes_allowing_literal() {
        const std::size_t n_cols = color_classes.all().size();
        const std::size_t nc_lit = 2 * n_concrete;
        std::vector<Bitset> classes_allowing_lit(nc_lit, Bitset(n_cols, true));
        for (std::size_t i = 0; i < n_cols; ++i) {
            const SharedDBPropagator& prop = color_classes[i];
            for (Lit l : prop.get_trail()) {
                if (l < nc_lit) {
                    l = lit::negate(l);
                    classes_allowing_lit[l][i] = false;
                }
            }
        }
        return classes_allowing_lit;
    }

    void p_prop_to_buffer(const SharedDBPropagator& prop, Bitset& buffer) {
        buffer.reset();
        for (Lit l : prop.get_trail()) {
            buffer[l] = true;
        }
    }

    bool p_find_first_suitable_candidate_one_new(Lit lmin, Lit lmax, Lit lnew) {
        for (Index cls_index : candidates) {
            auto& cls = color_classes[cls_index];
            if (!cls.is_false(lnew) &&
                color_classes.push_vertex(cls_index, lmin, lmax))
            {
                coloring_order.emplace_back(lmin, lmax);
                return true;
            }
        }
        return false;
    }

    bool p_find_first_suitable_candidate_from_all(Lit lmin, Lit lmax) {
        for (Index i = 0, nc = color_classes.num_classes(); i < nc; ++i) {
            auto& cls = color_classes[i];
            if (!cls.is_false(lmin) && !cls.is_false(lmax)) {
                if (color_classes.push_vertex(i, lmin, lmax)) {
                    coloring_order.emplace_back(lmin, lmax);
                    return true;
                }
            }
        }
        return false;
    }

    bool p_find_first_suitable_candidate(Lit lmin, Lit lmax) {
        for (Index cls_index : candidates) {
            auto& cls = color_classes[cls_index];
            if (!cls.is_false(lmin) && !cls.is_false(lmax)) {
                if (color_classes.push_vertex(cls_index, lmin, lmax)) {
                    return true;
                }
            }
        }
        return false;
    }

    ClauseDB* all_clauses;
    Lit n_concrete;
    ColorClasses<ColoringHeuristicSolver> color_classes;
    std::optional<PairInfeasibilityMap> internal_inf_map;
    PairInfeasibilityMap* inf_map;
    IndexedVertexHeap<MatrixIndexMap> vertex_queue;
    ClassCompleter<ColoringHeuristicSolver> class_completer;
    std::vector<std::vector<Lit>> explicit_partners_of;
    std::vector<Index> candidates;
    std::vector<Bitset> literals_in_class;
    Bitset m_old_true, m_new_true, m_buffer;
    mutable BitsetOperationsBuffer m_bitset_buffer;
    std::vector<Vertex> coloring_order;
    bool extracted_feasibilities = false;
    bool quiet = false;
};

} // namespace sammy

#endif
==> ./input.h <==
#ifndef SAMMY_INPUT_H_INCLUDED_
#define SAMMY_INPUT_H_INCLUDED_

#include "clause_db.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
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
    if (json_data.at("type") != "software configuration model") {
        throw std::runtime_error(
            "JSON data missing the 'software configuration model' type flag!");
    }
    auto clauses =
        json_data.at("cnf_clauses").get<std::vector<ExternalClause>>();
    auto n_all = json_data.at("num_variables").get<Var>();
    auto n_concrete = json_data.at("num_concrete_features").get<Var>();
    return InputData{ClauseDB{n_all, clauses}, n_concrete,
                     json_data.at("name").get<std::string>()};
}

} // namespace sammy

#endif
==> ./ortools_solver.h <==
#ifndef SAMMY_ORTOOLS_SOLVER_H_INCLUDED_
#define SAMMY_ORTOOLS_SOLVER_H_INCLUDED_
#ifdef SAMMY_ORTOOLS_SUPPORTED
#if SAMMY_ORTOOLS_SUPPORTED

#include <ortools/sat/cp_model.h>
#include <ortools/sat/cp_model.pb.h>
#include <ortools/sat/cp_model_solver.h>
#include <ortools/sat/model.h>
#include <ortools/sat/sat_parameters.pb.h>
#include <ortools/util/time_limit.h>

#include "clause_db.h"
#include "dynamic_bitset.h"
#include "output.h"
#include "pair_infeasibility_map.h"
#include "partial_solution.h"
#include "vertex_operations.h"

#include <atomic>

namespace sammy {

class ORToolsSolver {
  public:
    /**
     * Create an ORTools solver.
     */
    ORToolsSolver(std::vector<Vertex> considered_vertices,
                  PairInfeasibilityMap* infeasibility_map, ClauseDB& clause_db,
                  std::vector<Vertex> best_local_mes,
                  std::size_t /*best_global_mes*/, std::size_t best_global_lb,
                  std::vector<DynamicBitset> covering_assignments,
                  std::optional<DynamicBitset> remaining_configuration)
        : m_all_vertices(std::move(considered_vertices)),
          m_infeasibility_map(infeasibility_map), m_clauses(&clause_db),
          m_propagator(m_clauses), m_best_local_mes(std::move(best_local_mes)),
          m_covering_assignments(std::move(covering_assignments)),
          m_remaining_configuration(std::move(remaining_configuration)),
          m_lower_bound(m_best_local_mes.size()),
          m_initial_global_lower_bound(best_global_lb),
          m_partial_solution(m_clauses->num_vars(), m_infeasibility_map) {}

    /**
     * Turn on 'single improvement' mode, where we do not
     * minimize a objective function but simply ask for
     * an improvement by at least one configuration.
     */
    void set_find_single_improvement(bool find_single_improvement) {
        if (m_initialized_objective)
            throw std::logic_error(
                "Called set_find_single_improvement too late!");
        m_want_single_improvement = find_single_improvement;
    }

    /**
     * Abort the solve (on timeout, optimality, LNS with better solution, ...).
     */
    void abort() noexcept { m_abort_flag.store(true); }

    /**
     * The size of the graph (all vertices) we are working on.
     */
    std::size_t get_graph_size() const noexcept {
        return m_all_vertices.size();
    }

    /**
     * Enum describing the outcome of the solve process.
     */
    enum class SolveResult {
        IMPROVED_SOLUTION, //< We found a new, optimal solution.
        ABORTED,           //< We were aborted (timeout, ...) before completion.
        SOLUTION_WAS_OPTIMAL //< We completed but the initial solution was
                             // optimal.
    };

    /**
     * Run the solution process.
     * @return SolveResult describing the outcome.
     */
    SolveResult solve() {
        p_construct_model();
        if (m_abort_flag.load())
            return SolveResult::ABORTED;
        Model model = m_model.Build();
        if (m_abort_flag.load())
            return SolveResult::ABORTED;
        p_report_event("CPSAT_MODEL_BUILT", {});
        operations_research::sat::Model model_config;
        model_config.GetOrCreate<operations_research::TimeLimit>()
            ->RegisterExternalBooleanAsLimit(&m_abort_flag);
        p_report_event("CPSAT_SOLVE_BEGIN", {});
        operations_research::sat::CpSolverResponse response;
        response = operations_research::sat::SolveCpModel(model, &model_config);
        switch (response.status()) {
        case operations_research::sat::CpSolverStatus::OPTIMAL:
            p_report_event("CPSAT_SOLVE_OPTIMAL", {});
            break;

        case operations_research::sat::CpSolverStatus::FEASIBLE:
            p_report_event("CPSAT_SOLVE_FEASIBLE", {});
            break;

        default:
        case operations_research::sat::CpSolverStatus::UNKNOWN:
            p_report_event("CPSAT_SOLVER_ABORTED", {});
            return SolveResult::ABORTED;

        case operations_research::sat::CpSolverStatus::MODEL_INVALID:
            throw std::logic_error("Model invalid error!");

        case operations_research::sat::CpSolverStatus::INFEASIBLE:
            p_report_event("CPSAT_SOLVE_INFEASIBLE", {});
            m_lower_bound = m_covering_assignments.size();
            if (m_lower_bound > m_initial_global_lower_bound &&
                m_lower_bound_callback)
            {
                m_lower_bound_callback(m_lower_bound, m_all_vertices);
            }
            return SolveResult::SOLUTION_WAS_OPTIMAL;
        }
        if (!m_want_single_improvement) {
            std::size_t new_bound =
                std::size_t(response.best_objective_bound() + 0.01);
            if (new_bound > m_lower_bound) {
                m_lower_bound = new_bound;
                if (m_lower_bound > m_initial_global_lower_bound &&
                    m_lower_bound_callback)
                {
                    m_lower_bound_callback(m_lower_bound, m_all_vertices);
                }
            }
        }
        p_extract_partial_solution(response);
        return SolveResult::IMPROVED_SOLUTION;
    }

    /**
     * Get the best solution found. Only valid if
     * the SolveResult was IMPROVED_SOLUTION.
     *
     * @return The best solution found.
     */
    PartialSolution get_partial_solution() const { return m_partial_solution; }

    /**
     * @return A reference to the covering assignments given initially.
     */
    const std::vector<DynamicBitset>&
    get_covering_assignments() const noexcept {
        return m_covering_assignments;
    }

    /**
     * Check whether we may have found an improvement
     * on the global lower bound.
     */
    bool improved_global_bound() const noexcept {
        return m_lower_bound > m_initial_global_lower_bound;
    }

    /**
     * Check whether we may have found an improvement
     * on the global MES size. This method exists for
     * compatibility with the other solvers and always returns false.
     */
    bool improved_mes() const noexcept { return false; }

    /**
     * Get a copy of the best clique found.
     */
    std::vector<Vertex> get_best_mes() const { return m_best_local_mes; }

    /**
     * Get the best lower bound on the number of configurations
     * needed for the subgraph.
     */
    std::size_t get_best_bound() const noexcept { return m_lower_bound; }

    /**
     * Set the event recorder to report events to;
     * this is purely optional.
     */
    void set_event_recorder(EventRecorder* recorder) noexcept {
        m_local_recorder = recorder;
    }

    /**
     * Get the vertices inducing the subgraph on which
     * we established the lower bound.
     */
    std::vector<Vertex> get_best_bound_subgraph() const {
        return m_all_vertices;
    }

    /**
     * Set a callback to be invoked each time we try to search for a new clique.
     * It can return a clique (containing any vertices); the solver will filter
     * out all irrelevant vertices and only consider the ones in m_all_vertices.
     * This method exists for compatibility with other solvers and does nothing.
     */
    void set_clique_candidate_callback(std::function<std::vector<Vertex>()>) {}

    /**
     * Set a callback to be invoked each time this solver finds a new lower
     * bound. The callback is involved with the lower bound and the vertices of
     * the subgraph inducing it.
     */
    void set_lower_bound_callback(
        std::function<void(std::size_t, const std::vector<Vertex>&)> callback) {
        m_lower_bound_callback = std::move(callback);
    }

  private:
    /**
     * All vertices in the subproblem we are considering.
     */
    std::vector<Vertex> m_all_vertices;

    /**
     * The pair infeasibility map.
     */
    PairInfeasibilityMap* m_infeasibility_map;

    /**
     * The clause DB.
     */
    ClauseDB* m_clauses;

    /**
     * A propagator connected to m_clauses.
     */
    SharedDBPropagator m_propagator;

    /**
     * Symmetry-breaking MES.
     */
    std::vector<Vertex> m_best_local_mes;

    /**
     * The removed configurations constituting this subproblem.
     */
    std::vector<DynamicBitset> m_covering_assignments;

    /**
     * A single remaining configuration.
     * If it is available, can help breaking symmetries.
     */
    std::optional<DynamicBitset> m_remaining_configuration;

    /**
     * The best lower bound we have found.
     */
    std::size_t m_lower_bound;

    /**
     * The best lower bound we have found for the complete problem.
     */
    std::size_t m_initial_global_lower_bound;

    /**
     * Pointer to the event recorder to store our events in.
     */
    EventRecorder* m_local_recorder{nullptr};

    /**
     * A flag to control whether we do minimization
     * or just find any improvement.
     */
    bool m_want_single_improvement{false};

    /**
     * A flag to store whether we have initialized
     * the objective function.
     */
    bool m_initialized_objective{false};

    /**
     * Callback invoked with new lower bounds
     * and the subgraphs on which they were obtained.
     */
    std::function<void(std::size_t, const std::vector<Vertex>&)>
        m_lower_bound_callback;

    /**
     * A flag to abort the CP-SAT solve.
     */
    std::atomic<bool> m_abort_flag{false};

    /**
     * Key types from CP-SAT.
     */
    using ModelBuilder = operations_research::sat::CpModelBuilder;
    using Model = operations_research::sat::CpModelProto;
    using SolverResponse = operations_research::sat::CpSolverResponse;
    using BoolVar = operations_research::sat::BoolVar;
    using VarOrVal = std::variant<BoolVar, bool>;

    /**
     * The CP-SAT model (or rather, its builder).
     */
    ModelBuilder m_model;

    /**
     * The CP-SAT variables (or fixed values) for the configurations.
     * Index [c][x] has the value of variable x in configuration c.
     */
    std::vector<std::vector<VarOrVal>> m_configuration_variables;

    /**
     * The interaction coverage variables (or fixed values).
     * Index [c][i] checks for coverage of the ith interaction in configuration
     * c.
     */
    std::vector<std::vector<VarOrVal>> m_interaction_coverage_variables;

    /**
     * Variables for the usage/activity of configurations.
     */
    std::vector<VarOrVal> m_configuration_usage_variables;

    /**
     * Buffer for clause building.
     */
    std::vector<BoolVar> m_clause_buffer;

    /**
     * Partial solution (initially empty).
     */
    PartialSolution m_partial_solution;

    /**
     * Construct the CP-SAT model.
     */
    void p_construct_model() {
        std::size_t num_copies_needed = m_covering_assignments.size() - 1;
        if (!m_want_single_improvement) {
            // turn on single improvement if we are in a -1 situation
            if (num_copies_needed == m_lower_bound) {
                m_want_single_improvement = true;
            }
        }
        m_initialized_objective = true;
        p_report_event(
            "BEGIN_CONSTRUCT_CPSAT_MODEL",
            {{"max_num_configurations", num_copies_needed},
             {"best_lb", m_lower_bound},
             {"best_global_lb", m_initial_global_lower_bound},
             {"best_local_mes_size", m_best_local_mes.size()},
             {"num_interactions", m_all_vertices.size()},
             {"single_improvement_search", m_want_single_improvement}},
            "max_num_configurations", "num_interactions",
            "single_improvement_search", "best_lb", "best_global_lb",
            "best_local_mes_size");
        for (std::size_t copy_index : range(num_copies_needed)) {
            p_add_copy(copy_index);
        }
        p_add_coverage_constraints();
        p_add_objective();
        p_report_event("CPSAT_MODEL_BUILDER_POPULATED", {});
    }

    /**
     * Extract a partial solution from the CP-SAT solution.
     */
    void p_extract_partial_solution(const SolverResponse& response) {
        std::vector<DynamicBitset> new_configurations;
        std::size_t num_vars = m_clauses->num_vars();
        std::size_t num_copies = m_covering_assignments.size() - 1;
        for (std::size_t copy_index : range(num_copies)) {
            VarOrVal usage_var = m_configuration_usage_variables[copy_index];
            if (!p_get_value(usage_var, response))
                break;
            DynamicBitset next_configuration(num_vars, false);
            for (std::size_t v : range(num_vars)) {
                VarOrVal entry = m_configuration_variables[copy_index][v];
                if (p_get_value(entry, response)) {
                    next_configuration[v].set();
                }
            }
            new_configurations.emplace_back(std::move(next_configuration));
        }
        m_partial_solution = PartialSolution(Var(num_vars), m_infeasibility_map,
                                             new_configurations.begin(),
                                             new_configurations.end());
    }

    /**
     * Construct another configuration copy.
     */
    void p_add_copy(std::size_t copy_index) {
        p_add_copy_vars(copy_index);
        p_copy_add_clauses(copy_index);
        p_copy_add_usage(copy_index);
        p_copy_add_coverage(copy_index);
    }

    /**
     * Add configuration copy usage variables.
     */
    void p_copy_add_usage(std::size_t copy_index) {
        if (copy_index < m_best_local_mes.size() || m_want_single_improvement) {
            // definitely used
            m_configuration_usage_variables.emplace_back(
                std::in_place_type_t<bool>(), true);
            return;
        }
        m_configuration_usage_variables.push_back(m_model.NewBoolVar());
    }

    /**
     * Add variables for a configuration copy.
     */
    void p_add_copy_vars(std::size_t copy_index) {
        std::vector<VarOrVal> copy_vars;
        bool all_open = true;
        if (copy_index < m_best_local_mes.size()) {
            all_open = false;
            m_propagator.reset_to_zero();
            Vertex vmes = m_best_local_mes[copy_index];
            if (push_vertex(m_propagator, vmes) < 0) {
                throw std::logic_error("Invalid interaction in MES!");
            }
        }
        for (Var v = 0, nv = m_clauses->num_vars(); v < nv; ++v) {
            Lit p = lit::positive_lit(v);
            if (all_open || m_propagator.is_open(p)) {
                copy_vars.emplace_back(std::in_place_type_t<BoolVar>{},
                                       m_model.NewBoolVar());
            } else {
                bool val = m_propagator.is_true(p);
                copy_vars.emplace_back(std::in_place_type_t<bool>{}, val);
            }
        }
        m_configuration_variables.emplace_back(std::move(copy_vars));
    }

    VarOrVal p_value_of(std::size_t copy_index, Lit l) {
        Var v = lit::var(l);
        bool negated = lit::negative(l);
        return std::visit(
            overloaded{
                [&](BoolVar& b) -> VarOrVal { return negated ? b.Not() : b; },
                [&](bool& b) -> VarOrVal { return negated ? !b : b; }},
            m_configuration_variables[copy_index][v]);
    }

    /**
     * Add interaction coverage variables for the current copy.
     */
    void p_copy_add_coverage(std::size_t copy_index) {
        bool all_open = (copy_index >= m_best_local_mes.size());
        auto usage_var = m_configuration_usage_variables[copy_index];
        bool definitely_used = std::holds_alternative<bool>(usage_var) &&
                               std::get<bool>(usage_var);
        m_interaction_coverage_variables.emplace_back();
        auto& coverage_vars = m_interaction_coverage_variables[copy_index];
        for (std::size_t vi = 0, n = m_all_vertices.size(); vi < n; ++vi) {
            Vertex v = m_all_vertices[vi];
            if (all_open || (m_propagator.is_open(v.first) &&
                             m_propagator.is_open(v.second)))
            {
                BoolVar vcov = m_model.NewBoolVar();
                coverage_vars.push_back(vcov);
                BoolVar v1 = std::get<BoolVar>(p_value_of(copy_index, v.first));
                BoolVar v2 =
                    std::get<BoolVar>(p_value_of(copy_index, v.second));
                m_model.AddBoolOr({vcov.Not(), v1});
                m_model.AddBoolOr({vcov.Not(), v2});
                if (definitely_used || m_remaining_configuration) {
                    // if the class is definitely used, or if we have
                    // a remaining configuration covering none of
                    // m_all_vertices, we can just add the backwards clause
                    // ensuring vcov == v1 & v2
                    m_model.AddBoolOr({vcov, v1.Not(), v2.Not()});
                } else {
                    // otherwise, we ensure that vcov == u & v1 & v2
                    // where u is the variable indicating the class is used
                    BoolVar u = std::get<BoolVar>(usage_var);
                    m_model.AddBoolOr({vcov.Not(), u});
                    m_model.AddBoolOr({vcov, u.Not(), v1.Not(), v2.Not()});
                }
                continue;
            }
            // here, we are in a definitely used class!
            // otherwise, all_open would be true.
            if (m_propagator.is_false(v.first) ||
                m_propagator.is_false(v.second))
            {
                coverage_vars.push_back(false); // v is not covered here!
                continue;
            }
            if (m_propagator.is_true(v.first) && m_propagator.is_true(v.second))
            {
                coverage_vars.push_back(true); // v is covered here!
                continue;
            }
            // one true, one open: if the open variable is made true, v is
            // covered here!
            Lit vopen = m_propagator.is_open(v.first) ? v.first : v.second;
            coverage_vars.push_back(p_value_of(copy_index, vopen));
        }
    }

    /**
     * Add the relevant/remaining clauses for the given copy to the model.
     */
    void p_copy_add_clauses(std::size_t copy_index) {
        for (Lit u : m_clauses->unary_literals()) {
            const Lit ubuf[1] = {u};
            p_copy_add_clause(copy_index, ubuf, ubuf + 1);
        }
        for (auto b : m_clauses->binary_clauses()) {
            const Lit bbuf[2] = {b.first, b.second};
            p_copy_add_clause(copy_index, bbuf, bbuf + 2);
        }
        for (CRef c = 1, n = m_clauses->literal_db_size(); c < n;
             c = m_clauses->next_clause(c))
        {
            auto lits = m_clauses->lits_of(c);
            p_copy_add_clause(copy_index, lits.begin(), lits.end());
        }
    }

    /**
     * Add a clause attached to the given copy.
     * Only add literals that aren't set to false.
     * Only add clauses that aren't already satisfied.
     */
    template <typename LitIter>
    void p_copy_add_clause(std::size_t copy_index, LitIter lbegin,
                           LitIter lend) {
        bool is_already_satisfied = false;
        m_clause_buffer.clear();
        const auto& copy_vars = m_configuration_variables[copy_index];
        for (Lit l : IteratorRange{lbegin, lend}) {
            VarOrVal entry = copy_vars[lit::var(l)];
            if (std::holds_alternative<bool>(entry)) {
                if (std::get<bool>(entry) != lit::negative(l)) {
                    is_already_satisfied = true;
                    break;
                }
            } else {
                BoolVar bv = std::get<BoolVar>(entry);
                m_clause_buffer.push_back(lit::negative(l) ? bv.Not() : bv);
            }
        }
        if (!is_already_satisfied) {
            if (m_clause_buffer.empty()) {
                throw std::logic_error(
                    "Invalid MES vertex or unsatisfiable configuration model!");
            }
            m_model.AddBoolOr(m_clause_buffer);
        }
    }

    /**
     * Add constraints to make sure all interactions are covered.
     */
    void p_add_coverage_constraints() {
        std::size_t num_copies = m_covering_assignments.size() - 1;
        for (std::size_t vi : range(m_all_vertices.size())) {
            m_clause_buffer.clear();
            bool already_satisfied = false;
            for (std::size_t ci : range(num_copies)) {
                VarOrVal vcov = m_interaction_coverage_variables[ci][vi];
                if (std::holds_alternative<bool>(vcov)) {
                    if (std::get<bool>(vcov)) {
                        already_satisfied = true;
                        break;
                    }
                } else {
                    m_clause_buffer.push_back(std::get<BoolVar>(vcov));
                }
            }
            if (!already_satisfied) {
                if (m_clause_buffer.empty()) {
                    throw std::logic_error("Uncoverable interaction in model!");
                }
                m_model.AddBoolOr(m_clause_buffer);
            }
        }
    }

    /**
     * Add the objective function (if we are not doing single improvements).
     */
    void p_add_objective() {
        if (m_want_single_improvement)
            return;
        std::size_t num_copies = m_covering_assignments.size() - 1;
        std::optional<BoolVar> previous = std::nullopt;
        operations_research::sat::LinearExpr sum;
        for (std::size_t i = 0; i < num_copies; ++i) {
            VarOrVal usage = m_configuration_usage_variables[i];
            if (std::holds_alternative<bool>(usage)) {
                if (!std::get<bool>(usage)) {
                    throw std::logic_error("Configuration fixed to unused!");
                }
                sum += 1;
                continue;
            } else {
                BoolVar vu = std::get<BoolVar>(usage);
                if (previous) {
                    m_model.AddBoolOr({vu.Not(), *previous});
                }
                sum += vu;
                previous = vu;
                if (m_remaining_configuration) {
                    // if we have a remaining configuration, we can
                    // force unused classes to the remaining configuration
                    // (which does not cover anything in m_all_vertices)
                    p_compute_fixed_from_remaining(i);
                    m_model.AddBoolAnd(m_clause_buffer).OnlyEnforceIf(vu.Not());
                }
            }
        }
        m_model.Minimize(sum);
    }

    /**
     * Compute fixed setting for unused configurations
     * from the remaining configuration (whose covered interactions are
     * all not in m_all_vertices).
     */
    void p_compute_fixed_from_remaining(std::size_t copy_index) {
        m_clause_buffer.clear();
        const auto& copy_vars = m_configuration_variables[copy_index];
        const auto& fixed_conf = *m_remaining_configuration;
        for (Var v : range(Var(m_clauses->num_vars()))) {
            BoolVar xv = std::get<BoolVar>(copy_vars[v]);
            if (fixed_conf[v]) {
                m_clause_buffer.push_back(xv);
            } else {
                m_clause_buffer.push_back(xv.Not());
            }
        }
    }

    /**
     * Report events if we have a recorder.
     */
    template <typename... PrintArgs>
    void p_report_event(std::string event_name, OutputObject data,
                        PrintArgs&&... args) {
        if (m_local_recorder) {
            m_local_recorder->store_event(std::move(event_name),
                                          std::move(data),
                                          std::forward<PrintArgs>(args)...);
        }
    }

    /**
     * Get a boolean value in the given solution.
     */
    bool p_get_value(VarOrVal v, const SolverResponse& response) {
        return std::visit(
            overloaded{
                [&](BoolVar& b) -> bool {
                    return operations_research::sat::SolutionBooleanValue(
                        response, b);
                },
                [&](bool& b) -> bool { return b; }},
            v);
    }
};

} // namespace sammy

#endif
#endif // SAMMY_ORTOOLS_SUPPORTED
#endif
==> ./gurobi.h <==
#ifndef SAMMY_GUROBI_H_INCLUDED_
#define SAMMY_GUROBI_H_INCLUDED_

#include <gurobi_c++.h>
#include <mutex>

namespace sammy {

namespace detail {

inline void start_environment(GRBEnv& env, bool quiet) {
    static std::mutex start_environment_lock;
    std::unique_lock<std::mutex> lock{start_environment_lock};
    if (quiet) {
        env.set(GRB_IntParam_OutputFlag, 0);
        env.set(GRB_IntParam_LogToConsole, 0);
    }
    env.start();
}

} // namespace detail

inline GRBEnv& gurobi_environment(bool quiet) {
    try {
        thread_local bool initialized = false;
        thread_local GRBEnv environment{true};
        if (!initialized) {
            detail::start_environment(environment, quiet);
            initialized = true;
        }
        return environment;
    } catch (const GRBException& ex) {
        std::stringstream out;
        out << "Gurobi exception during initialization, likely "
               "license-related: "
            << ex.getMessage();
        throw std::runtime_error(out.str());
    }
}

} // namespace sammy

#endif
==> ./class_completer.h <==
#ifndef SAMMY_CLASS_COMPLETER_H_INCLUDED_
#define SAMMY_CLASS_COMPLETER_H_INCLUDED_

#include "literals.h"
#include "pair_infeasibility_map.h"

namespace sammy {

/**
 * A class that is used to complete 'color classes', i.e.,
 * turn partial configurations into complete valid configurations,
 * to finalize our partial assignments.
 * May fail to complete classes if the partial assignment turns out
 * to be invalid; in that case, we may have to cover vertices
 * that become uncovered during conflict resolution.
 */
template <typename EventHandler> class ClassCompleter {
  public:
    ClassCompleter(Lit n_concrete, Lit n_all, EventHandler* handler)
        : m_handler(handler), n_concrete(n_concrete), n_all(n_all),
          m_phases(n_all, false), m_given_literals_bits(2 * n_concrete, false),
          m_failed_literals_bits(2 * n_concrete, false), m_failed_qpos(0) {}

    /**
     * Try to turn the partial assignment represented by
     * the given propatator into a complete, valid configuration.
     * @return true if the assignment could be completed without
     *              changing any of the already-assigned literals.
     */
    bool complete_class(SharedDBPropagator& prop) {
        if (prop.get_trail().size() == n_all)
            return true;
        p_complete_init(prop);
        bool any_open = true;
        while (any_open) {
            if (p_push_from_failed_queue(prop)) {
                continue;
            }
            any_open = p_scan_all(prop);
        }
        for (Lit l : m_failed_concretes) {
            if (prop.is_false(l)) {
                return false;
            }
        }
        return true;
    }

    /**
     * Compute the complete memory usage, in bytes, that this class completer
     * uses.
     */
    std::size_t total_memory_usage() const noexcept {
        return m_phases.bytes_used() + m_given_literals_bits.bytes_used() +
               m_failed_literals_bits.bytes_used() +
               m_given_concretes.capacity() * sizeof(Lit) +
               m_failed_concretes.capacity() * sizeof(Lit);
    }

  private:
    bool p_push_from_failed_queue(SharedDBPropagator& prop) {
        while (m_failed_qpos < m_failed_concretes.size()) {
            Lit l = m_failed_concretes[m_failed_qpos++];
            if (prop.is_open(l)) {
                p_push(prop, l);
                return true;
            }
        }
        return false;
    }

    bool p_scan_all(SharedDBPropagator& prop) {
        bool any_open = false;
        for (Lit v = 0; v < n_all; ++v) {
            Lit l = m_phases[v] ? lit::positive_lit(v) : lit::negative_lit(v);
            if (prop.is_open(l)) {
                any_open = true;
                if (!p_push(prop, l))
                    break;
            }
        }
        return any_open;
    }

    bool p_push(SharedDBPropagator& prop, Lit l) {
        struct AssignmentHandler {
            void assignment_undone(Lit l) const {
                if (lit::var(l) < that->n_concrete &&
                    that->m_given_literals_bits[l])
                {
                    if (!that->m_failed_literals_bits[l]) {
                        that->m_failed_literals_bits[l] = true;
                        that->m_failed_concretes.push_back(l);
                    }
                }
            }
            void assignment_forced(Lit) const {}

            ClassCompleter* that;
        };
        if (prop.push_level(l))
            return true;
        m_phases[lit::var(l)] = lit::negative(l);
        AssignmentHandler handler{this};
        if (!prop.resolve_conflicts(handler)) {
            throw UNSATError();
        }
        m_failed_qpos = 0;
        return false;
    }

    void p_complete_init(SharedDBPropagator& prop) {
        m_given_concretes = prop.get_trail();
        auto new_end =
            std::remove_if(m_given_concretes.begin(), m_given_concretes.end(),
                           [&](Lit l) { return lit::var(l) >= n_concrete; });
        m_given_concretes.erase(new_end, m_given_concretes.end());
        p_prop_to_buffer(prop, m_given_literals_bits);
        m_phases.reset();
        m_failed_concretes.clear();
        m_failed_qpos = 0;
        m_failed_literals_bits.reset();
    }

    void p_prop_to_buffer(const SharedDBPropagator& prop, Bitset& buffer) {
        buffer.reset();
        for (Lit l : prop.get_trail()) {
            if (lit::var(l) < n_concrete) {
                buffer[l] = true;
            }
        }
    }

    EventHandler* m_handler;
    Lit n_concrete, n_all;
    Bitset m_phases;
    Bitset m_given_literals_bits;
    Bitset m_failed_literals_bits;
    std::vector<Lit> m_given_concretes;
    std::vector<Lit> m_failed_concretes;
    std::size_t m_failed_qpos;
};

} // namespace sammy

#endif
==> ./clause_db.h <==
#ifndef HS_CLAUSE_DB_H_INCLUDED_
#define HS_CLAUSE_DB_H_INCLUDED_

#include "clause_reduction.h"
#include "literals.h"
#include "range.h"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace sammy {

/**
 * Clause reference type.
 * An unsigned integer type indexing a clause in
 * a clause database.
 */
using CRef = std::uint32_t;

/**
 * A structure keeping a snapshot of
 * the number of clauses in a clause database.
 */
struct ClauseCounts {
    CRef unary_clause_end;
    CRef binary_clause_end;
    CRef long_clause_end;
};

/**
 * @brief A database of clauses in internal format
 * (i.e., unsigned 0-based even-odd-encoding of literals)
 * that can be shared by multiple SharedDBPropagator instances.
 */
class ClauseDB {
  private:
    // the large array containing all literals of all clauses
    // that are not unary and not binary
    std::vector<Lit> m_literals;

    // array of all literals in unary clauses
    std::vector<Lit> m_unaries;

    // array of arrays of binary watches (other literal)
    std::vector<std::vector<Lit>> m_binaries;

    // array of binary clause literals:
    // duplicated to make keeping track of learnt binaries possible,
    // but separated from regular clauses because accesses
    // should be much rarer
    std::vector<std::pair<Lit, Lit>> m_binary_clauses;

    // number of variables
    std::uint32_t m_num_vars, m_num_clauses;

    // marked as frozen?
    // frozen: must not add new clauses.
    bool m_frozen = false;

    // turn an external clause into an internal one
    void p_internalize_external(const ExternalClause& c) {
        m_literals.push_back(c.size());
        for (auto l : c) {
            m_literals.push_back(lit::internalize(l));
        }
    }

    // turn an array of external clauses (that should be reduced)
    // into internal clauses
    void
    p_internalize_reduced_externals(const std::vector<ExternalClause>& cs) {
        m_num_clauses = cs.size();
        for (const auto& c : cs) {
            if (c.size() == 1) {
                m_unaries.push_back(lit::internalize(c[0]));
            } else if (c.size() == 2) {
                std::uint32_t l1 = lit::internalize(c[0]);
                std::uint32_t l2 = lit::internalize(c[1]);
                m_binaries[l1].push_back(l2);
                m_binaries[l2].push_back(l1);
                m_binary_clauses.emplace_back(l1, l2);
            } else {
                p_internalize_external(c);
            }
        }
    }

    // turn an external clause that we known nothing about
    // (might be tautological/contain duplicates) into an internal one
    void p_internalize_internal(const std::vector<std::vector<Lit>>& cs) {
        m_num_clauses = cs.size();
        std::vector<Lit> buffer;
        for (const auto& c : cs) {
            bool tautology = (std::find_if(c.begin(), c.end(), [&](Lit l) {
                                  return std::find(c.begin(), c.end(),
                                                   lit::negate(l)) != c.end();
                              }) != c.end());
            if (tautology) {
                --m_num_clauses;
                continue;
            }
            buffer = c;
            std::sort(buffer.begin(), buffer.end());
            buffer.erase(std::unique(buffer.begin(), buffer.end()),
                         buffer.end());
            if (c.size() == 1) {
                m_unaries.push_back(buffer[0]);
            } else if (c.size() == 2) {
                std::uint32_t l1 = buffer[0];
                std::uint32_t l2 = buffer[1];
                m_binaries[l1].push_back(l2);
                m_binaries[l2].push_back(l1);
                m_binary_clauses.emplace_back(l1, l2);
            } else {
                m_literals.push_back(c.size());
                for (auto l : c) {
                    m_literals.push_back(l);
                }
            }
        }
    }

  public:
    using LitIterator = const Lit*;
    using Lits = IteratorRange<LitIterator>;

    ClauseCounts get_clause_counts() const noexcept {
        return {num_unaries(), num_binaries(), literal_db_size()};
    }

    static ClauseDB import_from(std::istream& input) {
        std::uint32_t num_vars = 0;
        input >> num_vars;
        std::vector<ExternalClause> clauses;
        while (true) {
            ExternalClause cl;
            std::uint32_t s = 0;
            std::int32_t l;
            if (!(input >> s))
                break;
            cl.reserve(s);
            for (std::uint32_t i = 0; i < s; ++i) {
                input >> l;
                cl.push_back(l);
            }
            clauses.emplace_back(std::move(cl));
        }
        return ClauseDB(num_vars, clauses);
    }

    static ClauseDB import_from(const std::filesystem::path& path) {
        std::ifstream input;
        input.exceptions(std::ios::failbit | std::ios::badbit);
        input.open(path, std::ios::in);
        input.exceptions(std::ios::badbit);
        return import_from(input);
    }

    explicit ClauseDB(std::uint32_t num_vars)
        : m_binaries(2 * num_vars, std::vector<Lit>{}), m_num_vars(num_vars),
          m_num_clauses(0) {}

    explicit ClauseDB(std::uint32_t num_vars,
                      const std::vector<ExternalClause>& clauses)
        : m_binaries(2 * num_vars, std::vector<Lit>{}), m_num_vars(num_vars),
          m_num_clauses(0) {
        auto reduced = reduce_external_clauses(clauses);
        p_internalize_reduced_externals(reduced);
    }

    explicit ClauseDB(std::uint32_t num_vars,
                      const std::vector<std::vector<Lit>>& clauses)
        : m_binaries(2 * num_vars, std::vector<Lit>{}), m_num_vars(num_vars),
          m_num_clauses(0) {
        p_internalize_internal(clauses);
    }

    void export_to(std::ostream& output) const {
        output << m_num_vars << '\n';
        auto all_clauses = export_all_clauses();
        for (const auto& ec : all_clauses) {
            output << ec.size();
            for (auto l : ec) {
                output << ' ' << l;
            }
            output << '\n';
        }
    }

    /**
     * Mark the clause database as frozen.
     */
    void mark_frozen() noexcept { m_frozen = true; }

    /**
     * Check if the clause database is marked as frozen.
     */
    bool is_frozen() const noexcept { return m_frozen; }

    void export_to(const std::filesystem::path& path) const {
        std::ofstream output(path, std::ios::out | std::ios::trunc);
        output.exceptions(std::ios::badbit | std::ios::failbit);
        export_to(output);
    }

    void export_to_dimacs(std::ostream& output) const {
        output << "p cnf " << m_num_vars << ' ' << m_num_clauses << std::endl;
        auto all_clauses = export_all_clauses();
        for (const auto& ec : all_clauses) {
            for (auto l : ec) {
                output << l << ' ';
            }
            output << '0' << '\n';
        }
    }

    void export_to_dimacs(const std::filesystem::path& path) const {
        std::ofstream output(path, std::ios::out | std::ios::trunc);
        output.exceptions(std::ios::badbit | std::ios::failbit);
        export_to_dimacs(output);
    }

    std::vector<ExternalClause> export_all_clauses() const {
        std::vector<ExternalClause> result;
        result.reserve(m_num_clauses);
        for (Lit l : m_unaries) {
            result.emplace_back(
                std::initializer_list<ExternalLit>{lit::externalize(l)});
        }
        for (auto b : m_binary_clauses) {
            result.emplace_back(std::initializer_list<ExternalLit>{
                lit::externalize(b.first), lit::externalize(b.second)});
        }
        for (CRef current = 1, s = literal_db_size(); current < s;
             current = next_clause(current))
        {
            auto ls = lits_of(current);
            result.emplace_back();
            auto& out = result.back();
            out.reserve(ls.end() - ls.begin());
            for (Lit l : ls) {
                out.push_back(lit::externalize(l));
            }
        }
        return result;
    }

    template <typename Iterator> CRef add_clause(Iterator ibeg, Iterator iend) {
        if (m_frozen) {
            throw std::logic_error("Cannot add clauses to frozen ClauseDB!");
        }
        CRef res = NIL;
        auto length = std::distance(ibeg, iend);
        if (length == 1) {
            m_unaries.push_back(*ibeg);
        } else if (length == 2) {
            Lit l1 = ibeg[0];
            Lit l2 = ibeg[1];
            m_binaries[l1].push_back(l2);
            m_binaries[l2].push_back(l1);
            m_binary_clauses.emplace_back(l1, l2);
        } else {
            m_literals.push_back(length);
            res = m_literals.size();
            m_literals.insert(m_literals.end(), ibeg, iend);
        }
        ++m_num_clauses;
        return res;
    }

    CRef add_clause(Lits literals) {
        return add_clause(literals.begin(), literals.end());
    }

    CRef cref_of(Lits lits) const noexcept {
        return CRef(lits.begin() - m_literals.data());
    }

    Lits lits_of(CRef clause) const noexcept {
        LitIterator begin = m_literals.data() + clause;
        return {begin, begin + begin[-1]};
    }

    std::uint32_t clause_length(CRef clause) const noexcept {
        return m_literals[clause - 1];
    }

    CRef next_clause(CRef clause) const noexcept {
        return clause + m_literals[clause - 1] + 1;
    }

    std::uint32_t num_vars() const noexcept { return m_num_vars; }

    std::uint32_t num_clauses() const noexcept { return m_num_clauses; }

    std::uint32_t num_unaries() const noexcept { return m_unaries.size(); }

    std::uint32_t num_binaries() const noexcept {
        return m_binary_clauses.size();
    }

    std::size_t total_clause_size() const noexcept {
        return num_unaries() + 2 * num_binaries() + literal_db_size();
    }

    Lits unary_literals(std::uint32_t begin, std::uint32_t end) const noexcept {
        return {m_unaries.data() + begin, m_unaries.data() + end};
    }

    Lits unary_literals() const noexcept {
        return unary_literals(0, m_unaries.size());
    }

    using BinaryClauseIterator = const std::pair<Lit, Lit>*;
    using BinaryClauses = IteratorRange<BinaryClauseIterator>;
    BinaryClauses binary_clauses(std::uint32_t begin,
                                 std::uint32_t end) const noexcept {
        return {m_binary_clauses.data() + begin, m_binary_clauses.data() + end};
    }
    BinaryClauses binary_clauses() const noexcept {
        return binary_clauses(0, m_binary_clauses.size());
    }

    const std::vector<Lit>& binary_partners_of(Lit lit) const noexcept {
        return m_binaries[lit];
    }

    std::uint32_t literal_db_size() const noexcept { return m_literals.size(); }

    std::size_t total_memory_usage() const noexcept {
        std::size_t total_binary_usage = 0;
        for (const auto& v : m_binaries) {
            total_binary_usage += v.capacity();
        }
        total_binary_usage *= sizeof(Lit);
        return m_literals.capacity() * sizeof(Lit) +
               m_unaries.capacity() * sizeof(Lit) +
               m_binaries.capacity() * sizeof(std::vector<Lit>) +
               total_binary_usage +
               m_binary_clauses.capacity() * sizeof(std::pair<Lit, Lit>) +
               sizeof(ClauseDB);
    }
};

class ClauseDBView {
    ClauseDB* db;
    std::vector<std::uint32_t> seen_binary_partners_of;
    std::uint32_t seen_unary_count;
    std::uint32_t seen_binary_count;
    CRef seen_literal_end;

    template <typename NewClauseHandler>
    bool p_handle_new_unaries(NewClauseHandler& handler) {
        if (seen_unary_count != db->num_unaries()) {
            for (Lit l :
                 db->unary_literals(seen_unary_count, db->num_unaries()))
            {
                ++seen_unary_count;
                if (!handler.new_unary(l)) {
                    return false;
                }
            }
        }
        return true;
    }

    template <typename NewClauseHandler>
    bool p_handle_new_larger_clauses(NewClauseHandler& handler) {
        while (seen_literal_end != db->literal_db_size()) {
            CRef clause = seen_literal_end + 1;
            seen_literal_end = db->next_clause(clause) - 1;
            if (!handler.new_clause(db->lits_of(clause))) {
                return false;
            }
        }
        return true;
    }

    template <typename NewClauseHandler>
    bool p_handle_new_binaries(NewClauseHandler& handler) {
        if (seen_binary_count != db->num_binaries()) {
            for (const auto& bc :
                 db->binary_clauses(seen_binary_count, db->num_binaries()))
            {
                ++seen_binary_count;
                seen_binary_partners_of[bc.first]++;
                seen_binary_partners_of[bc.second]++;
                if (!handler.new_binary(bc.first, bc.second)) {
                    return false;
                }
            }
        }
        return true;
    }

  public:
    explicit ClauseDBView(ClauseDB* db) noexcept
        : db(db), seen_binary_partners_of(2 * db->num_vars(), 0),
          seen_unary_count(0), seen_binary_count(0), seen_literal_end(0) {}

    template <typename NewClauseHandler>
    void handle_new_unaries(NewClauseHandler& handler) {
        p_handle_new_unaries(handler);
    }

    template <typename NewClauseHandler>
    void handle_new_long_clauses(NewClauseHandler& handler) {
        p_handle_new_larger_clauses(handler);
    }

    template <typename NewClauseHandler>
    void handle_new_binary_clauses(NewClauseHandler& handler) {
        p_handle_new_binaries(handler);
    }

    template <typename NewClauseHandler>
    void handle_new_clauses(NewClauseHandler& handler) {
        if (!p_handle_new_unaries(handler))
            return;
        if (!p_handle_new_binaries(handler))
            return;
        if (!p_handle_new_larger_clauses(handler))
            return;
    }

    ClauseDB& database() noexcept { return *db; }

    const ClauseDB& database() const noexcept { return *db; }

    ClauseDB::Lits binary_partners_of(Lit l) const noexcept {
        std::size_t lim = seen_binary_partners_of[l];
        const auto& partners = db->binary_partners_of(l);
        return {partners.data(), partners.data() + lim};
    }

    std::size_t total_memory_usage() const noexcept {
        return 4 * seen_binary_partners_of.capacity() + sizeof(ClauseDBView);
    }
};

/**
 * Turn the given ClauseDB into a list of clauses (internal representation,
 * i.e., 0-based indexing and even/odd encoding of true/false literals).
 */
template <typename ClauseType>
inline std::vector<ClauseType> to_clause_list(const ClauseDB& clauses) {
    std::vector<ClauseType> result;
    for (Lit l : clauses.unary_literals()) {
        Lit a[1] = {l};
        result.emplace_back(+a, a + 1);
    }
    for (auto [l1, l2] : clauses.binary_clauses()) {
        Lit a[2] = {l1, l2};
        result.emplace_back(+a, a + 2);
    }
    for (CRef c = 1, ndb = clauses.literal_db_size(); c < ndb;
         c = clauses.next_clause(c))
    {
        auto lits = clauses.lits_of(c);
        result.emplace_back(lits.begin(), lits.end());
    }
    return result;
}

} // namespace sammy

#endif
==> ./bounded_variable_elimination.h <==
#ifndef SAMMY_BOUNDED_VARIABLE_ELIMINATION_H_INCLUDED_
#define SAMMY_BOUNDED_VARIABLE_ELIMINATION_H_INCLUDED_

#include "literals.h"
#include "simplify_datastructure.h"

namespace sammy {

/**
 * Can we use BVE (bounded variable elimination) on non-concrete features?
 *  - Let (l1,l2) be a feasible interaction.
 *  - Assume we eliminate some non-concrete y by VE.
 *  - This adds the resolvents of all clauses with y and all clauses with -y.
 *  - As such, these clauses preserve logical equivalence: any solution to
 *    the old formula is a solution (including y) to the new formula.
 *  - Dropping y transforms this into a solution of the new formula, since
 *    removal of clauses maintains the solution in that direction.
 *  - Let (p1, p2) be feasible after the transformation.
 *  - Let Q be the model including (p1, p2) on the new formula.
 *  - Let Q+ be Q extended by y = True.
 *  - Analogously, let Q- be the model Q extended by y = False.
 *  - If Q+ is not a model of the original formula, there is a clause
 *    that is not satisfied; this must be one of the dropped clauses containing
 * -y.
 *  - Thus, all resolvents of this clause with clauses containing y must be
 *    satisfied in Q by satisfying the parts coming from clauses containing y.
 *  - Therefore, Q- must be a model of the original formula.
 * We should be allowed to use BVE on non-concrete features!
 * Do we introduce another dummy? Or some other way to represent elimination?
 */

class BoundedVariableEliminator {
  public:
    explicit BoundedVariableEliminator(SimplifyDatastructure* simplify_ds)
        : simplify_ds(simplify_ds),
          clauses_with_literal(2 * simplify_ds->original_num_vars()),
          stats(&stats_buffer) {
        p_init_cl_with_lit();
    }

    void set_stats(SimplificationStats* stats) noexcept { this->stats = stats; }

    bool is_clause_eliminated(CRef c) {
        std::size_t nc = simplify_ds->clauses().size();
        std::size_t nec = clause_eliminated.size();
        if (nec < nc) {
            clause_eliminated.reserve(std::max(2 * nec, nc + 32));
            clause_eliminated.resize(nc, false);
        }
        return clause_eliminated[c];
    }

    void set_clause_eliminated(CRef c) {
        std::size_t nc = simplify_ds->clauses().size();
        std::size_t nec = clause_eliminated.size();
        if (nec < nc) {
            clause_eliminated.reserve(std::max(2 * nec, nc + 32));
            clause_eliminated.resize(nc, false);
        }
        clause_eliminated[c] = true;
    }

    bool run_elimination(std::size_t max_gap) {
        bool changed = false, once_changed = false;
        do {
            changed = false;
            std::vector<std::pair<Var, std::size_t>> best_scores;
            changed = p_score_and_eliminate_pure(best_scores);
            if (!changed) {
                p_compress_scored(best_scores);
                bool eliminated_pure = false;
                for (const auto& e : best_scores) {
                    if (e.second == 0) {
                        changed = eliminated_pure = true;
                        p_eliminate_pure(e.first);
                        continue;
                    }
                    if (eliminated_pure)
                        break;
                    auto [score, baseline] = p_actual_score(e.first);
                    if (score <= baseline + max_gap) {
                        p_eliminate_variable(e.first);
                        changed = true;
                        break;
                    }
                }
            }
            once_changed |= changed;
        } while (changed);
        return once_changed;
    }

  private:
    bool p_score_and_eliminate_pure(
        std::vector<std::pair<Var, std::size_t>>& best_scores) {
        bool changed = false;
        const Var nv = simplify_ds->original_num_vars();
        for (Var v = 0; v < nv; ++v) {
            if (simplify_ds->is_eliminated(v) || simplify_ds->is_concrete(v))
                continue;
            Lit p = lit::positive_lit(v), n = lit::negative_lit(v);
            const auto& plist = clauses_with_literal[p];
            const auto& nlist = clauses_with_literal[n];
            if (plist.empty() || nlist.empty()) {
                p_eliminate_pure(v);
                changed = true;
                continue;
            }
            if (!changed)
                p_add_score(best_scores, v, plist.size() * nlist.size(), 20);
        }
        return changed;
    }

    void p_eliminate_variable(Var v) {
        Lit p = lit::positive_lit(v), n = lit::negative_lit(v);
        const auto& plist = clauses_with_literal[p];
        const auto& nlist = clauses_with_literal[n];
        auto& clauses = simplify_ds->clauses();
        for (CRef c1 : plist) {
            p_push_reconstruct(clauses[c1], p);
            set_clause_eliminated(c1);
        }
        for (CRef c2 : nlist) {
            p_push_reconstruct(clauses[c2], n);
            set_clause_eliminated(c2);
        }
        auto not_p = [=](Lit l) { return l != p; };
        auto not_n = [=](Lit l) { return l != n; };
        for (CRef c1 : plist) {
            for (CRef c2 : nlist) {
                SCVec buffer;
                std::copy_if(clauses[c1].begin(), clauses[c1].end(),
                             std::back_inserter(buffer), not_p);
                std::copy_if(clauses[c2].begin(), clauses[c2].end(),
                             std::back_inserter(buffer), not_n);
                if (buffer.empty())
                    throw UNSATError();
                std::sort(buffer.begin(), buffer.end());
                buffer.erase(std::unique(buffer.begin(), buffer.end()),
                             buffer.end());
                if (find_pair_if(buffer.begin(), buffer.end(),
                                 [](Lit pr, Lit c) {
                                     return pr == lit::negate(c);
                                 }) != buffer.end())
                {
                    continue;
                }
                CRef nclause = clauses.size();
                for (Lit l : buffer) {
                    clauses_with_literal[l].push_back(nclause);
                }
                clauses.emplace_back(std::move(buffer));
            }
        }
        simplify_ds->mark_eliminated(v);
        stats->variables_eliminated_by_resolution += 1;
    }

    std::pair<std::size_t, std::size_t> p_actual_score(Var v) {
        Lit p = lit::positive_lit(v), n = lit::negative_lit(v);
        std::size_t actual_score = 0;
        const auto& plist = clauses_with_literal[p];
        const auto& nlist = clauses_with_literal[n];
        const auto& clauses = simplify_ds->clauses();
        for (CRef c1 : plist) {
            for (CRef c2 : nlist) {
                actual_score +=
                    (p_resolvent_is_tautology(clauses[c1], clauses[c2], v) ? 0
                                                                           : 1);
            }
        }
        return {actual_score, plist.size() + nlist.size()};
    }

    bool p_resolvent_is_tautology(const SCVec& cpos, const SCVec& cneg,
                                  Var v) const {
        Lit p = lit::positive_lit(v), n = lit::negative_lit(v);
        for (Lit lp : cpos) {
            if (lp == p)
                continue;
            for (Lit ln : cneg) {
                if (ln == n)
                    continue;
                if (lp == lit::negate(ln)) {
                    return true;
                }
            }
        }
        return false;
    }

    void p_compress_scored(std::vector<std::pair<Var, std::size_t>>& scores) {
        for (auto& e : scores) {
            e.second = p_compress_lists(e.first);
        }
        auto compare = [](const std::pair<Var, std::size_t>& v1,
                          const std::pair<Var, std::size_t>& v2) {
            return v1.second < v2.second;
        };
        std::sort(scores.begin(), scores.end(), compare);
    }

    void p_add_score(std::vector<std::pair<Var, std::size_t>>& scores, Var v,
                     std::size_t score, std::size_t size_cap) {
        if (scores.size() < size_cap) {
            scores.emplace_back(v, score);
            return;
        }
        if (scores.back().second <= score)
            return;
        std::pair<Var, std::size_t> val{v, score};
        auto compare = [](const std::pair<Var, std::size_t>& v1,
                          const std::pair<Var, std::size_t>& v2) {
            return v1.second < v2.second;
        };
        auto insert_pos =
            std::lower_bound(scores.begin(), scores.end(), val, compare);
        std::move_backward(insert_pos, scores.end() - 1, scores.end());
        *insert_pos = std::pair<Var, std::size_t>{v, score};
    }

    void p_push_reconstruct(const SCVec& clause, Lit nelit) {
        SCVec buffer(clause);
        std::swap(*std::find(buffer.begin(), buffer.end(), nelit),
                  buffer.front());
        simplify_ds->push_reconstruction_clause(std::move(buffer));
    }

    void p_eliminate_pure(Var v) {
        Lit p = lit::positive_lit(v), n = lit::negative_lit(v);
        const auto& plist = clauses_with_literal[p];
        const auto& nlist = clauses_with_literal[n];
        const auto& nelist = (plist.empty() ? nlist : plist);
        Lit nelit = (plist.empty() ? n : p);
        simplify_ds->push_reconstruction_clause(SCVec(1, nelit));
        simplify_ds->mark_eliminated(v);
        stats->pure_literals_eliminated += 1;
        for (CRef c : nelist) {
            set_clause_eliminated(c);
        }
    }

    void p_init_cl_with_lit() {
        const auto& clauses = simplify_ds->clauses();
        for (CRef ci = 0, cn = clauses.size(); ci < cn; ++ci) {
            const auto& clause = clauses[ci];
            for (Lit l : clause) {
                clauses_with_literal[l].push_back(ci);
            }
        }
    }

    std::size_t p_compress_lists(Var v) {
        Lit p = lit::positive_lit(v);
        Lit n = lit::negative_lit(v);
        p_compress_list(p);
        p_compress_list(n);
        return clauses_with_literal[p].size() * clauses_with_literal[n].size();
    }

    void p_compress_list(Lit l) {
        auto is_eliminated = [&](CRef c) { return is_clause_eliminated(c); };
        auto& list = clauses_with_literal[l];
        list.erase(std::remove_if(list.begin(), list.end(), is_eliminated),
                   list.end());
    }

    SimplifyDatastructure* simplify_ds;
    std::vector<bool> clause_eliminated;
    std::vector<std::vector<CRef>> clauses_with_literal;
    SimplificationStats stats_buffer;
    SimplificationStats* stats;
};

inline bool bounded_variable_elimination(SimplifyDatastructure& simplifier,
                                         std::size_t max_gain_clauses,
                                         SimplificationStats* stats = nullptr) {
    BoundedVariableEliminator eliminator{&simplifier};
    if (stats)
        eliminator.set_stats(stats);
    bool result = eliminator.run_elimination(max_gain_clauses);
    if (result)
        simplifier.remove_eliminated_clauses();
    return result;
}

} // namespace sammy

#endif
==> ./thread_group.h <==
#ifndef SAMMY_THREAD_GROUP_H_INCLUDED_
#define SAMMY_THREAD_GROUP_H_INCLUDED_

#include <cassert>
#include <condition_variable>
#include <future>
#include <list>
#include <mutex>
#include <thread>

namespace sammy {

/**
 * A group/pool of threads that handle tasks
 * which all return results of a specific type R.
 */
template <class R> class ThreadGroup {
  public:
    using ThreadIndex = std::size_t;
    using FunctionType = R(ThreadIndex);
    using TaskType = std::packaged_task<FunctionType>;
    using ResultType = R;

    explicit ThreadGroup(std::size_t desired_extra_threads) noexcept
        : m_thread_count(desired_extra_threads) {}

    explicit ThreadGroup() noexcept
        : ThreadGroup(std::thread::hardware_concurrency() - 1) {}

    // no copy or move allowed
    ThreadGroup(const ThreadGroup&) = delete;
    ThreadGroup(ThreadGroup&&) = delete;
    ThreadGroup& operator=(const ThreadGroup&) = delete;
    ThreadGroup& operator=(ThreadGroup&&) = delete;

    ~ThreadGroup() { stop(); }

    void stop() noexcept {
        {
            LockType l{m_mutex};
            m_stopped = true;
            m_condvar.notify_all();
        }
        std::for_each(m_threads.begin(), m_threads.end(),
                      [](std::thread& t) { t.join(); });
        m_threads.clear();
        m_stopped = false;
    }

    std::future<ResultType> post(TaskType&& t) {
        assert(t.valid());
        p_ensure_running();
        std::future<ResultType> res = t.get_future();
        {
            LockType l{m_mutex};
            m_tasks.emplace_back(std::move(t));
            m_condvar.notify_one();
        }
        return res;
    }

    template <typename TaskInputIterator, typename FutureOutputIterator>
    void post(TaskInputIterator task_begin, TaskInputIterator task_end,
              FutureOutputIterator out) {
        LockType l{m_mutex};
        p_ensure_running();
        for (; task_begin != task_end; ++task_begin) {
            std::future<ResultType> res = task_begin->get_future();
            m_tasks.emplace_back(std::move(*task_begin));
            *out = std::move(res);
            ++out;
        }
        m_condvar.notify_all();
    }

    template <typename RandomAccessIterator, typename Callable>
    void parallel_foreach(RandomAccessIterator begin, RandomAccessIterator end,
                          Callable&& callable) {
        parallel_foreach_iterator(begin, end,
                                  [&](const RandomAccessIterator& iter) {
                                      std::forward<Callable>(callable)(*iter);
                                  });
    }

    template <typename RandomAccessIterator, typename Callable>
    void parallel_foreach_iterator(RandomAccessIterator begin,
                                   RandomAccessIterator end,
                                   Callable&& callable) {
        if (begin == end)
            return;
        RandomAccessIterator last_beg = begin;
        std::vector<std::future<ResultType>>& futures = p_make_futures();
        {
            LockType l{m_mutex};
            p_ensure_running();
            split_range(
                begin, end, num_threads() + 1,
                [&](RandomAccessIterator s_beg, RandomAccessIterator s_end) {
                    if (s_beg == end)
                        return;
                    if (s_end == end) {
                        last_beg = s_beg;
                    } else {
                        m_tasks.emplace_back([s_beg, s_end,
                                              &callable](ThreadIndex) {
                            for (auto s_cur = s_beg; s_cur != s_end; ++s_cur) {
                                callable(s_cur);
                            }
                        });
                        futures.emplace_back(
                            std::move(m_tasks.back().get_future()));
                    }
                });
            m_condvar.notify_all();
        }
        for (RandomAccessIterator iter = last_beg; iter != end; ++iter) {
            callable(iter);
        }
        std::for_each(futures.begin(), futures.end(),
                      [](std::future<ResultType>& f) { f.get(); });
        futures.clear();
    }

    template <typename RandomAccessIterator, typename ContextType,
              typename Callable>
    void parallel_foreach(RandomAccessIterator begin, RandomAccessIterator end,
                          const ContextType& context, Callable&& callable) {
        parallel_foreach_iterator(
            begin, end, context,
            [&](ContextType& local_ctx, RandomAccessIterator iter) {
                std::forward<Callable>(callable)(local_ctx, *iter);
            });
    }

    template <typename RandomAccessIterator, typename ContextType,
              typename Callable>
    void parallel_foreach_iterator(RandomAccessIterator begin,
                                   RandomAccessIterator end,
                                   const ContextType& context,
                                   Callable&& callable) {
        if (begin == end)
            return;
        std::vector<std::future<ResultType>>& futures = p_make_futures();
        RandomAccessIterator last_beg = begin;
        {
            LockType l{m_mutex};
            p_ensure_running();
            split_range(
                begin, end, num_threads() + 1,
                [&](RandomAccessIterator s_beg, RandomAccessIterator s_end) {
                    if (s_beg == end)
                        return;
                    if (s_end == end) {
                        last_beg = s_beg;
                    } else {
                        m_tasks.emplace_back([s_beg, s_end, &callable,
                                              &context](ThreadIndex) {
                            ContextType local_context{context};
                            for (auto s_cur = s_beg; s_cur != s_end; ++s_cur) {
                                callable(local_context, s_cur);
                            }
                        });
                        futures.emplace_back(
                            std::move(m_tasks.back().get_future()));
                    }
                });
            m_condvar.notify_all();
        }
        ContextType local_context{context};
        for (RandomAccessIterator iter = last_beg; iter != end; ++iter) {
            callable(local_context, iter);
        }
        std::for_each(futures.begin(), futures.end(),
                      [](std::future<ResultType>& f) { f.get(); });
        futures.clear();
    }

    template <typename RandomAccessIterator,
              typename ContextFunction /*(ThreadIndex,RandomAccessIterator)*/,
              typename Callable, /*(Context&, ThreadIndex,
                                    RandomAccessIterator)*/
              typename FinalContextFunction /*(Context&, ThreadIndex)*/>
    void context_function_parallel_foreach_iterator(
        RandomAccessIterator begin, RandomAccessIterator end,
        ContextFunction&& context_fn, Callable&& callable,
        FinalContextFunction&& final_ctx_fn) {
        if (begin == end)
            return;
        std::vector<std::future<ResultType>>& futures = p_make_futures();
        RandomAccessIterator last_beg = begin;
        {
            LockType l{m_mutex};
            p_ensure_running();
            split_range(
                begin, end, num_threads() + 1,
                [&](RandomAccessIterator s_beg, RandomAccessIterator s_end) {
                    if (s_beg == end)
                        return;
                    if (s_end == end) {
                        last_beg = s_beg;
                    } else {
                        m_tasks.emplace_back(
                            [s_beg, s_end, &context_fn, &callable,
                             &final_ctx_fn](ThreadIndex index) {
                                auto s_cur = s_beg;
                                auto& local_context = context_fn(index, s_cur);
                                for (++s_cur; s_cur != s_end; ++s_cur) {
                                    callable(local_context, index, s_cur);
                                }
                                final_ctx_fn(local_context, index);
                            });
                        futures.emplace_back(
                            std::move(m_tasks.back().get_future()));
                    }
                });
            m_condvar.notify_all();
        }
        auto& local_context = context_fn(0, last_beg);
        for (++last_beg; last_beg != end; ++last_beg) {
            callable(local_context, 0, last_beg);
        }
        final_ctx_fn(local_context, 0);
        std::for_each(futures.begin(), futures.end(),
                      [](std::future<ResultType>& f) { f.get(); });
        futures.clear();
    }

    template <typename Callable>
    void run_n_copies(std::size_t n, Callable&& callable) {
        std::vector<std::future<ResultType>>& futures = p_make_futures();
        futures.reserve(n - 1);
        {
            LockType l{m_mutex};
            p_ensure_running();
            for (std::size_t i = 0; i < n - 1; ++i) {
                m_tasks.emplace_back([&callable](ThreadIndex) {
                    std::forward<Callable>(callable)();
                });
                futures.emplace_back(std::move(m_tasks.back().get_future()));
            }
            m_condvar.notify_all();
        }
        std::forward<Callable>(callable)();
        std::for_each(futures.begin(), futures.end(),
                      [](std::future<ResultType>& f) { f.get(); });
        futures.clear();
    }

    std::size_t num_threads() const noexcept { return m_thread_count; }

    std::size_t num_threads_running() const noexcept {
        return m_threads.size();
    }

  private:
    using LockType = std::unique_lock<std::mutex>;

    std::vector<std::future<ResultType>>& p_make_futures() {
        m_futures_cache.reserve(num_threads());
        return m_futures_cache;
    }

    void p_worker_main() {
        const ThreadIndex my_index = m_thread_index_counter++;
        for (;;) {
            TaskType task;
            {
                LockType l{m_mutex};
                if (!p_await_task(l, task))
                    return;
            }
            task(my_index);
        }
    }

    void p_ensure_running() {
        while (m_threads.size() < m_thread_count) {
            m_threads.emplace_back([this]() { p_worker_main(); });
        }
    }

    bool p_should_wake() const {
        if (m_stopped)
            return true;
        return !m_tasks.empty();
    }

    bool p_await_task(LockType& held_lock, TaskType& out) {
        m_condvar.wait(held_lock, [this]() -> bool { return p_should_wake(); });
        if (m_tasks.empty())
            return false;
        out = std::move(m_tasks.back());
        m_tasks.pop_back();
        return true;
    }

    std::vector<std::thread> m_threads;
    std::size_t m_thread_count;
    std::mutex m_mutex;
    std::condition_variable m_condvar;
    std::vector<TaskType> m_tasks;
    std::vector<std::future<ResultType>> m_futures_cache;
    bool m_stopped = false;
    std::atomic<ThreadIndex> m_thread_index_counter{1};
};

/**
 * A buffer used by ResultToLocalContext to store
 * each thread's local context objects in.
 */
template <typename ResultType> class ResultToLocalContextBuffer {
  public:
    template <typename... ResultConsArgs>
    explicit ResultToLocalContextBuffer(ResultConsArgs&&... args) {
        m_results.emplace_back(std::forward<ResultConsArgs>(args)...);
    }

    ResultType& front() noexcept { return m_results.front(); }
    const ResultType& front() const noexcept { return m_results.front(); }

    ResultType* create_pointer() {
        std::unique_lock<std::mutex> l{m_mutex};
        m_results.emplace_back(m_results.front());
        return &m_results.back();
    }

    template <typename ResultInitType, typename ReduceFn>
    ResultType reduce(ResultInitType&& init, ReduceFn&& reduce) {
        std::unique_lock<std::mutex> l{m_mutex};
        return std::reduce(std::next(m_results.begin()), m_results.end(),
                           std::forward<ResultInitType>(init),
                           std::forward<ReduceFn>(reduce));
    }

    template <typename ForEachFn> void for_each(ForEachFn&& callback) {
        std::for_each(std::next(m_results.begin()), m_results.end(),
                      std::forward<ForEachFn>(callback));
    }

  private:
    std::mutex m_mutex;
    std::list<ResultType> m_results;
};

/**
 * A local context structure that creates a new copy
 * of a given object in a std::list for each copy
 * made of the local context structure itself.
 * This is useful to allow a void ThreadGroup to
 * return values that need to be joined/reduced to
 * a single result later on.
 */
template <typename ResultType> class ResultToLocalContext {
  public:
    explicit ResultToLocalContext(
        ResultToLocalContextBuffer<ResultType>& buffer) noexcept
        : m_buffer(&buffer), m_local_pointer(&buffer.front()) {}

    explicit ResultToLocalContext(
        const ResultToLocalContext<ResultType>& other) noexcept
        : m_buffer(other.m_buffer),
          m_local_pointer(m_buffer->create_pointer()) {}

    ResultType& get() noexcept { return *m_local_pointer; }
    const ResultType& get() const noexcept { return *m_local_pointer; }

  private:
    ResultToLocalContextBuffer<ResultType>* m_buffer;
    ResultType* m_local_pointer;
};

} // namespace sammy

#endif
==> ./lns_destroy.h <==
#ifndef SAMMY_LNS_DESTROY_H_INCLUDED_
#define SAMMY_LNS_DESTROY_H_INCLUDED_

#include "barrage_lns_subproblem.h"
#include "fast_clique.h"
#include "implied_vertices.h"
#include "output.h"
#include "partial_solution.h"
#include "thread_interrupt.h"

namespace sammy {

/**
 * When destroying a solution S by removal of
 * C configurations D, we call a destruction
 * 'doomed' if there is a MES of size C that
 * becomes uncovered by the removal.
 * We can argue that such MES M are uniquely covered,
 * i.e., consist only of interactions covered
 * by a single configuration in S.
 * This is due to the following:
 * - if a vertex of M is covered by two configurations in D,
 *   that leads to a contradiction with the fact that
 *   the vertices in M are mutually exclusive and M has size C.
 * - if a vertex of M is covered by a configuration in D and
 *   another configuration not in D, that that vertex is not uncovered
 *   in the destroyed solution.
 * So, for avoidance of doomed destructions, looking for
 * uniquely covered MES is a good idea.
 */

/**
 * Info on a mutually exclusive set on
 * uniquely covered vertices
 * used to identify doomed destructions
 * during the LNS destroy procedure.
 */
struct LNSDestroyUniquelyCoveredMESInfo {
    LNSDestroyUniquelyCoveredMESInfo(std::vector<Vertex> vertices,
                                     const PartialSolution& current_solution)
        : m_uniquely_covered_mes(std::move(vertices)) {
        for (const Vertex& v : m_uniquely_covered_mes) {
            std::optional<Index> cls;
            current_solution.find_covering_classes(
                v.first, v.second, [&](Index i) {
                    if (cls) {
                        throw std::logic_error("Uniquely covered MES vertex "
                                               "covered by two classes.");
                    }
                    cls.emplace(i);
                });
            if (!cls) {
                m_uncovered_vertices.push_back(v);
            } else {
                m_configurations_with_mes_vertex[*cls] = v;
            }
        }
    }

    /**
     * The list of all vertices in the MES; if the solution is a full solution,
     * these vertices are all uniquely covered.
     */
    std::vector<Vertex> m_uniquely_covered_mes;

    /**
     * The set of configuration indices which contain a vertex
     * from m_uniquely_covered_mes.
     */
    HashMap<Index, Vertex> m_configurations_with_mes_vertex;

    /**
     * If the solution is not a full solution, contains a list of
     * vertices from m_uniquely_covered_mes that are not currently covered.
     * Necessary to undo the removal of configurations.
     */
    std::vector<Vertex> m_uncovered_vertices;

    /**
     * On removal of configurations, update the indices.
     */
    void update_indices_on_removal(std::vector<Index> old_to_new) {
        HashMap<Index, Vertex> new_configs_with_vertex;
        for (const auto& entry : m_configurations_with_mes_vertex) {
            Index old_index = entry.first;
            if (old_to_new[old_index] == NIL) {
                m_uncovered_vertices.push_back(entry.second);
            } else {
                new_configs_with_vertex[old_to_new[old_index]] = entry.second;
            }
        }
        m_configurations_with_mes_vertex = std::move(new_configs_with_vertex);
    }

    /**
     * After the removal is undone on failure to find
     * an improved solution, restore the uncovered vertices.
     */
    void undo_removal(const PartialSolution& restored) {
        for (const Vertex& v : m_uncovered_vertices) {
            m_configurations_with_mes_vertex[restored.find_covering_class(
                v.first, v.second)] = v;
        }
        m_uncovered_vertices.clear();
    }

    /**
     * Check if this MES dooms the given destruction to fail.
     */
    bool dooms(const std::vector<Index>& destruction_indices) const {
        return std::all_of(destruction_indices.begin(),
                           destruction_indices.end(), [this](Index i) -> bool {
                               return m_configurations_with_mes_vertex.count(i);
                           });
    }

    void remove_included(std::vector<Index>& indices) const {
        indices.erase(
            std::remove_if(indices.begin(), indices.end(),
                           [this](Index i) {
                               return m_configurations_with_mes_vertex.count(i);
                           }),
            indices.end());
    }
};

/**
 * LNS destroy operation of different types.
 */
class LNSDestroy {
  public:
    LNSDestroy(EventRecorder* local_recorder, std::size_t worker_index,
               const ImpliedVertexCache* implied_cache,
               PartialSolution full_solution, std::vector<Vertex> global_mes,
               SharedDBPropagator propagator)
        : m_local_recorder(local_recorder), m_worker_index(worker_index),
          m_implied_cache(implied_cache),
          m_clique_builder(std::move(propagator)),
          m_current_subproblem(std::move(full_solution)),
          m_best_global_mes(std::move(global_mes)) {}

    std::size_t total_num_configs() const noexcept {
        return m_current_subproblem.size() + m_currently_removed.size();
    }

    void update_global_mes_if_better(const std::vector<Vertex>& global_mes) {
        if (global_mes.size() > m_best_global_mes.size()) {
            m_best_global_mes = global_mes;
        }
    }

    LNSSubproblem move_out_subproblem() const {
        Lit n_concrete =
            static_cast<Lit>(m_current_subproblem.get_n_concrete());
        return LNSSubproblem{std::move(m_filtered_uncovered),
                             std::move(m_best_initial_mes),
                             std::move(m_currently_removed), n_concrete};
    }

    /**
     * Get the current remaining assignments.
     */
    const PartialSolution& get_partial() const { return m_current_subproblem; }

    /**
     * Update the full solution that we want to destroy.
     * Among other side-effects, this invalidates our constraining MESs.
     */
    void
    update_full_solution_if_better(const PartialSolution& new_full_solution) {
        if (new_full_solution.size() <
            m_current_subproblem.size() + m_currently_removed.size())
        {
            m_currently_removed.clear();
            m_current_subproblem = new_full_solution;
            m_constraining_unique_mes.clear();
            m_superset_constraints.clear();
        }
    }

    void update_full_solution_if_better(PartialSolution&& new_full_solution) {
        if (new_full_solution.size() <
            m_current_subproblem.size() + m_currently_removed.size())
        {
            m_currently_removed.clear();
            m_current_subproblem = std::move(new_full_solution);
            m_constraining_unique_mes.clear();
            m_superset_constraints.clear();
        }
    }

    /**
     * Integrate a new MES into the constraints.
     * Must be a 'uniquely covered' MES.
     */
    void integrate_mes(const std::vector<Vertex>& mes) {
        m_constraining_unique_mes.emplace_back(mes, m_current_subproblem);
    }

    /**
     * Called on success to return the subproblem and the improved solution.
     */
    void return_subproblem_on_success(LNSSubproblem&& returned,
                                      PartialSolution&& improved) {
        m_currently_removed = std::move(returned.removed_configurations);
        m_filtered_uncovered = std::move(returned.uncovered_universe);
        m_best_initial_mes = std::move(returned.mutually_exclusive_set);
        m_currently_removed.clear();
        m_current_subproblem = std::move(improved);
        m_constraining_unique_mes.clear();
        m_superset_constraints.clear();
    }

    /**
     * Called on abort (or lost race against other solver) to return the
     * subproblem and undo the destruction.
     */
    void return_subproblem_on_abort(LNSSubproblem&& returned) noexcept {
        m_currently_removed = std::move(returned.removed_configurations);
        m_filtered_uncovered = std::move(returned.uncovered_universe);
        m_best_initial_mes = std::move(returned.mutually_exclusive_set);
        p_undo_removal();
    }

    std::size_t worker_index() const noexcept { return m_worker_index; }

    /**
     * Called on impossible improvement to
     * return the subproblem and undo the destruction.
     */
    void improvement_impossible(LNSSubproblem&& returned,
                                const std::vector<Vertex>& mes) {
        m_currently_removed = std::move(returned.removed_configurations);
        m_filtered_uncovered = std::move(returned.uncovered_universe);
        m_best_initial_mes = std::move(returned.mutually_exclusive_set);
        m_best_initial_mes.clear();
        if (mes.size() >= m_currently_removed.size()) {
            integrate_mes(mes);
        } else {
            std::vector<Index> indices(m_currently_removed.size(), 0);
            std::iota(indices.begin(), indices.end(),
                      Index(m_current_subproblem.size()));
            m_superset_constraints.push_back(std::move(indices));
        }
        p_undo_removal();
    }

    /**
     * Trigger the destroy operation with a given number
     * of goal configurations to remove.
     * Throws an interrupt exception if the operation
     * detects the interruption flag being set.
     * After that exception, the class remains in a usable state.
     * Returns number of removed configurations;
     * return 0 iff the solution is provably optimal and cannot be improved.
     */
    std::size_t destroy(std::size_t goal_removed);

    /**
     * After finding a replacement for the destroyed
     * part of the solution, make a new full solution.
     */
    PartialSolution
    improve_destroyed(const std::vector<DynamicBitset>& replacement) {
        PartialSolution result = m_current_subproblem;
        for (const auto& bs : replacement) {
            result.add_assignment(bs);
        }
        return result;
    }

    const std::vector<Vertex>& best_global_mes() const noexcept {
        return m_best_global_mes;
    }

  private:
    /**
     * Event recorder for destruction events.
     */
    EventRecorder* m_local_recorder;

    /**
     * The worker index to add to events for identification.
     */
    std::size_t m_worker_index;

    /**
     * Probability for destruction by random + least-uniquely-covered removal.
     */
    double m_prob_random_least_unique = 0.5;

    /**
     * Probability for purely random destruction.
     */
    double m_prob_random_destruction = 0.2;

    /**
     * Probability for destruction avoiding doom.
     */
    // double m_prob_avoid_doom = 0.3;

    /**
     * Cache storing implied vertices.
     */
    const ImpliedVertexCache* m_implied_cache;

    /**
     * Local clique heuristic to build MESs.
     */
    FastCliqueBuilder m_clique_builder;

    /**
     * PartialSolution encoding the (remaining)
     * subproblem (before destruction, a solution
     * that we want to destroy).
     */
    PartialSolution m_current_subproblem;

    /**
     * List of currently removed configurations
     * (i.e., configurations we removed to obtain m_current_subproblem from a
     * full solution).
     */
    std::vector<DynamicBitset> m_currently_removed;

    /**
     * Uniquely-covered MESs constrain our destructions
     * to avoid doomed destructions.
     */
    std::vector<LNSDestroyUniquelyCoveredMESInfo> m_constraining_unique_mes;

    /**
     * Filtered uncovered vertices after destruction.
     */
    std::vector<Vertex> m_filtered_uncovered;

    /**
     * Best initial MES found for the current destruction.
     */
    std::vector<Vertex> m_best_initial_mes;

    /**
     * The best global MES.
     */
    std::vector<Vertex> m_best_global_mes;

    /**
     * Sets for which we have established that no improvement is possible.
     */
    std::vector<std::vector<Index>> m_superset_constraints;

    /**
     * Attempt to destroy actual_to_remove classes purely at random
     * for at most num_tries tries.
     * Returns true if the destruction looks successful, i.e., does not
     * violate any constraints.
     */
    inline bool p_destroy_random(Index actual_to_remove, std::size_t num_tries);

    /**
     * Destroy mostly randomly but use the MESs to avoid doomed destructions
     * as much as possible.
     */
    inline bool p_destroy_avoid_doom(Index actual_to_remove,
                                     std::size_t num_tries);

    /**
     * Destroy actual_to_remove classes by a mix of random and
     * least-unique-covered.
     */
    inline bool p_destroy_random_and_least_unique(Index actual_to_remove,
                                                  std::size_t num_tries);

    /**
     * Check if the given list of indices is a doomed removal.
     */
    inline bool
    p_removal_doomed(const std::vector<Index>& removal_indices) const;

    /**
     * Check if the current removal is doomed by a MES;
     * if so, undo it and return true, else return false.
     */
    inline bool p_undo_if_doomed();

    /**
     * Actually perform the removal of the given indices.
     * This may be undone later using the m_currently_removed list,
     * but the indices may change.
     * Updates all indices of configurations accordingly.
     */
    inline void p_perform_removal(const std::vector<Index>& removal_indices,
                                  bool enum_uncovered = true);

    /**
     * Undo a removal if the removal did not lead to an improved solution.
     */
    inline void p_undo_removal();

    /**
     * Run the heuristic MES after removal of the given indices.
     */
    inline bool p_removal_heuristic_mes_dooms();

    /**
     * Consider the given MES (which may include covered interactions if
     * all_known_uncovered = false).
     */
    inline bool p_heuristic_mes_from_dooms(std::vector<Vertex>& mes,
                                           bool all_known_uncovered,
                                           bool integrate_unextended);

    template <typename... Args>
    inline void p_store_event(const char* event_name, OutputObject object,
                              Args&&... args) {
        if (!m_local_recorder)
            return;
        object["worker_index"] = m_worker_index;
        m_local_recorder->store_event(event_name, std::move(object),
                                      std::forward<Args>(args)...,
                                      "worker_index");
    }
};

// -------------------- IMPLEMENTATION --------------------
std::size_t LNSDestroy::destroy(std::size_t goal_removed) {
    if (goal_removed < 3) {
        goal_removed = 3;
    }
    p_store_event("LNS_DESTROY_BEGIN", {{"goal_removed", goal_removed}},
                  "goal_removed");
    auto& rng = sammy::rng();
    if (!m_currently_removed.empty()) {
        throw std::logic_error(
            "Entering destroy with partially-destroyed solution!");
    }
    if (goal_removed >= m_current_subproblem.size()) {
        goal_removed = m_current_subproblem.size();
        p_store_event("LNS_DESTROY_WHOLE_SOLUTION", {});
        std::vector<Index> indices = vector(range(Index(goal_removed)));
        if (p_removal_doomed(indices)) {
            p_store_event("LNS_DESTROY_SOLUTION_OPTIMAL", {});
            return 0;
        }
        p_perform_removal(indices);
        if (p_removal_heuristic_mes_dooms()) {
            p_store_event("LNS_DESTROY_SOLUTION_OPTIMAL", {});
            return 0;
        }
        p_store_event("LNS_DESTROY_WHOLE_SOLUTION_DONE", {});
        return goal_removed;
    }
    double x = std::uniform_real_distribution<double>(0.0, 1.0)(rng);
    if (x < m_prob_random_least_unique) {
        if (p_destroy_random_and_least_unique(goal_removed, 5)) {
            return goal_removed;
        }
    } else if (x < m_prob_random_destruction + m_prob_random_least_unique) {
        if (p_destroy_random(goal_removed, 10)) {
            return goal_removed;
        }
    }
    if (p_destroy_avoid_doom(goal_removed, 10)) {
        return goal_removed;
    }
    // TODO: current fallback is increasing the subproblem size
    return destroy(goal_removed + 1);
    // return p_destroy_fallback(goal_removed);
}

bool LNSDestroy::p_destroy_random(Index actual_to_remove,
                                  std::size_t num_tries) {
    for (std::size_t i = 0; i < num_tries; ++i) {
        throw_if_interrupted();
        p_store_event("LNS_DESTROY_RANDOM_BEGIN_TRIAL",
                      {{"trial", i + 1}, {"goal_removed", actual_to_remove}},
                      "goal_removed", "trial");
        std::vector<Index> indices = vector(range(actual_to_remove));
        std::shuffle(indices.begin(), indices.end(), sammy::rng());
        indices.resize(actual_to_remove);
        std::sort(indices.begin(), indices.end());
        if (p_removal_doomed(indices))
            continue;
        p_perform_removal(indices);
        if (p_removal_heuristic_mes_dooms()) {
            p_undo_removal();
            continue;
        }
        p_store_event("LNS_DESTROY_RANDOM_DONE",
                      {{"goal_removed", actual_to_remove}}, "goal_removed");
        return true;
    }
    p_store_event(
        "LNS_DESTROY_RANDOM_FAILED",
        {{"goal_removed", actual_to_remove}, {"num_trials", num_tries}},
        "goal_removed");
    return false;
}

bool LNSDestroy::p_destroy_avoid_doom(Index actual_to_remove,
                                      std::size_t num_tries) {
    auto& rng = sammy::rng();
    DynamicBitset available(m_current_subproblem.size(), true);
    DynamicBitset superset_buffer(m_current_subproblem.size(), false);
    std::vector<Index> buffer;
    std::vector<Index> current_removal;

    auto available_to_buffer = [&](const auto& constraint) {
        buffer.clear();
        for (Index i : range(Index(m_current_subproblem.size()))) {
            if (available[i] &&
                !constraint.m_configurations_with_mes_vertex.count(i))
            {
                buffer.push_back(i);
            }
        }
    };

    auto superset_available_to_buffer = [&](const auto& constraint) {
        buffer.clear();
        for (Index i : range(Index(m_current_subproblem.size()))) {
            if (available[i] && !superset_buffer[i]) {
                buffer.push_back(i);
            }
        }
    };

    auto superset_dooms = [&](const auto& constraint) {
        superset_buffer.reset();
        for (Index i : constraint) {
            superset_buffer[i].set();
        }
        for (Index i : current_removal) {
            if (!superset_buffer[i]) {
                return false;
            }
        }
        return true;
    };

    auto random_extend = [&]() -> bool {
        if (buffer.empty()) {
            return false;
        }
        Index next_to_remove = buffer[std::uniform_int_distribution<Index>(
            0, buffer.size() - 1)(rng)];
        current_removal.push_back(next_to_remove);
        available[next_to_remove] = false;
        if (current_removal.size() > actual_to_remove) {
            return false;
        }
        return true;
    };

    auto run_try = [&](std::size_t trial) -> bool {
        throw_if_interrupted();
        p_store_event("LNS_DESTROY_AVOID_DOOM_BEGIN_TRIAL",
                      {{"goal_removed", actual_to_remove}, {"trial", trial}},
                      "goal_removed", "trial");
        current_removal.clear();
        available.set();
        for (const auto& constraint : m_constraining_unique_mes) {
            if (constraint.dooms(current_removal)) {
                available_to_buffer(constraint);
                if (!random_extend())
                    return false;
            }
        }
        for (const auto& constraint : m_superset_constraints) {
            if (superset_dooms(constraint)) {
                superset_available_to_buffer(constraint);
                if (!random_extend())
                    return false;
            }
        }
        std::size_t remaining_to_remove =
            actual_to_remove - current_removal.size();
        if (remaining_to_remove > 0) {
            buffer.clear();
            for (Index i : range(Index(m_current_subproblem.size()))) {
                if (available[i]) {
                    buffer.push_back(i);
                }
            }
            std::shuffle(buffer.begin(), buffer.end(), rng);
            buffer.resize(remaining_to_remove);
            current_removal.insert(current_removal.end(), buffer.begin(),
                                   buffer.end());
        }
        std::sort(current_removal.begin(), current_removal.end());
        p_store_event("LNS_DESTROY_AVOID_DOOM_DONE",
                      {{"goal_removed", actual_to_remove}, {"trial", trial}},
                      "goal_removed", "trial");
        return true;
    };

    for (std::size_t i = 0; i < num_tries; ++i) {
        if (run_try(i + 1)) {
            p_perform_removal(current_removal);
            if (p_removal_heuristic_mes_dooms()) {
                p_undo_removal();
                continue;
            }
            return true;
        }
    }
    return false;
}

bool LNSDestroy::p_undo_if_doomed() {
    for (const auto& info : m_constraining_unique_mes) {
        if (info.m_uncovered_vertices.size() >= m_currently_removed.size()) {
            p_undo_removal();
            return true;
        }
    }
    return false;
}

bool LNSDestroy::p_destroy_random_and_least_unique(Index goal_removed,
                                                   std::size_t num_tries) {
    std::size_t goal_removed_random =
        (std::min)(goal_removed, goal_removed / 3 + 1);
    std::size_t goal_removed_least_unique = goal_removed - goal_removed_random;
    if (goal_removed_least_unique == 0) {
        return p_destroy_random(goal_removed_random, num_tries);
    }
    p_store_event("BEGIN_DESTROY_RANDOM_AND_LEAST_UNIQUE",
                  {{"goal_removed", goal_removed},
                   {"random", goal_removed_random},
                   {"least_unique", goal_removed_least_unique}},
                  "goal_removed");
    for (std::size_t t = 0; t < num_tries; ++t) {
        throw_if_interrupted();
        p_store_event("BEGIN_DESTROY_RANDOM_AND_LEAST_UNIQUE_TRIAL",
                      {{"goal_removed", goal_removed},
                       {"random", goal_removed_random},
                       {"least_unique", goal_removed_least_unique},
                       {"trial", t + 1}},
                      "goal_removed", "trial");
        std::vector<Index> indices =
            vector(range(Index(m_current_subproblem.size())));
        std::shuffle(indices.begin(), indices.end(), sammy::rng());
        indices.resize(goal_removed_random);
        std::sort(indices.begin(), indices.end());
        p_perform_removal(indices, /*enum_uncovered=*/false);
        std::vector<std::size_t> unique_per_class(m_current_subproblem.size(),
                                                  0);
        m_current_subproblem.iterate_all_uniquely_covered_with_class(
            [&](Index i, Lit, Lit) { unique_per_class[i] += 1; });
        indices.clear();
        auto new_range = range(Index(m_current_subproblem.size()));
        indices.assign(new_range.begin(), new_range.end());
        std::nth_element(indices.begin(),
                         indices.begin() + goal_removed_least_unique,
                         indices.end(), [&](Index a, Index b) {
                             return unique_per_class[a] < unique_per_class[b];
                         });
        indices.resize(goal_removed_least_unique);
        std::sort(indices.begin(), indices.end());
        p_perform_removal(indices, /*enum_uncovered=*/true);
        if (p_undo_if_doomed()) {
            continue;
        }
        if (p_removal_heuristic_mes_dooms()) {
            p_undo_removal();
            continue;
        }
        p_store_event("DONE_DESTROY_RANDOM_AND_LEAST_UNIQUE_DONE",
                      {{"goal_removed", goal_removed},
                       {"random", goal_removed_random},
                       {"least_unique", goal_removed_least_unique},
                       {"trial", t + 1},
                       {"filtered_uncovered", m_filtered_uncovered.size()}},
                      "goal_removed", "trial", "filtered_uncovered");
        return true;
    }
    return false;
}

bool LNSDestroy::p_removal_doomed(
    const std::vector<Index>& removal_indices) const {
    if (std::any_of(m_constraining_unique_mes.begin(),
                    m_constraining_unique_mes.end(),
                    [&](const LNSDestroyUniquelyCoveredMESInfo& info) -> bool {
                        return info.dooms(removal_indices);
                    }))
    {
        return true;
    }
    auto is_subseteq = [&](const std::vector<Index>& superset) {
        return std::includes(superset.begin(), superset.end(),
                             removal_indices.begin(), removal_indices.end());
    };
    return std::any_of(m_superset_constraints.begin(),
                       m_superset_constraints.end(), is_subseteq);
}

void LNSDestroy::p_perform_removal(const std::vector<Index>& removal_indices,
                                   bool enum_uncovered) {
    std::vector<Index> old_to_new(m_current_subproblem.size(), NIL);
    auto current = removal_indices.begin(), end = removal_indices.end();
    Index new_index = 0;
    for (Index config_index : range(Index(m_current_subproblem.size()))) {
        if (config_index == *current) {
            if (++current == end) {
                for (Index c2 = config_index + 1,
                           n = m_current_subproblem.size();
                     c2 < n; ++c2)
                {
                    old_to_new[c2] = new_index++;
                }
                break;
            }
            continue;
        } else {
            old_to_new[config_index] = new_index++;
        }
    }
    for (auto& info : m_constraining_unique_mes) {
        info.update_indices_on_removal(old_to_new);
    }
    // handle already-removed indices
    for (Index x = 0; x < m_currently_removed.size(); ++x) {
        old_to_new.push_back(new_index++);
    }
    for (Index i : removal_indices) {
        old_to_new[i] = new_index++;
    }
    for (auto& sup : m_superset_constraints) {
        std::transform(sup.begin(), sup.end(), sup.begin(),
                       [&](Index i) { return old_to_new[i]; });
        std::sort(sup.begin(), sup.end());
    }
    if (!std::is_sorted(removal_indices.begin(), removal_indices.end())) {
        throw std::logic_error("removal_indices not sorted");
    }
    auto cp = removal_indices;
    if (std::unique(cp.begin(), cp.end()) != cp.end()) {
        throw std::logic_error("removal_indices not unique");
    }
    if (std::any_of(removal_indices.begin(), removal_indices.end(),
                    [&](Index i) { return i >= m_current_subproblem.size(); }))
    {
        throw std::logic_error("removal_indices out of bounds");
    }
    for (auto& r : m_current_subproblem.remove_assignments(removal_indices)) {
        m_currently_removed.push_back(std::move(r));
    }
    if (enum_uncovered) {
        try {
            m_filtered_uncovered.clear();
            m_current_subproblem.iterate_all_uncovered(
                [this](Lit lmin, Lit lmax) {
                    m_filtered_uncovered.emplace_back(lmin, lmax);
                },
                /*interruptible=*/true);
            m_implied_cache->remove_implied(m_filtered_uncovered,
                                            m_clique_builder.propagator());
        } catch (const InterruptError& e) {
            p_undo_removal();
            throw e;
        }
    }
}

void LNSDestroy::p_undo_removal() {
    m_filtered_uncovered.clear();
    for (auto& assignment : m_currently_removed) {
        m_current_subproblem.add_assignment(std::move(assignment));
    }
    m_currently_removed.clear();
    for (auto& info : m_constraining_unique_mes) {
        info.undo_removal(m_current_subproblem);
    }
}

bool LNSDestroy::p_heuristic_mes_from_dooms(std::vector<Vertex>& mes,
                                            bool all_known_uncovered,
                                            bool integrate_unextended) {
    if (!all_known_uncovered) {
        // if necessary, remove covered interactions
        auto is_covered = [this](Vertex v) {
            return !std::binary_search(m_filtered_uncovered.begin(),
                                       m_filtered_uncovered.end(), v);
        };
        mes.erase(std::remove_if(mes.begin(), mes.end(), is_covered),
                  mes.end());
    }
    if (mes.size() >= m_currently_removed.size()) {
        if (integrate_unextended) {
            // if this mes is not already in place, we can add it
            integrate_mes(mes);
        }
        return true;
    }
    mes = m_clique_builder.compute_clique_known_valid(mes.begin(), mes.end(),
                                                      m_filtered_uncovered);
    if (mes.size() >= m_currently_removed.size()) {
        integrate_mes(mes);
        return true;
    }
    if (mes.size() > m_best_initial_mes.size()) {
        m_best_initial_mes = mes;
    }
    return false;
}

bool LNSDestroy::p_removal_heuristic_mes_dooms() {
    auto mes = m_clique_builder.random_multistart_best_clique_known_valid(
        10, m_filtered_uncovered);
    if (mes.size() >= m_currently_removed.size()) {
        integrate_mes(mes);
        p_store_event("INITIAL_MES_PRECLUDES_IMPROVEMENT",
                      {{"removed_configs", m_currently_removed.size()},
                       {"mes_size", mes.size()},
                       {"mes_source", "random_multistart_clique"}},
                      "removed_configs", "mes_size", "mes_source");
        return true;
    }
    m_best_initial_mes = mes;
    mes.clear();
    for (const auto& info : m_constraining_unique_mes) {
        if (info.m_uncovered_vertices.size() > mes.size()) {
            mes = info.m_uncovered_vertices;
        }
    }
    if (!mes.empty() && p_heuristic_mes_from_dooms(mes, true, false)) {
        p_store_event("UNIQUELY_COVERED_MES_PRECLUDES_IMPROVEMENT",
                      {{"removed_configs", m_currently_removed.size()},
                       {"mes_size", mes.size()},
                       {"mes_source", "cached"}},
                      "removed_configs", "mes_size", "mes_source");
        return true;
    }
    std::vector<Vertex> tmp = m_best_global_mes;
    if (!m_best_global_mes.empty() &&
        p_heuristic_mes_from_dooms(tmp, false, true))
    {
        p_store_event("UNIQUELY_COVERED_MES_PRECLUDES_IMPROVEMENT",
                      {{"removed_configs", m_currently_removed.size()},
                       {"mes_source", "global best locally extended"}},
                      "removed_configs", "mes_source");
        return true;
    }
    return false;
}

} // namespace sammy

#endif
==> ./greedysat.h <==
#ifndef SAMMY_GREEDYSAT_H_INCLUDED_
#define SAMMY_GREEDYSAT_H_INCLUDED_

#include "literals.h"
#include "shared_db_propagator.h"

namespace sammy {

/**
 * A very simple greedy SAT solver based on SharedDBPropagator.
 * Whenever a concrete feature is unassigned, it picks an
 * arbitrary unassigned literal l and assigns it to the preferred
 * polarity.
 * The preferred polarity is given by PreferredAssignmentFn::operator[](l).
 * When all concrete features are assigned, any open non-concrete features
 * (of which there usually should be none) are greedily set to false
 * one after another, with propagation in between.
 */
template <typename PreferredAssignmentFn> class GreedySAT {
  public:
    GreedySAT(SharedDBPropagator* prop, Lit n_concrete,
              PreferredAssignmentFn&& assignment)
        : prop(prop),
          preferred_assignment(std::forward<PreferredAssignmentFn>(assignment)),
          n_concrete(n_concrete) {}

    bool solve() {
        while (prop->get_trail().size() < prop->db().num_vars()) {
            bool changed = false;
            for (Lit v = 0; v < n_concrete; ++v) {
                Lit l = preferred_assignment[v] ? lit::positive_lit(v)
                                                : lit::negative_lit(v);
                if (prop->is_open(l)) {
                    changed = true;
                    if (!prop->push_level(l) && !prop->resolve_conflicts()) {
                        throw UNSATError();
                    }
                }
            }
            if (changed)
                continue;
            Lit n_all = prop->db().num_vars();
            for (Lit v = n_concrete; v < n_all; ++v) {
                Lit l = lit::negative_lit(v);
                if (prop->is_open(l)) {
                    if (!prop->push_level(l)) {
                        prop->resolve_or_throw();
                    }
                }
            }
        }
        return true;
    }

  private:
    SharedDBPropagator* prop;
    PreferredAssignmentFn preferred_assignment;
    Lit n_concrete;
};

template <bool B> struct PreferValue {
    bool operator[](Lit /*var*/) const noexcept { return B; }
};
using PreferFalse = PreferValue<false>;
using PreferTrue = PreferValue<true>;

} // namespace sammy

#endif
==> ./algorithm_ex.h <==
#ifndef SAMMY_ALGORITHM_EX_H_INCLUDED_
#define SAMMY_ALGORITHM_EX_H_INCLUDED_

#include "range.h"
#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>

namespace sammy {

template <typename ForwardIterator, typename Callable>
static inline ForwardIterator
find_pair_if(ForwardIterator begin, ForwardIterator end, Callable&& callable) {
    ForwardIterator last = begin;
    if (begin != end) {
        for (++begin; begin != end; ++begin) {
            if (callable(*last, *begin)) {
                return last;
            }
            last = begin;
        }
    }
    return end;
}

template <typename InputIterator, typename OutputIterator,
          typename Transformation, typename Predicate>
static inline OutputIterator
copy_transformed_if(InputIterator begin, InputIterator end, OutputIterator out,
                    Transformation&& transform, Predicate&& predicate) {
    begin = std::find_if(begin, end, std::forward<Predicate>(predicate));
    while (begin != end) {
        *out = transform(*begin);
        ++begin;
        ++out;
        begin = std::find_if(begin, end, std::forward<Predicate>(predicate));
    }
    return out;
}

template <typename Container>
static inline void sort_unique(Container& container) {
    std::sort(container.begin(), container.end());
    container.erase(std::unique(container.begin(), container.end()),
                    container.end());
}

template <typename SetType, typename Container>
static inline void nosort_unique(Container& container) {
    using Reference = typename Container::const_reference;
    SetType set;
    auto new_end = std::remove_if(
        container.begin(), container.end(),
        [&set](Reference elem) { return !set.insert(elem).second; });
    container.erase(new_end, container.end());
}

template <typename Container, typename InputIterator>
static inline void add_and_make_unique(Container& container,
                                       InputIterator begin, InputIterator end) {
    container.insert(container.end(), begin, end);
    std::sort(container.begin(), container.end());
    container.erase(std::unique(container.begin(), container.end()),
                    container.end());
}

template <typename ForwardIterator, typename IndexForwardIterator>
static inline ForwardIterator remove_indices(ForwardIterator begin,
                                             ForwardIterator end,
                                             IndexForwardIterator remove_begin,
                                             IndexForwardIterator remove_end) {
    if (remove_begin == remove_end)
        return end;
    auto first_removed = *remove_begin++;
    ForwardIterator out = begin;
    std::advance(out, first_removed);
    ForwardIterator in = out;
    ++in;
    auto current_index = first_removed;
    ++current_index;
    while (remove_begin != remove_end) {
        for (auto next_removed = *remove_begin; current_index < next_removed;
             ++current_index, ++in, ++out)
        {
            *out = std::move(*in);
        }
        ++in;
        ++current_index;
        ++remove_begin;
    }
    for (; in != end; ++in, ++out) {
        *out = std::move(*in);
    }
    return out;
}

template <typename Range> static inline auto vector(const Range& range) {
    using std::begin;
    using std::end;
    using Iterator = decltype(begin(range));
    using ValueType = typename std::iterator_traits<Iterator>::value_type;
    return std::vector<ValueType>{begin(range), end(range)};
}

template <typename Iterator>
static inline auto vector(Iterator begin, Iterator end) {
    using ValueType = typename std::iterator_traits<Iterator>::value_type;
    return std::vector<ValueType>{begin, end};
}

template <typename IntType>
static inline auto range(IntType begin, IntType end) {
    end = (std::max)(begin, end);
    return IteratorRange{boost::counting_iterator<IntType>(begin),
                         boost::counting_iterator<IntType>(end)};
}

template <typename IntType> static inline auto range(IntType end) {
    return range(IntType(0), end);
}

// helpers for std::variant visitation
template <class... Ts> struct overloaded : Ts... {
    using Ts::operator()...;
};
template <class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

} // namespace sammy

#endif
==> ./clique_sat_dsatur.h <==
#ifndef SAMMY_CLIQUE_SAT_DSATUR_H_INCLUDED_
#define SAMMY_CLIQUE_SAT_DSATUR_H_INCLUDED_

#include "algorithm_ex.h"
#include "best_k.h"
#include "class_completer.h"
#include "dynamic_bitset.h"
#include "gurobi.h"
#include "literals.h"
#include "output.h"
#include "pair_infeasibility_map.h"
#include "partial_solution.h"
#include "thread_interrupt.h"
#include "vertex_operations.h"
#include <boost/iterator/transform_iterator.hpp>
#include <variant>

namespace sammy {

/**
 * The maximum number of new constraints added in
 * each call of p_prohibit_violated_nonedges.
 */
static constexpr std::size_t CAP_PROHIBIT_VIOLATED_NUM_NEW_CONSTRAINTS = 20;

/**
 * The number of cheap cut rounds to do before
 * performing a check to see whether the gap
 * has shrunk enough to continue cheap cuts.
 */
static constexpr std::size_t CHEAP_CUT_ROUNDS_PER_GAP_CHECK = 10;

/**
 * Gap: (objective - m_best_clique.size()) / m_best_clique.size()
 * The gap reduction required between two groups of
 * CHEAP_CUT_ROUNDS_PER_GAP_CHECK cheap cut rounds to
 * continue cheap cuts.
 */
static constexpr double CHEAP_CUT_GAP_REDUCTION_REQUIRED = 0.02;

/**
 * The maximum number of cheap cut rounds to do
 * before giving up on cheap cuts for the time being.
 */
static constexpr std::size_t CHEAP_CUT_ROUNDS_HARD_CAP = 100;

/**
 * The number of greedy cuts to add in each round of
 * greedy cut generation.
 */
static constexpr std::size_t GREEDY_CUT_CANDIDATE_SET_SIZE = 7;

template <typename IncrementalSATSolver> class CliqueSatDSaturSolver {
  public:
    /**
     * @brief Constructs a new CliqueSatDSaturSolver object.
     * @param considered_vertices All vertices to consider.
     * @param infeasibility_map The pair infeasibility map.
     *                          Is not changed. learn_infeasibilities should be
     * called before.
     * @param clause_db The clause database. May be changed by learning conflict
     * clauses.
     * @param best_local_mes The best MES involving vertices from
     * considered_vertices.
     * @param best_global_mes The best MES involving all vertices.
     * @param best_global_lb The best lower bound on the number of configuration
     * needed.
     * @param covering_assignments A set of complete configurations covering all
     * vertices in considered_vertices.
     */
    CliqueSatDSaturSolver(std::vector<Vertex> considered_vertices,
                          PairInfeasibilityMap* infeasibility_map,
                          ClauseDB& clause_db,
                          const std::vector<Vertex>& best_local_mes,
                          std::size_t best_global_mes,
                          std::size_t best_global_lb,
                          std::vector<DynamicBitset> covering_assignments,
                          bool enable_incremental_clique = true);

    /**
     * Abort the solve (on timeout, optimality, LNS with better solution, ...).
     */
    void abort();

    /**
     * The size of the graph (all vertices) we are working on.
     */
    std::size_t get_graph_size() const noexcept {
        return m_all_vertices.size();
    }

    /**
     * Enum describing the outcome of the solve process.
     */
    enum class SolveResult {
        IMPROVED_SOLUTION, //< We found a new, optimal solution.
        ABORTED,           //< We were aborted (timeout, ...) before completion.
        SOLUTION_WAS_OPTIMAL //< We completed but the initial solution was
                             // optimal.
    };

    /**
     * Run the solution process.
     * @return SolveResult describing the outcome.
     */
    SolveResult solve();

    /**
     * Get the best solution found. Only valid if
     * the SolveResult was IMPROVED_SOLUTION.
     *
     * @return The best solution found.
     */
    PartialSolution get_partial_solution() const;

    /**
     * @return A reference to the covering assignments given initially.
     */
    const std::vector<DynamicBitset>&
    get_covering_assignments() const noexcept {
        return m_covering_assignments;
    }

    /**
     * Check whether we may have found an improvement
     * on the global lower bound.
     */
    bool improved_global_bound() const noexcept {
        return m_lower_bound > m_initial_global_lower_bound;
    }

    /**
     * Check whether we may have found an improvement
     * on the global MES size.
     */
    bool improved_mes() const noexcept {
        return m_best_clique.size() > m_initial_global_clique_size;
    }

    /**
     * Get a copy of the best clique found.
     */
    std::vector<Vertex> get_best_mes() const {
        std::vector<Vertex> result;
        result.reserve(m_best_clique.size());
        std::transform(m_best_clique.begin(), m_best_clique.end(),
                       std::back_inserter(result),
                       [&](std::size_t vi) { return m_all_vertices[vi]; });
        return result;
    }

    /**
     * Get the best lower bound on the number of configurations
     * needed for the subgraph.
     */
    std::size_t get_best_bound() const noexcept { return m_lower_bound; }

    /**
     * Set the event recorder to report events to;
     * this is purely optional.
     */
    void set_event_recorder(EventRecorder* recorder) noexcept {
        m_local_recorder = recorder;
    }

    /**
     * Get the vertices inducing the subgraph on which
     * we established the lower bound.
     */
    std::vector<Vertex> get_best_bound_subgraph() const {
        std::vector<Vertex> result;
        result.reserve(m_lower_bound_subgraph.size());
        std::transform(m_lower_bound_subgraph.begin(),
                       m_lower_bound_subgraph.end(), std::back_inserter(result),
                       [&](std::size_t vi) { return m_all_vertices[vi]; });
        return result;
    }

    /**
     * Set a callback to be invoked each time we try to search for a new clique.
     * It can return a clique (containing any vertices); the solver will filter
     * out all irrelevant vertices and only consider the ones in m_all_vertices.
     */
    void set_clique_candidate_callback(
        std::function<std::vector<Vertex>()> callback) {
        m_clique_candidate_callback = std::move(callback);
    }

    /**
     * Set a callback to be invoked each time this solver finds a new lower
     * bound. The callback is involved with the lower bound and the vertices of
     * the subgraph inducing it.
     */
    void set_lower_bound_callback(
        std::function<void(std::size_t, const std::vector<Vertex>&)> callback) {
        m_lower_bound_callback = std::move(callback);
    }

  private:
    /**
     * The literal type of the SAT solver.
     */
    using SatLit = typename IncrementalSATSolver::Lit;

    /**
     * A SAT literal or a fixed boolean value.
     */
    using LitOrFixed = std::variant<SatLit, bool>;

    /**
     * The pair infeasibility map.
     */
    PairInfeasibilityMap* m_inf_map;

    /**
     * The number of features.
     */
    Index m_n_all;

    /**
     * The number of concrete features.
     */
    Index m_n_concrete;

    /**
     * A propagator we can reuse.
     * After use, it must be emptied again.
     */
    SharedDBPropagator m_empty_propagator;

    /**
     * Optional callback that is called to generate clique
     * candidates (may contain vertices not in m_all_vertices,
     * which are then filtered out).
     */
    std::function<std::vector<Vertex>()> m_clique_candidate_callback;

    /**
     * Optional callback that is called each time we find
     * a new lower bound.
     */
    std::function<void(std::size_t, const std::vector<Vertex>&)>
        m_lower_bound_callback;

    /**
     * The set of all vertices we need to consider.
     * These can be reduced by removing implied vertices.
     */
    std::vector<Vertex> m_all_vertices;

    /**
     * For each concrete literal, the list of vertices containing that literal.
     */
    std::vector<std::vector<std::size_t>>
        m_vertices_containing_concrete_literal;

    /**
     * For each concrete literal, the list of vertex indices
     * into m_existing_clique_vars containing that literal
     * (that are present in the clique).
     */
    std::vector<std::vector<std::size_t>>
        m_clique_model_vertices_containing_concrete_literal;

    /**
     * For each (concrete or non-concrete) literal, the set of vertices
     * _implying_ that literal by propagation.
     */
    std::vector<DynamicBitset> m_vertices_implying_literal;

    /**
     * The adjacency matrix of the graph.
     * A 'true' means that there is an edge (i.e., conflict)
     * between two vertices (i.e., interactions).
     * However, a 'false' only means the edge is not there
     * if the corresponding entry in m_definitive_nonedges is true.
     * Initialized from conflicts according to m_vertices_implying_literal,
     * i.e., initially corresponding to graph G_1.
     */
    std::vector<DynamicBitset> m_adjacency_matrix;

    /**
     * A 'true' means that there is definitely no edge,
     * i.e., propagation of interactions m_all_vertices[i]
     * and m_all_vertices[j] at the same time found no conflict.
     */
    std::vector<DynamicBitset> m_definitive_nonedges;

    /**
     * Whether the incremental clique solver should be used.
     * If false (not the default), the solver will not be
     * run and the initial clique is used throughout the entire solve.
     */
    bool m_incremental_clique_turned_off;

    /**
     * The Gurobi clique solver model.
     */
    GRBModel m_clique_model;

    /**
     * Array of all (potential) clique vertex variables.
     * Index i corresponds to index i in m_all_vertices.
     * The variable is present if and only if the
     * corresponding bit in m_present_in_clique_model is set.
     */
    std::vector<GRBVar> m_all_clique_vars;

    /**
     * Array of all existing clique vertex variables.
     * Necessary for querying all values at once.
     */
    std::vector<GRBVar> m_existing_clique_vars;

    /**
     * Buffer for variable values from the last solution
     * of the clique solver model.
     */
    std::unique_ptr<double[]> m_existing_clique_var_values;
    std::size_t m_existing_clique_var_value_count = 0;
    std::string m_last_vertex_addition_source = "not added any";

    /**
     * The last relaxation objective value.
     */
    double m_last_clique_objective = -std::numeric_limits<double>::infinity();

    /**
     * The last relaxation objective value rounded to
     * an integer upper bound (up to some numeric slack).
     */
    std::size_t m_last_rounded_objective = 0;

    /**
     * Struct containing clique model index and value
     * for sorting vertices by their value in the last solution.
     */
    struct OrderedSolutionValue {
        std::size_t clique_model_index;
        double value;

        bool operator<(const OrderedSolutionValue& other) const noexcept {
            return value > other.value;
        }
    };

    /**
     * The non-zero elements of the last solution of the clique solver model,
     * sorted by the variable value in the solution.
     */
    std::vector<OrderedSolutionValue> m_clq_ordered_solution;

    using OrderedSolutionIter =
        typename std::vector<OrderedSolutionValue>::const_iterator;

    /**
     * Collection/buffer for violated non-edges
     * in the relaxation of the clique model.
     */
    std::vector<std::pair<OrderedSolutionIter, OrderedSolutionIter>>
        m_clq_violated_nonedges;

    /**
     * Index i has the vertex index of the vertex
     * referred to by m_existing_clique_vars[i].
     */
    std::vector<std::size_t> m_existing_clique_var_indices;

    /**
     * The set of all vertices present in the clique solver model.
     * Indices into m_all_vertices.
     */
    DynamicBitset m_present_in_clique_model;

    /**
     * The vertices that are currently colored.
     */
    std::vector<std::size_t> m_currently_colored;

    /**
     * The set of vertices that are in
     * m_currently_colored, i.e., are currently
     * explicitly colored.
     */
    DynamicBitset m_in_currently_colored;

    /**
     * The order of vertices in the SAT model.
     * I.e., m_sat_vertex_order[i] is the vertex index
     * in m_all_vertices corresponding to the i-th vertex
     * in the SAT model.
     * Reverse mapping via m_vertex_info[vi].sat_index.
     * Over time, this set only grows unless the model
     * is reset on finding a new clique.
     */
    std::vector<std::size_t> m_sat_vertex_order;

    /**
     * Information on all vertices for the DSatur-like
     * configuration assignment algorithm.
     */
    struct VertexInfo {
        DynamicBitset
            open_classes; //< a superset of the open classes for this vertex
        std::size_t num_open_classes; //< number of set bits in open_classes
        bool in_some_class;    //< number of classes containing this vertex
        std::size_t degree;    //< degree of the vertex in the graph
        std::size_t sat_index; //< index of the vertex in the SAT solver

        VertexInfo(std::size_t num_initial_classes, std::size_t ub_num_classes,
                   std::size_t degree)
            : num_open_classes{num_initial_classes}, in_some_class{false},
              degree{degree},
              sat_index{std::numeric_limits<std::size_t>::max()} {
            open_classes.reserve(ub_num_classes);
            open_classes.assign(num_initial_classes, true);
        }

        bool close_class(std::size_t class_index) {
            if (open_classes[class_index]) {
                open_classes[class_index].reset();
                --num_open_classes;
                return true;
            }
            return false;
        }
    };

    /**
     * Compare vertex indices by the VertexInfo.
     * Used for finding the best K
     * vertices according to their info.
     */
    struct CompareVertexIndexByInfo {
        inline bool operator()(std::size_t v1, std::size_t v2) const noexcept;
        CliqueSatDSaturSolver* that;
    };

    /**
     * The vertex info for each vertex.
     */
    std::vector<VertexInfo> m_vertex_info;

    /**
     * The best K vertex indices according to their info.
     */
    BestK<std::size_t, CompareVertexIndexByInfo> m_best_infos;

    /**
     * Information on a single color class (partial configuration).
     */
    struct ColorClass {
        SharedDBPropagator propagator;
        std::size_t index;

        ColorClass(const SharedDBPropagator& empty, std::size_t index)
            : propagator(empty), index(index) {}

        void add_vertex(CliqueSatDSaturSolver* that, std::size_t vertex);
        bool try_add_vertex(CliqueSatDSaturSolver* that, std::size_t vertex);

      private:
        void p_added_vertex(CliqueSatDSaturSolver* that,
                            std::size_t old_trail_end);
    };

    /**
     * The color classes.
     */
    std::vector<ColorClass> m_color_classes;

    /**
     * If we found any, the set of vertices
     * with zero possible configurations.
     */
    std::vector<std::size_t> m_zero_candidates;

    /**
     * If we found any, the set of vertices with only
     * one possible configuration. Last priority
     * before returning to the DSatur-like choice.
     */
    std::vector<std::size_t> m_one_candidates;

    /**
     * Gurobi callback that only checks for aborts.
     */
    class GurobiCallbackCheckAbort : public GRBCallback {
      public:
        explicit GurobiCallbackCheckAbort(CliqueSatDSaturSolver* that) noexcept
            : that(that) {}

        void callback() override;

      private:
        CliqueSatDSaturSolver* that;
    };
    std::unique_ptr<GurobiCallbackCheckAbort> m_check_abort_cb;

    /**
     * The incremental SAT solver.
     */
    IncrementalSATSolver m_incremental_solver;

    /**
     * SAT model variables.
     * If present, m_class_literals[c][i]
     * is the i-th variable of class c.
     */
    std::vector<std::vector<LitOrFixed>> m_class_literals;

    /**
     * Local recorder of events.
     */
    EventRecorder* m_local_recorder{nullptr};

    /**
     * SAT model variables.
     * If present, m_vertex_in_class[c][si]
     * is the variable indicating that the vertex
     * with vertex index m_sat_vertex_order[si] is in class c.
     */
    std::vector<std::vector<LitOrFixed>> m_vertex_in_class;

    /**
     * SAT model variable.
     * A variable used to disable the constraint that
     * a vertex has to be in one of classes 0, ... c - 1.
     * Fixed as soon as c + 1 colors are available.
     */
    std::optional<SatLit> m_not_enough_colors;

    /**
     * The number of colors m_not_enough_colors is for.
     */
    std::size_t m_not_enough_colors_is_for = 0;

    /**
     * Buffer for building SAT clauses.
     */
    std::vector<SatLit> m_clause_buffer;

    /**
     * Buffer for building vertex sets.
     */
    std::vector<std::size_t> m_vertex_set_buffer;

    /**
     * Buffer for building vertex bitsets.
     */
    DynamicBitset m_vertex_bitset_buffer;

    /**
     * Buffer for GRBLinExprs.
     */
    GRBLinExpr m_buffer_expr;

    /**
     * Buffer for GRBColumns.
     */
    GRBColumn m_column_buffer;

    /**
     * Buffer for building vertex index sets.
     */
    std::vector<std::size_t> m_vertices_buffer;

    /**
     * The best clique found on the vertex set.
     * Indices into m_all_vertices.
     */
    std::vector<std::size_t> m_best_clique;

    /**
     * The clique size last used by SAT.
     * If the clique size increases, we have to
     * rebuild the SAT model.
     */
    std::size_t m_last_sat_clique_size;

    /**
     * If the clique model contains vertices that are
     * not yet colored, there may be vertices in the
     * best clique found so far that are not colored yet.
     * This is the set of those vertices.
     *
     * These vertices are added to the SAT calls (to be used
     * to break symmetries).
     * An UNSAT result leaves these vertices uncolored.
     * At any time, these vertices have the highest priority,
     * disregarding the DSatur-like choice.
     */
    std::vector<std::size_t> m_in_clique_uncolored;

    /**
     * Vertices with zero leftover choices.
     */
    std::vector<std::size_t> m_zero_choice_vertices;

    /**
     * The initial size of the clique on m_all_vertices.
     */
    std::size_t m_initial_clique_size;

    /**
     * The size of the best-known clique (overall)
     * when we started the search.
     */
    std::size_t m_initial_global_clique_size;

    /**
     * The initial lower bound (overall) when
     * we started the search.
     */
    std::size_t m_initial_global_lower_bound;

    /**
     * The current lower bound.
     */
    std::size_t m_lower_bound;

    /**
     * The vertices inducing the subgraph on which
     * we established the lower bound.
     */
    std::vector<std::size_t> m_lower_bound_subgraph;

    /**
     * The upper bound (configurations covering all vertices
     * in m_all_vertices).
     */
    std::vector<DynamicBitset> m_covering_assignments;

    /**
     * Flag to abort the search.
     */
    std::atomic<bool> m_aborted{false};

    /**
     * Constraints in the Gurobi clique model:
     * complete assignments give us independent
     * sets that we use as constraints.
     */
    std::vector<DynamicBitset> m_clq_constr_complete_assignments;

    /**
     * Constraints in the Gurobi clique model
     * from the assignments in m_clq_constr_complete_assignments.
     */
    std::vector<GRBConstr> m_clq_constr_complete_constraints;

    /**
     * Dual values of the constraints in the Gurobi clique model
     * from the assignments in m_clq_constr_complete_assignments.
     */
    std::unique_ptr<double[]> m_clq_constr_complete_dual_values;

    /**
     * Constraints in the Gurobi clique model:
     * partial assignments give us independent
     * sets that we use as constraints.
     */
    std::vector<SharedDBPropagator> m_clq_constr_partial_assignments;

    /**
     * Constraints in the Gurobi clique model
     * from the assignments in m_clq_constr_partial_assignments.
     */
    std::vector<GRBConstr> m_clq_constr_partial_constraints;

    /**
     * Dual values of the constraints in the Gurobi clique model
     * from the assignments in m_clq_constr_partial_assignments.
     */
    std::unique_ptr<double[]> m_clq_constr_partial_dual_values;

    /**
     * Buffer for tracking the set of violated non-edges that we already
     * covered.
     */
    DynamicBitset m_covered_violated_nonedges;

    /**
     * Candidate for a greedily generated cutting plane.
     */
    struct GreedyCutCandidate {
        std::vector<std::size_t> vertex_indices;
        double cut_value;

        // is-better-than operator
        bool operator<(const GreedyCutCandidate& other) const noexcept {
            return cut_value > other.cut_value;
        }

        // is-better-than-score operator
        bool operator<(double cut_value) const noexcept {
            return this->cut_value > cut_value;
        }

        template <typename... VArgs>
        GreedyCutCandidate(double cut_value, VArgs&&... args)
            : vertex_indices(std::forward<VArgs>(args)...),
              cut_value(cut_value) {}
    };

    /**
     * Buffer for the best greedy cut candidates during greedy cut generation.
     */
    BestK<GreedyCutCandidate> m_greedy_cut_candidates;

    /**
     * PricingEntry: vertex and dual price + 1.0.
     */
    struct PricingEntry {
        PricingEntry(std::size_t v, double dw) noexcept
            : vertex_index(v), dual_weight(dw) {}
        std::size_t vertex_index;
        double dual_weight;
    };

    /**
     * The subset of 'good' vertices identified during pricing.
     * Good vertices are vertices for which the dual constraint
     * is violated by a margin of at least 0.01, i.e.,
     * the corresponding dual constraint's LHS is at most 0.99.
     */
    std::vector<PricingEntry> m_pricing_good_vertices;

    /**
     * The subset of 'possible' vertices identified during pricing.
     * Possible vertices are vertices for which the dual constraint
     * is close to being violated, i.e., the dual constraint's LHS
     * is between 0.99 and 1.01.
     */
    std::vector<PricingEntry> m_pricing_possible_vertices;

    // ---------- PRIVATE METHODS ----------
    /**
     * Call the lower bound callback with the current LB.
     */
    void p_call_lb_callback();

    /**
     * Initialize the adjacency matrices.
     * Called by the constructor.
     */
    void p_initialize_graph();

    /**
     * Initialize m_vertices_containing_concrete_literal.
     */
    void p_initialize_vertices_containing_concrete_literal();

    /**
     * Initialize m_vertices_implying_literal.
     */
    void p_initialize_vertices_implying_literal();

    /**
     * Initialize m_adjacency_matrix from
     * m_vertices_impliying_literal.
     */
    void p_initialize_matrix_from_implied_literals();

    /**
     * Transform a set of vertices to a set of indices.
     */
    std::vector<std::size_t> p_clique_to_indices(const std::vector<Vertex>& cvs,
                                                 bool must_match = true);

    /**
     * Extract the vertex indices into m_all_vertices that are
     * covered the given assignment into the vector vertices.
     */
    void
    p_extract_vertices_from_configuration(const DynamicBitset& assignment,
                                          std::vector<std::size_t>& vertices);

    /**
     * Extract the vertex indices into m_all_vertices that are
     * covered by the given assignment into m_vertices_buffer,
     * clearing m_vertices_buffer first.
     */
    void p_buffer_vertices_from_configuration(const DynamicBitset& assignment) {
        m_vertices_buffer.clear();
        p_extract_vertices_from_configuration(assignment, m_vertices_buffer);
    }

    enum class VertexChoiceType {
        ALL_COVERED,       //< all vertices are colored!
        FROM_CLIQUE,       //< from m_in_clique_uncolored
        WITH_ZERO_CHOICES, //< an element with zero choices
        WITH_ONE_CHOICE,   //< an element with one choice
        REGULAR_CHOICE     //< a regular DSatur-like choice
    };

    /**
     * A struct describing the next vertex to color.
     */
    struct VertexChoice {
        /**
         * The vertex to color.
         */
        std::size_t vertex;

        /**
         * Characterises the situation/choice.
         */
        VertexChoiceType type;
    };

    /**
     * Find the next vertex to color.
     */
    VertexChoice p_choose_next_vertex();

    /**
     * Check the clique for vertices to color.
     * Keeps the invariant that any vertex in m_best_clique
     * is either in m_currently_colored or in m_in_clique_uncolored.
     *
     * @param vc the vertex choice to fill in.
     * @return true if a vertex was found, false otherwise.
     */
    bool p_choose_next_vertex_color_clique(VertexChoice& vc);

    /**
     * Check m_one_candidates for vertices to color.
     *
     * @param vc the vertex choice to fill in.
     * @return true if a vertex was found, false otherwise.
     */
    bool p_choose_next_vertex_one_color_left(VertexChoice& vc);

    /**
     * Check all vertices for the best candidates
     * according to m_vertex_info.
     *
     * @return true if any uncolored vertex was found, false otherwise.
     */
    bool p_identify_best_candidates();

    /**
     * Check the best-looking candidates from m_best_infos
     * for their precise number of open classes.
     * @return the best vertex index according to our priority.
     */
    std::size_t p_precise_best_candidates();

    /**
     * Initialize the vertex info (m_vertex_info).
     */
    void p_initialize_vertex_info();

    /**
     * Initialize the color classes (m_color_classes).
     * Assigns the clique vertices to a color class each.
     */
    void p_initialize_color_classes();

    /**
     * Complete all color classes to full assignments.
     * @return true if all classes could be extended
     *         false if not all classes could be extended.
     *
     * If not all classes could be extended, the
     * SAT solver is called to find a complete assignment.
     * This may result in the lower bound increasing or
     * in a complete assignment that however does not
     * incorporate all vertices that were implicitly covered before.
     */
    bool p_complete_all_classes();

    /**
     * Handle the situation where the given vertex
     * has no possible choices.
     * @return true if continuing is sensible (in this case, the given vertex
     * will be added to m_currently_colored), false if we have proved that no
     * improvement is possible or were aborted (check m_aborted after return).
     */
    bool p_handle_zero_choice(std::size_t vertex);

    /**
     * Handle a given vertex.
     * Color it (and add it to m_currently_colored),
     * if that is possible.
     *
     * @return true if the vertex could be colored
     *         false if it turned out to have zero choices.
     */
    bool p_handle_vertex(std::size_t vertex);

    /**
     * Remove the given vertex from the stack
     * of vertices to be covered next, or the
     * stack of vertices with one candidate.
     */
    void p_colored_vertex(std::size_t vertex);

    /**
     * When we have a new clique, we completely
     * reset the SAT model.
     */
    void p_reset_sat_model();

    /**
     * Add a new color class, initially containing
     * the given vertex.
     */
    void p_add_color_class(std::size_t vertex);

    /**
     * When we need to run a new SAT/Clique model solve,
     * we use this method to compute the current set.
     */
    std::vector<std::size_t>
    p_compute_current_vertex_set(std::size_t new_vertex);

    /**
     * Compute/update the set of vertices that are
     * in the clique but not in m_currently_colored.
     */
    void p_compute_in_clique_uncolored();

    /**
     * Extend the given sorted unique set by the given elements;
     * the result is sorted & unique again.
     */
    void p_extend_sorted_unique_set(std::vector<std::size_t>& sorted_unique,
                                    const std::vector<std::size_t>& to_add);

    /**
     * Rebuild the SAT model.
     */
    void p_rebuild_sat_model(const std::vector<std::size_t>& current_vertices);

    /**
     * Rebuild a color class of the SAT model.
     */
    void p_rebuild_sat_color_class(std::optional<std::size_t> clique_vertex);

    /**
     * Append a vector of variables to the sat model to
     * obtain  the values of the features in a new color class.
     */
    std::vector<LitOrFixed>& p_append_unconstrained_class_lits();

    /**
     * Append a vector of variables to the sat model
     * to obtain the values of the features in a new color class,
     * initialized with the given vertex.
     */
    void p_append_constrained_class_lits(std::size_t clique_vertex);

    /**
     * Append a vector of variables to the sat model
     * to obtain the presence of vertices in a new color class.
     */
    std::vector<LitOrFixed>& p_append_vertex_in_class();

    /**
     * Lookup a sat literal or fixed value in a color class.
     */
    LitOrFixed p_lookup_sat_literal_in_class(Lit l, std::size_t class_index);

    /**
     * Negate a LitOrFixed.
     */
    static LitOrFixed p_negate(LitOrFixed l);

    /**
     * Return a LitOrFixed corresponding
     * to l1 & l2 (may add a new variable).
     */
    LitOrFixed p_and(LitOrFixed l1, LitOrFixed l2);

    /**
     * Add clauses to enforce that the given
     * class variables corresponds to a satisfying
     * assignment/valid configuration.
     */
    void p_enforce_formula_on(std::size_t class_index);

    /**
     * Add a clause for the given class variables.
     */
    void p_add_clause(const Lit* begin, const Lit* end,
                      std::size_t class_index);

    /**
     * Run SAT model.
     * @return true on SAT, false on UNSAT, nullopt on timeout/abort.
     */
    std::optional<bool> p_run_sat_model();

    /**
     * Extend the SAT model to all vertices in
     * current_vertices that are not yet part of the SAT model.
     */
    void p_extend_sat_model(const std::vector<std::size_t>& current_vertices);

    /**
     * Extend the SAT model to m_lower_bound color classes.
     */
    void p_extend_sat_model_to_classes();

    /**
     * Extract a satisfying assignment and
     * turn it into a partial coloring.
     * This may increase the number of members
     * in m_color_class to m_lower_bound.
     */
    template <typename ModelMapType>
    void p_extract_sat_model(const ModelMapType& model);

    /**
     * Add the constraint that the given vertex
     * must be in some color class.
     */
    void p_add_vertex_color_constraint(std::size_t sat_index,
                                       SatLit not_enough);

    /**
     * Check that m_not_enough_colors (and the vertex color constraints)
     * are set up for all vertices in m_sat_vertex_order.
     */
    void p_check_not_enough_colors();

    /**
     * Setup the clique model if necessary,
     * and solve it to hopefully find a larger clique.
     * Update m_lower_bound and m_best_clique if a better clique is found.
     *
     * @param current_vertices the vertices that are under consideration for the
     * SAT model
     * @return true if a larger clique was found, false otherwise
     */
    bool
    p_zero_choice_run_clique(const std::vector<std::size_t>& current_vertices);

    /**
     * Run the SAT model, adding colors until it is satisfiable.
     * @return true if satisfiable before hitting upper bound, false if
     * unsatisfiable or if we were aborted/interrupted.
     */
    bool p_run_until_satisfiable();

    /**
     * Reset the color classes to an empty state,
     * making sure that num_classes exist.
     * This is called when we replace the current coloring
     * with one from the SAT solver.
     */
    void p_reset_color_classes(std::size_t num_classes);

    /**
     * Look up a LitOrFixed value in a satisfying assignment of the SAT model.
     */
    template <typename ModelMap>
    bool p_lookup_in_model(const ModelMap& map, LitOrFixed l) const;

    /**
     * Completely reset the coloring in the rare case we failed to
     * extend all color classes to valid configurations.
     */
    void p_reset();

    /**
     * Create a LitOrFixed from a SatLit.
     */
    static inline LitOrFixed p_lit(SatLit sl) noexcept {
        return LitOrFixed(std::in_place_type<SatLit>, sl);
    }

    /**
     * Create a LitOrFixed from a fixed boolean value.
     */
    static inline LitOrFixed p_fixed(bool b) noexcept {
        return LitOrFixed(std::in_place_type<bool>, b);
    }

    /**
     * Called to initialize the clique model
     * on the first attempt to improve our clique.
     */
    bool
    p_first_run_clique_model(const std::vector<std::size_t>& current_vertices);

    /**
     * Called to add a complete assignment as
     * a <=1-constraint to the clique model.
     * The constraint RHS contains all vertices
     * that are present in the clique model and
     * covered by the complete assignment.
     * When new vertices are added to the model, they
     * are automatically added to all complete
     * assignment constraints that cover them.
     */
    void
    p_clq_add_complete_constraint(const DynamicBitset& complete_assignment);

    /**
     * Called to add a partial assignment as
     * a <=1-constraint to the clique model.
     * The constraint RHS contains all vertices
     * that are present in the clique model and
     * covered by the partial assignment.
     * When new vertices are added to the model, they
     * are automatically added to all partial
     * assignment constraints that cover them.
     * When the partial assignment is extended,
     * vertices present in the clique model
     * that are newly covered are added to the
     * constraint.
     */
    template <typename Propagator>
    void p_clq_add_partial_constraint(Propagator&& partial_assignment,
                                      bool expr_in_buffer = false);

    /**
     * Called to add a vertex as 0-1-variable to the clique model.
     */
    void p_clq_add_vertex(std::size_t vertex);

    /**
     * Notify that the partial assignment indicated by partial_index
     * has been extended (previous to the extension, the trail should
     * have had length old_trail_size).
     */
    void p_clq_extended_partial(std::size_t partial_index,
                                std::size_t old_trail_size);

    /**
     * Extend a partial assignment by adding a vertex.
     */
    bool p_clq_extend_partial_by_vertex(std::size_t partial_index,
                                        std::size_t vertex_index);

    /**
     * Setup model parameters.
     */
    void p_setup_clique_model();

    /**
     * Solve the current clique model relaxation.
     * @return true if the relaxation was solved, false on abort.
     */
    bool p_solve_relaxation();

    /**
     * Identify non-edges that violated the <=1-constraint
     * in the current clique model relaxation solution.
     * Actually checks for at least a violation by 0.01.
     * The violated nonedges are placed into m_clq_violated_nonedges.
     *
     * @return true if some violated non-edges were found, false otherwise.
     */
    bool p_identify_violated_nonedges();

    /**
     * Solve the 'full' relaxation, i.e.,
     * the relaxation including all lazily generated
     * non-edge constraints.
     * @return true if the full relaxation was solve, false on abort.
     */
    bool p_solve_full_relaxation();

    /**
     * Check if the indicated edge is actually a non-edge.
     * @return true if the edge is a non-edge, false otherwise.
     */
    bool p_is_nonedge(std::size_t vertex1, std::size_t vertex2);

    /**
     * Check if the indicated edge is actually an edge.
     * @return true if the edge is present in G_2, false otherwise.
     */
    bool p_is_edge(std::size_t vertex1, std::size_t vertex2) {
        return !p_is_nonedge(vertex1, vertex2);
    }

    /**
     * We lazily prohibit violated non-edges in the clique model.
     * The violated nonedges are stored in m_clq_violated_nonedges.
     * Each call to p_prohibit_violated_nonedges generates at most
     * one new constraint per vertex incident to a violated non-edge;
     * in practice, much fewer constraints are generated.
     * The method prefers strengthening existing partial constraints to
     * cover new violated non-edges over generating new constraints.
     */
    void p_prohibit_violated_nonedges();

    /**
     * Greedily extend the current partial assignment
     * in m_empty_propagator by adding non-edges from
     * m_clq_violated_nonedges[begin_index, end_index).
     */
    void p_extend_nonedge_cover_greedy(std::size_t begin_index,
                                       std::size_t end_index);

    /**
     * Try to handle the given violated non-edge by
     * extending an existing partial assignment constraint.
     */
    bool p_prohibit_extend_existing(OrderedSolutionIter vit,
                                    OrderedSolutionIter wit);

    /**
     * Use the clique model to try and improve the best known clique.
     * @return false if the search was aborted or the best clique was not
     * improved, true otherwise.
     */
    bool p_use_clique_model(const std::vector<std::size_t>& current_vertices);

    /**
     * Greedily round the current clique relaxation
     * to try and find an improved clique.
     */
    bool p_round_greedy();

    /**
     * Randomly turn the clique in m_vertex_set_buffer
     * with compatible vertices marked in m_vertex_bitset_buffer
     * into a maximal clique.
     */
    void p_randomly_extend_clique();

    /**
     * Ensure that the best clique found so far
     * is contained in the clique model.
     */
    void p_ensure_best_clique_in_model();

    /**
     * Try to add cuts to the clique model
     * coming from the current color classes.
     * @return true if cuts were added, false otherwise.
     */
    bool p_clique_model_cclass_cuts();

    /**
     * Run cheap cutting plane rounds.
     * @return false on abortion, true otherwise.
     */
    bool p_cheap_cut_rounds();

    /**
     * Run a single cheap cut round.
     * @return false on abortion or failure, true otherwise.
     */
    bool p_cheap_cut_round();

    /**
     * Run expensive cutting plane rounds.
     * @return false on abortion, true otherwise.
     */
    bool p_expensive_cut_rounds() {
        // TODO: implement expensive cuts (might not be very successful)
        return true;
    }

    /**
     * The gap in the current clique model.
     */
    double p_clique_gap() const noexcept {
        return (m_last_clique_objective - m_best_clique.size()) /
               m_best_clique.size();
    }

    /**
     * Greedily add vertices from the current relaxation
     * solution to existing partial-assignment cutting planes.
     * @return true if a constraint became violated, false otherwise.
     */
    bool p_greedy_add_to_cuts();

    /**
     * Greedily generate cutting planes from the current relaxation solution
     * and add them to the clique model.
     * @return true if a violated constraint was added, false otherwise.
     */
    bool p_greedy_generate_cuts();

    /**
     * Extend a greedy cut being built (partial assignment in
     * m_empty_propagator, vertices in m_vertex_set_buffer) by adding vertices
     * from the current relaxation solution.
     * @return the cut value of the extended cut.
     */
    double p_extend_greedy_cut(OrderedSolutionIter begin,
                               OrderedSolutionIter end);

    /**
     * Do a limited pricing round (limited means that there is a limit
     * on the number of vertices that can be added to the clique model).
     */
    bool p_limited_price_vertices(const std::vector<std::size_t>& vertices);

    /**
     * Do an unlimited pricing round (that adds all vertices
     * that could be an improvement to the clique model).
     */
    bool p_unlimited_price_vertices(const std::vector<std::size_t>& vertices);

    /**
     * Part of the pricing routine: Go through the given set of vertices,
     * and add good and possible vertices to the vectors m_pricing_good_vertices
     * and m_pricing_possible_vertices.
     * @return true if vertices were added, and false if there either were
     *         no possible vertices or if we were aborted.
     */
    bool
    p_price_collect_good_and_possible(const std::vector<std::size_t>& vertices);

    /**
     * For a vertex that is not present in the clique model,
     * compute the value of the LHS of its corresponding dual constraint.
     */
    double p_weight_in_dual(Vertex v);

    /**
     * Select a number of vertices from the set of
     * good vertices (or possible vertices if no good vertices exist).
     * The vertices are collected in m_vertex_set_buffer.
     */
    void p_pricing_select_with_goal(std::size_t goal_vertices);

    /**
     * Add the vertices in m_vertex_set_buffer to the clique model.
     */
    void p_pricing_add_vertices();

    /**
     * Report events if we have a recorder.
     */
    template <typename... PrintArgs>
    void p_report_event(std::string event_name, OutputObject data,
                        PrintArgs&&... args) {
        if (m_local_recorder) {
            m_local_recorder->store_event(std::move(event_name),
                                          std::move(data),
                                          std::forward<PrintArgs>(args)...);
        }
    }
};

template <typename IncrementalSATSolver>
CliqueSatDSaturSolver<IncrementalSATSolver>::CliqueSatDSaturSolver(
    std::vector<Vertex> considered_vertices,
    PairInfeasibilityMap* infeasibility_map, ClauseDB& clause_db,
    const std::vector<Vertex>& best_local_mes, std::size_t best_global_mes,
    std::size_t best_global_lb, std::vector<DynamicBitset> covering_assignments,
    bool enable_incremental_clique)
    : m_inf_map(infeasibility_map), m_n_all(clause_db.num_vars()),
      m_n_concrete(infeasibility_map->get_n_concrete()),
      m_empty_propagator(&clause_db),
      m_all_vertices(std::move(considered_vertices)),
      m_vertices_containing_concrete_literal(2 * m_n_concrete),
      m_clique_model_vertices_containing_concrete_literal(2 * m_n_concrete),
      m_vertices_implying_literal(2 * m_n_all,
                                  DynamicBitset(m_all_vertices.size())),
      m_adjacency_matrix(m_all_vertices.size(),
                         DynamicBitset(m_all_vertices.size())),
      m_definitive_nonedges(m_all_vertices.size(),
                            DynamicBitset(m_all_vertices.size())),
      m_incremental_clique_turned_off(!enable_incremental_clique),
      m_clique_model(sammy::gurobi_environment(true)),
      m_all_clique_vars(m_all_vertices.size()),
      m_present_in_clique_model(m_all_vertices.size()),
      m_in_currently_colored(m_all_vertices.size(), false),
      m_best_infos(20, CompareVertexIndexByInfo{this}),
      m_best_clique(p_clique_to_indices(best_local_mes)),
      m_last_sat_clique_size(0), m_initial_clique_size(best_local_mes.size()),
      m_initial_global_clique_size(best_global_mes),
      m_initial_global_lower_bound(best_global_lb),
      m_lower_bound(m_best_clique.size()),
      m_lower_bound_subgraph(m_best_clique),
      m_covering_assignments(std::move(covering_assignments)),
      m_greedy_cut_candidates(GREEDY_CUT_CANDIDATE_SET_SIZE) {
    p_initialize_graph();
    p_initialize_vertex_info();
    p_initialize_color_classes();
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_initialize_graph() {
    p_initialize_vertices_containing_concrete_literal();
    p_initialize_vertices_implying_literal();
    // TODO: the next function can take way too long.
    p_initialize_matrix_from_implied_literals();
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_initialize_vertices_containing_concrete_literal() {
    assert(m_vertices_containing_concrete_literal.size() == 2 * m_n_concrete);
    for (std::size_t vi = 0, vn = m_all_vertices.size(); vi < vn; ++vi) {
        Vertex v = m_all_vertices[vi];
        assert(lit::var(v.first) < m_n_concrete);
        assert(lit::var(v.second) < m_n_concrete);
        m_vertices_containing_concrete_literal[v.first].push_back(vi);
        m_vertices_containing_concrete_literal[v.second].push_back(vi);
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_initialize_vertices_implying_literal() {
    assert(m_vertices_implying_literal.size() == 2 * m_n_all);
    for (std::size_t vi = 0, n = m_all_vertices.size(); vi < n; ++vi) {
        assert(m_empty_propagator.get_current_level() == 0);
        Vertex v = m_all_vertices[vi];
        if (push_vertex(m_empty_propagator, v) < 0)
            throw std::logic_error("Infeasible interaction in m_all_vertices!");
        for (Lit l : m_empty_propagator.get_trail()) {
            m_vertices_implying_literal[l][vi].set();
        }
        m_empty_propagator.reset_to_zero();
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_initialize_matrix_from_implied_literals() {
    assert(m_adjacency_matrix.size() == m_all_vertices.size());
    assert(m_adjacency_matrix.empty() ||
           m_adjacency_matrix[0].size() == m_all_vertices.size());
    for (std::size_t vi = 0, n = m_all_vertices.size(); vi < n; ++vi) {
        assert(m_empty_propagator.get_current_level() == 0);
        Vertex v = m_all_vertices[vi];
        DynamicBitset& row = m_adjacency_matrix[vi];
        push_vertex(m_empty_propagator, v);
        for (Lit lpos : m_empty_propagator.get_trail()) {
            row |= m_vertices_implying_literal[lit::negate(lpos)];
        }
        m_empty_propagator.reset_to_zero();
        if (vi % 16 == 15)
            throw_if_interrupted();
    }
}

template <typename IncrementalSATSolver>
std::vector<std::size_t>
CliqueSatDSaturSolver<IncrementalSATSolver>::p_clique_to_indices(
    const std::vector<Vertex>& cvs, bool must_match) {
    PairHashSet<Vertex> in_clique(cvs.begin(), cvs.end());
    std::vector<std::size_t> result;
    result.reserve(cvs.size());
    for (std::size_t vi = 0, n = m_all_vertices.size(); vi < n; ++vi) {
        Vertex v = m_all_vertices[vi];
        if (in_clique.count(v)) {
            result.push_back(vi);
        }
    }
    if (must_match && (result.size() != cvs.size())) {
        throw std::logic_error("Clique contains missing vertex!");
    }
    return result;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::
    p_extract_vertices_from_configuration(const DynamicBitset& assignment,
                                          std::vector<std::size_t>& vertices) {
    Lit block = NIL;
    for (std::size_t i = 0, n = m_all_vertices.size(); i < n; ++i) {
        Vertex v = m_all_vertices[i];
        if (v.first == block)
            continue;
        if (!lit::is_true_in(v.first, assignment)) {
            block = v.first;
            continue;
        }
        if (lit::is_true_in(v.second, assignment)) {
            vertices.push_back(i);
        }
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::abort() {
    m_aborted.store(true);
    m_incremental_solver.terminate();
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::CompareVertexIndexByInfo::
operator()(std::size_t v1, std::size_t v2) const noexcept {
    // return true if 'v1' has higher priority than 'v2'
    const VertexInfo& vinfo1 = that->m_vertex_info[v1];
    const VertexInfo& vinfo2 = that->m_vertex_info[v2];
    if (vinfo1.in_some_class != vinfo2.in_some_class) {
        return !vinfo1.in_some_class;
    }
    if (vinfo1.num_open_classes != vinfo2.num_open_classes) {
        return vinfo1.num_open_classes < vinfo2.num_open_classes;
    }
    return vinfo1.degree > vinfo2.degree;
}

template <typename IncrementalSATSolver>
typename CliqueSatDSaturSolver<IncrementalSATSolver>::VertexChoice
CliqueSatDSaturSolver<IncrementalSATSolver>::p_choose_next_vertex() {
    VertexChoice result;
    // first priority: color leftover vertices from the clique
    if (p_choose_next_vertex_color_clique(result))
        return result;
    // second priority: color the vertices that only have one color left
    if (p_choose_next_vertex_one_color_left(result))
        return result;
    // third priority: use the DSatur-style choice
    if (!p_identify_best_candidates()) {
        result.vertex = std::numeric_limits<std::size_t>::max();
        result.type = VertexChoiceType::ALL_COVERED;
        return result;
    }
    std::size_t v = p_precise_best_candidates();
    result.vertex = v;
    if (!m_zero_candidates.empty()) {
        result.type = VertexChoiceType::WITH_ZERO_CHOICES;
    } else if (m_vertex_info[v].num_open_classes == 1) {
        result.type = VertexChoiceType::WITH_ONE_CHOICE;
    } else {
        result.type = VertexChoiceType::REGULAR_CHOICE;
    }
    return result;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_choose_next_vertex_color_clique(VertexChoice& vc) {
    while (!m_in_clique_uncolored.empty()) {
        std::size_t v = m_in_clique_uncolored.back();
        if (m_vertex_info[v].in_some_class) {
            m_in_clique_uncolored.pop_back();
            m_currently_colored.push_back(v);
            m_in_currently_colored[v].set();
        } else {
            vc = VertexChoice{v, VertexChoiceType::FROM_CLIQUE};
            return true;
        }
    }
    return false;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::
    p_choose_next_vertex_one_color_left(VertexChoice& vc) {
    while (!m_one_candidates.empty()) {
        std::size_t v = m_one_candidates.back();
        m_one_candidates.pop_back();
        if (!m_vertex_info[v].in_some_class) {
            vc = VertexChoice{v, VertexChoiceType::WITH_ONE_CHOICE};
            return true;
        }
    }
    return false;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_identify_best_candidates() {
    m_zero_candidates.clear();
    m_one_candidates.clear();
    m_best_infos.clear();
    for (std::size_t i = 0, n = m_all_vertices.size(); i < n; ++i) {
        if (m_vertex_info[i].in_some_class)
            continue;
        if (m_vertex_info[i].num_open_classes == 0) {
            m_zero_candidates.push_back(i);
        } else if (m_vertex_info[i].num_open_classes == 1) {
            m_one_candidates.push_back(i);
        }
        m_best_infos.push(i);
    }
    return !m_best_infos.elements().empty();
}

template <typename IncrementalSATSolver>
std::size_t
CliqueSatDSaturSolver<IncrementalSATSolver>::p_precise_best_candidates() {
    for (std::size_t candidate : m_best_infos.elements()) {
        Vertex vc = m_all_vertices[candidate];
        auto& vi = m_vertex_info[candidate];
        std::size_t open_before = vi.num_open_classes;
        if (open_before == 0)
            continue;
        for (std::size_t candidate_class : vi.open_classes.ones()) {
            auto& cclass = m_color_classes[candidate_class];
            SharedDBPropagator& prop = cclass.propagator;
            if (!can_push(prop, vc)) {
                vi.num_open_classes -= 1;
                vi.open_classes[candidate_class].reset();
            }
        }
        if (open_before != vi.num_open_classes) {
            if (vi.num_open_classes == 0) {
                m_zero_candidates.push_back(candidate);
            } else if (vi.num_open_classes == 1) {
                m_one_candidates.push_back(candidate);
            }
        }
    }
    const auto& elms = m_best_infos.elements();
    std::size_t best_v = *(std::min_element(elms.begin(), elms.end(),
                                            CompareVertexIndexByInfo{this}));
    return best_v;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_initialize_vertex_info() {
    const auto n = m_all_vertices.size();
    m_vertex_info.reserve(n);
    std::size_t ub_num_classes = m_covering_assignments.size();
    std::size_t initial_num_classes = m_best_clique.size();
    for (std::size_t vi = 0; vi < n; ++vi) {
        const DynamicBitset& row = m_adjacency_matrix[vi];
        std::size_t degree = row.count();
        m_vertex_info.emplace_back(initial_num_classes, ub_num_classes, degree);
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_initialize_color_classes() {
    if (m_empty_propagator.get_current_level() != 0)
        throw std::logic_error("Non-empty propagator!");
    m_color_classes.reserve(m_covering_assignments.size());
    for (std::size_t i = 0; i < m_best_clique.size(); ++i) {
        m_color_classes.emplace_back(m_empty_propagator, i);
        m_color_classes.back().add_vertex(this, m_best_clique[i]);
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::ColorClass::add_vertex(
    CliqueSatDSaturSolver* that, std::size_t vertex) {
    Vertex v = that->m_all_vertices[vertex];
    const auto& trail = propagator.get_trail();
    std::size_t old_trail_length = trail.size();
    int p = push_vertex(propagator, v);
    if (p < 0) {
        throw std::logic_error("Called add_vertex on incompatible class!");
    }
    if (p > 0) {
        p_added_vertex(that, old_trail_length);
    }
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::ColorClass::try_add_vertex(
    CliqueSatDSaturSolver* that, std::size_t vertex) {
    Vertex v = that->m_all_vertices[vertex];
    const auto& trail = propagator.get_trail();
    std::size_t old_trail_length = trail.size();
    int p = push_vertex(propagator, v);
    if (p < 0) {
        return false;
    }
    if (p > 0) {
        p_added_vertex(that, old_trail_length);
    }
    return true;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::ColorClass::p_added_vertex(
    CliqueSatDSaturSolver* that, std::size_t old_trail_end) {
    // TODO: unnecessary if we do the thing below (but that may be expensive)
    // except for maybe 'extra edges' (from extended graph G_2)?
    //    for(std::size_t vother : that->m_adjacency_matrix[vertex].ones()) {
    //        that->m_vertex_info[vother].close_class(index);
    //    }
    const auto& trail = propagator.get_trail();
    Lit nclit = that->m_n_concrete * 2;
    for (Lit lnew : IteratorRange{trail.begin() + old_trail_end, trail.end()}) {
        Lit lneg = lit::negate(lnew);
        if (lnew < nclit) {
            for (std::size_t vpot :
                 that->m_vertices_containing_concrete_literal[lnew])
            {
                Vertex vo = that->m_all_vertices[vpot];
                if (propagator.is_true(vo.first) &&
                    propagator.is_true(vo.second))
                {
                    that->m_vertex_info[vpot].in_some_class = true;
                }
            }
        }
        // TODO the thing below:
        for (std::size_t vout : that->m_vertices_implying_literal[lneg].ones())
        {
            that->m_vertex_info[vout].close_class(index);
        }
    }
}

template <typename IncrementalSATSolver>
auto CliqueSatDSaturSolver<IncrementalSATSolver>::solve() -> SolveResult {
    for (;;) {
        if (m_aborted.load()) {
            return SolveResult::ABORTED;
        }
        VertexChoice vc = p_choose_next_vertex();
        switch (vc.type) {
        case VertexChoiceType::ALL_COVERED:
            if (!p_complete_all_classes()) {
                if (m_aborted.load())
                    return SolveResult::ABORTED;
                if (m_lower_bound == m_covering_assignments.size())
                    return SolveResult::SOLUTION_WAS_OPTIMAL;
                continue;
            }
            if (m_covering_assignments.size() <= m_color_classes.size())
                return SolveResult::SOLUTION_WAS_OPTIMAL;
            return SolveResult::IMPROVED_SOLUTION;

        case VertexChoiceType::WITH_ZERO_CHOICES: {
            p_report_event("BEGIN_HANDLE_ZERO_CHOICE",
                           {{"vertices_handled", m_currently_colored.size()}},
                           "vertices_handled");
            bool handled = p_handle_zero_choice(vc.vertex);
            p_report_event("DONE_HANDLE_ZERO_CHOICE",
                           {{"vertices_handled", m_currently_colored.size()}},
                           "vertices_handled");
            if (handled) {
                p_colored_vertex(vc.vertex);
                continue;
            }
            return m_aborted.load() ? SolveResult::ABORTED
                                    : SolveResult::SOLUTION_WAS_OPTIMAL;
        }

        default: {
            if (!p_handle_vertex(vc.vertex)) {
                p_report_event(
                    "BEGIN_HANDLE_ZERO_CHOICE",
                    {{"vertices_handled", m_currently_colored.size()}},
                    "vertices_handled");
                bool handled = p_handle_zero_choice(vc.vertex);
                p_report_event(
                    "DONE_HANDLE_ZERO_CHOICE",
                    {{"vertices_handled", m_currently_colored.size()}},
                    "vertices_handled");
                if (!handled) {
                    return m_aborted.load() ? SolveResult::ABORTED
                                            : SolveResult::SOLUTION_WAS_OPTIMAL;
                }
            }
            p_colored_vertex(vc.vertex);
            break;
        }
        }
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_colored_vertex(
    std::size_t vertex) {
    if (!m_in_clique_uncolored.empty() &&
        m_in_clique_uncolored.back() == vertex)
    {
        m_in_clique_uncolored.pop_back();
    }
    if (!m_one_candidates.empty() && m_one_candidates.back() == vertex) {
        m_one_candidates.pop_back();
    }
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_handle_vertex(
    std::size_t vertex) {
    auto& vinfo = m_vertex_info[vertex];
    if (vinfo.in_some_class) {
        m_currently_colored.push_back(vertex);
        m_in_currently_colored[vertex].set();
        return true;
    }
    for (std::size_t candidate : vinfo.open_classes.ones()) {
        auto& cclass = m_color_classes[candidate];
        if (cclass.try_add_vertex(this, vertex)) {
            m_currently_colored.push_back(vertex);
            m_in_currently_colored[vertex].set();
            return true;
        } else {
            vinfo.close_class(candidate);
        }
    }
    return false;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_handle_zero_choice(
    std::size_t vertex) {
    auto current_vertices = p_compute_current_vertex_set(vertex);
    bool rclique = false;
    if (!m_incremental_clique_turned_off) {
        p_report_event("BEGIN_ZERO_CHOICE_CLIQUE", {});
        rclique = p_zero_choice_run_clique(current_vertices);
        p_report_event("DONE_ZERO_CHOICE_CLIQUE",
                       {{"best_clique_size", m_best_clique.size()}},
                       "best_clique_size");
    }
    if (rclique) {
        // found a better clique (m_best_clique, m_lower_bound updated)!
        if (m_best_clique.size() > m_color_classes.size()) {
            if (m_best_clique.size() == m_covering_assignments.size()) {
                return false;
            }
            p_extend_sorted_unique_set(current_vertices, m_best_clique);
            p_reset_sat_model();
            p_rebuild_sat_model(current_vertices);
            return p_run_until_satisfiable();
        }
    }
    if (m_aborted.load())
        return false;
    p_report_event("BEGIN_ZERO_CHOICE_SAT_MODELING", {});
    if (m_best_clique.size() != m_last_sat_clique_size) {
        p_reset_sat_model();
        p_rebuild_sat_model(current_vertices);
    } else {
        p_extend_sat_model_to_classes();
        p_extend_sat_model(current_vertices);
    }
    p_report_event("DONE_ZERO_CHOICE_SAT_MODELING", {});
    p_report_event("BEGIN_ZERO_CHOICE_SAT_SOLVE", {});
    std::optional<bool> result = p_run_sat_model();
    p_report_event("DONE_ZERO_CHOICE_SAT_SOLVE", {});
    if (!result)
        return false;
    if (*result)
        return true;
    m_lower_bound += 1;
    m_lower_bound_subgraph = std::move(current_vertices);
    p_call_lb_callback();
    if (m_lower_bound == m_covering_assignments.size())
        return false;
    p_add_color_class(vertex);
    p_compute_in_clique_uncolored();
    return true;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_reset_sat_model() {
    if (m_class_literals.empty())
        return;
    m_incremental_solver.reset();
    m_class_literals.clear();
    m_vertex_in_class.clear();
    for (std::size_t vi : m_sat_vertex_order) {
        m_vertex_info[vi].sat_index = std::numeric_limits<std::size_t>::max();
    }
    m_sat_vertex_order.clear();
    m_not_enough_colors.reset();
    m_not_enough_colors_is_for = 0;
}

template <typename IncrementalSATSolver>
std::vector<std::size_t>
CliqueSatDSaturSolver<IncrementalSATSolver>::p_compute_current_vertex_set(
    std::size_t new_vertex) {
    std::vector<std::size_t> cover_vertex_set = m_in_clique_uncolored;
    cover_vertex_set.insert(cover_vertex_set.end(), m_zero_candidates.begin(),
                            m_zero_candidates.end());
    cover_vertex_set.push_back(new_vertex);
    cover_vertex_set.insert(cover_vertex_set.end(), m_currently_colored.begin(),
                            m_currently_colored.end());
    std::sort(cover_vertex_set.begin(), cover_vertex_set.end());
    cover_vertex_set.erase(
        std::unique(cover_vertex_set.begin(), cover_vertex_set.end()),
        cover_vertex_set.end());
    return cover_vertex_set;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_add_color_class(
    std::size_t vertex) {
    std::size_t new_class = m_color_classes.size();
    for (std::size_t vi = 0, n = m_all_vertices.size(); vi < n; ++vi) {
        m_vertex_info[vi].open_classes.push_back(true);
        m_vertex_info[vi].num_open_classes += 1;
    }
    m_color_classes.emplace_back(m_empty_propagator, new_class);
    m_color_classes.back().add_vertex(this, vertex);
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_compute_in_clique_uncolored() {
    m_in_clique_uncolored.clear();
    for (std::size_t v : m_best_clique) {
        if (!m_in_currently_colored[v]) {
            m_in_clique_uncolored.push_back(v);
        }
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_extend_sorted_unique_set(
    std::vector<std::size_t>& sorted_unique,
    const std::vector<std::size_t>& to_add) {
    auto ins_begin =
        sorted_unique.insert(sorted_unique.end(), to_add.begin(), to_add.end());
    std::sort(ins_begin, sorted_unique.end());
    std::inplace_merge(sorted_unique.begin(), ins_begin, sorted_unique.end());
    sorted_unique.erase(std::unique(sorted_unique.begin(), sorted_unique.end()),
                        sorted_unique.end());
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_rebuild_sat_model(
    const std::vector<std::size_t>& current_vertices) {
    // on calling this, the model should be reset/empty
    assert(m_sat_vertex_order.empty());
    assert(m_vertex_in_class.empty());
    assert(m_class_literals.empty());

    // rebuild SAT order
    m_sat_vertex_order.reserve(current_vertices.size());
    for (std::size_t vi : current_vertices) {
        m_vertex_info[vi].sat_index = m_sat_vertex_order.size();
        m_sat_vertex_order.push_back(vi);
    }

    // rebuild classes
    m_last_sat_clique_size = m_best_clique.size();
    std::size_t goal_classes = m_lower_bound;
    for (std::size_t ci = 0; ci < goal_classes; ++ci) {
        std::optional<std::size_t> clique_vertex = std::nullopt;
        if (ci < m_best_clique.size()) {
            clique_vertex = m_best_clique[ci];
        }
        p_rebuild_sat_color_class(clique_vertex);
    }

    // rebuild vertex constraints
    SatLit not_enough = m_incremental_solver.new_var();
    m_not_enough_colors = not_enough;
    m_not_enough_colors_is_for = goal_classes;
    for (std::size_t v : range(m_sat_vertex_order.size())) {
        p_add_vertex_color_constraint(v, not_enough);
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_rebuild_sat_color_class(
    std::optional<std::size_t> clique_vertex) {
    if (clique_vertex) {
        p_append_constrained_class_lits(*clique_vertex);
    } else {
        p_append_unconstrained_class_lits();
    }
    p_append_vertex_in_class();
    p_enforce_formula_on(m_class_literals.size() - 1);
}

template <typename IncrementalSATSolver>
auto CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_append_unconstrained_class_lits()
    -> std::vector<LitOrFixed>& {
    m_class_literals.emplace_back();
    auto& cl = m_class_literals.back();
    cl.reserve(m_n_all);
    for (Var v = 0; v < m_n_all; ++v) {
        Lit l = lit::positive_lit(v);
        if (m_empty_propagator.is_open(l)) {
            cl.emplace_back(p_lit(m_incremental_solver.new_var()));
        } else {
            cl.emplace_back(p_fixed(m_empty_propagator.is_true(l)));
        }
    }
    return cl;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::
    p_append_constrained_class_lits(std::size_t clique_vertex) {
    Vertex v = m_all_vertices[clique_vertex];
    push_vertex(m_empty_propagator, v);
    m_class_literals.emplace_back();
    auto& cl = m_class_literals.back();
    cl.reserve(m_n_all);
    for (Var v = 0; v < m_n_all; ++v) {
        Lit l = lit::positive_lit(v);
        if (m_empty_propagator.is_open(l)) {
            cl.emplace_back(m_incremental_solver.new_var());
        } else {
            cl.emplace_back(m_empty_propagator.is_true(l));
        }
    }
    m_empty_propagator.reset_to_zero();
}

template <typename IncrementalSATSolver>
auto CliqueSatDSaturSolver<IncrementalSATSolver>::p_append_vertex_in_class()
    -> std::vector<LitOrFixed>& {
    std::size_t class_index = m_vertex_in_class.size();
    m_vertex_in_class.emplace_back();
    assert(m_class_literals.size() == m_vertex_in_class.size());
    auto& v_in_c = m_vertex_in_class.back();
    v_in_c.reserve(m_sat_vertex_order.size());
    for (std::size_t vi : m_sat_vertex_order) {
        auto [l1, l2] = m_all_vertices[vi];
        LitOrFixed cl1 = p_lookup_sat_literal_in_class(l1, class_index);
        LitOrFixed cl2 = p_lookup_sat_literal_in_class(l2, class_index);
        v_in_c.push_back(p_and(cl1, cl2));
    }
    return v_in_c;
}

template <typename IncrementalSATSolver>
auto CliqueSatDSaturSolver<IncrementalSATSolver>::p_lookup_sat_literal_in_class(
    Lit l, std::size_t class_index) -> LitOrFixed {
    Var v = lit::var(l);
    bool neg = lit::negative(l);
    LitOrFixed v_in_c = m_class_literals[class_index][v];
    return neg ? p_negate(v_in_c) : v_in_c;
}

template <typename IncrementalSATSolver>
auto CliqueSatDSaturSolver<IncrementalSATSolver>::p_negate(LitOrFixed l)
    -> LitOrFixed {
    return std::visit(overloaded{[](bool& b) { return p_fixed(!b); },
                                 [](SatLit& l) { return p_lit(-l); }},
                      l);
}

template <typename IncrementalSATSolver>
auto CliqueSatDSaturSolver<IncrementalSATSolver>::p_and(LitOrFixed l1,
                                                        LitOrFixed l2)
    -> LitOrFixed {
    auto one_fixed = [&](SatLit b1, bool b2) {
        if (!b2)
            return p_fixed(false);
        return p_lit(b1);
    };
    return std::visit(
        overloaded{[&](bool& b1, bool& b2) { return p_fixed(b1 & b2); },
                   [&](bool& b1, SatLit& b2) { return one_fixed(b2, b1); },
                   [&](SatLit& b1, bool& b2) { return one_fixed(b1, b2); },
                   [&](SatLit& b1, SatLit& b2) {
                       SatLit nv = m_incremental_solver.new_var();
                       m_incremental_solver.add_short_clause(-nv, b1);
                       m_incremental_solver.add_short_clause(-nv, b2);
                       m_incremental_solver.add_short_clause(nv, -b1, -b2);
                       return p_lit(nv);
                   }},
        l1, l2);
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_enforce_formula_on(
    std::size_t class_index) {
    const auto& db = m_empty_propagator.db();
    for (auto [l1, l2] : db.binary_clauses()) {
        const Lit l[2] = {l1, l2};
        p_add_clause(l + 0, l + 2, class_index);
    }
    for (CRef c = 1, cn = db.literal_db_size(); c < cn; c = db.next_clause(c)) {
        auto r = db.lits_of(c);
        p_add_clause(r.begin(), r.end(), class_index);
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_add_clause(
    const Lit* begin, const Lit* end, std::size_t class_index) {
    m_clause_buffer.clear();
    for (Lit l : IteratorRange{begin, end}) {
        LitOrFixed lf = p_lookup_sat_literal_in_class(l, class_index);
        bool clause_satisfied =
            std::visit(overloaded{[&](bool& b) -> bool { return b; },
                                  [&](SatLit& s) -> bool {
                                      m_clause_buffer.push_back(s);
                                      return false;
                                  }},
                       lf);
        if (clause_satisfied)
            return;
    }
    if (m_clause_buffer.empty())
        throw std::logic_error("UNSAT class!");
    m_incremental_solver.add_clause(m_clause_buffer.begin(),
                                    m_clause_buffer.end());
}

template <typename IncrementalSATSolver>
std::optional<bool>
CliqueSatDSaturSolver<IncrementalSATSolver>::p_run_sat_model() {
    std::vector<SatLit> assumptions;
    assumptions.push_back(-*m_not_enough_colors);
    std::optional<bool> result = m_incremental_solver.solve(assumptions);
    if (!result || !*result) {
        return result;
    }
    p_extract_sat_model(m_incremental_solver.get_model());
    return true;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_add_vertex_color_constraint(
    std::size_t sat_index, SatLit not_enough) {
    m_clause_buffer.clear();
    for (std::size_t c = 0; c < m_vertex_in_class.size(); ++c) {
        LitOrFixed v_in_c = m_vertex_in_class[c][sat_index];
        bool is_satisfied =
            std::visit(overloaded{[&](bool& b) -> bool { return b; },
                                  [&](SatLit& l) -> bool {
                                      m_clause_buffer.push_back(l);
                                      return false;
                                  }},
                       v_in_c);
        if (is_satisfied)
            return;
    }
    m_clause_buffer.push_back(not_enough);
    m_incremental_solver.add_clause(m_clause_buffer.begin(),
                                    m_clause_buffer.end());
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_extend_sat_model(
    const std::vector<std::size_t>& current_vertices) {
    std::vector<std::size_t> new_vertices;
    std::size_t old_sat_count = m_sat_vertex_order.size();
    for (std::size_t v : current_vertices) {
        if (m_vertex_info[v].sat_index ==
            std::numeric_limits<std::size_t>::max())
        {
            std::size_t sindex = m_sat_vertex_order.size();
            m_sat_vertex_order.push_back(v);
            m_vertex_info[v].sat_index = sindex;
            new_vertices.push_back(v);
        }
    }
    std::size_t num_classes = m_vertex_in_class.size();
    for (std::size_t c = 0; c < num_classes; ++c) {
        auto& v_in_c = m_vertex_in_class[c];
        for (std::size_t new_vertex : new_vertices) {
            auto [l1, l2] = m_all_vertices[new_vertex];
            LitOrFixed sl1 = p_lookup_sat_literal_in_class(l1, c);
            LitOrFixed sl2 = p_lookup_sat_literal_in_class(l2, c);
            v_in_c.push_back(p_and(sl1, sl2));
        }
    }
    for (std::size_t new_vertex :
         range(old_sat_count, old_sat_count + new_vertices.size()))
    {
        p_add_vertex_color_constraint(new_vertex, *m_not_enough_colors);
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_extend_sat_model_to_classes() {
    std::size_t old_classes = m_vertex_in_class.size();
    if (old_classes == m_lower_bound)
        return;
    for (std::size_t c = m_vertex_in_class.size(); c < m_lower_bound; ++c) {
        std::optional<std::size_t> clique_vertex = std::nullopt;
        if (c < m_best_clique.size()) {
            clique_vertex = m_best_clique[c];
        }
        p_rebuild_sat_color_class(clique_vertex);
    }
    p_check_not_enough_colors();
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_check_not_enough_colors() {
    if (!m_not_enough_colors || m_not_enough_colors_is_for != m_lower_bound) {
        if (m_not_enough_colors)
            m_incremental_solver.fix(*m_not_enough_colors);
        m_not_enough_colors = m_incremental_solver.new_var();
        m_not_enough_colors_is_for = m_lower_bound;
        for (std::size_t v : range(m_sat_vertex_order.size())) {
            p_add_vertex_color_constraint(v, *m_not_enough_colors);
        }
    }
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_zero_choice_run_clique(
    const std::vector<std::size_t>& current_vertices) {
    std::size_t old_clique_size = m_best_clique.size();
    if (m_clique_candidate_callback) {
        auto candidate = m_clique_candidate_callback();
        if (candidate.size() > old_clique_size) {
            auto candidate_indices = p_clique_to_indices(candidate, false);
            if (candidate_indices.size() > old_clique_size) {
                m_best_clique = std::move(candidate_indices);
                if (m_best_clique.size() > m_lower_bound) {
                    m_lower_bound = m_best_clique.size();
                    m_lower_bound_subgraph = m_best_clique;
                    p_call_lb_callback();
                }
                if (!m_existing_clique_vars.empty()) {
                    p_ensure_best_clique_in_model();
                }
            }
        }
    }
    if (m_existing_clique_vars.empty()) {
        if (!p_first_run_clique_model(current_vertices))
            return false;
    } else {
        if (m_existing_clique_var_value_count !=
            m_existing_clique_var_indices.size())
        {
            if (!p_solve_full_relaxation())
                return false;
        }
        if (p_clique_model_cclass_cuts()) {
            if (!p_solve_full_relaxation())
                return false;
        }
    }
    if (!p_use_clique_model(current_vertices))
        return false;
    return old_clique_size < m_best_clique.size();
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_run_until_satisfiable() {
    for (;;) {
        auto res = p_run_sat_model();
        if (!res)
            return false;
        if (!*res) {
            m_lower_bound += 1;
            m_lower_bound_subgraph = m_sat_vertex_order;
            p_call_lb_callback();
            if (m_lower_bound == m_covering_assignments.size()) {
                return false;
            }
            p_extend_sat_model_to_classes();
        } else {
            return true;
        }
    }
}

template <typename IncrementalSATSolver>
template <typename ModelMap>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_extract_sat_model(
    const ModelMap& model_map) {
    const std::size_t sn = m_sat_vertex_order.size();
    const std::size_t cn = m_vertex_in_class.size();
    p_reset_color_classes(cn);
    std::vector<std::size_t> class_with_vertex;
    for (std::size_t si = 0; si < sn; ++si) {
        for (std::size_t c = 0; c < cn; ++c) {
            LitOrFixed l = m_vertex_in_class[c][si];
            if (p_lookup_in_model(model_map, l)) {
                class_with_vertex.push_back(c);
                break;
            }
        }
    }
    if (class_with_vertex.size() != sn) {
        std::cerr << "SAT model does not assign a color to each vertex!"
                  << std::endl;
        std::cerr << "SAT solver name: " << m_incremental_solver.name()
                  << std::endl;
        throw std::logic_error(
            "SAT model does not assign a color to each vertex!");
    }
    for (std::size_t si = 0; si < sn; ++si) {
        std::size_t vi = m_sat_vertex_order[si];
        m_color_classes[class_with_vertex[si]].add_vertex(this, vi);
        if (!m_in_currently_colored[vi]) {
            m_in_currently_colored[vi].set();
            m_currently_colored.push_back(vi);
        }
    }
}

template <typename IncrementalSATSolver>
template <typename ModelMap>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_lookup_in_model(
    const ModelMap& map, LitOrFixed l) const {
    return std::visit(overloaded{[&](bool& b) -> bool { return b; },
                                 [&](SatLit& sl) -> bool { return map[sl]; }},
                      l);
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_reset_color_classes(
    std::size_t num_classes) {
    const std::size_t cold = m_color_classes.size();
    if (num_classes > cold) {
        for (std::size_t c = 0; c < cold; ++c) {
            m_color_classes[c].propagator.reset_to_zero();
        }
        for (std::size_t c = cold; c < num_classes; ++c) {
            m_color_classes.emplace_back(m_empty_propagator, c);
        }
        for (VertexInfo& v : m_vertex_info) {
            v.open_classes.resize(num_classes);
            v.open_classes.set();
            v.num_open_classes = num_classes;
            v.in_some_class = false;
        }
    } else {
        for (std::size_t c = 0; c < num_classes; ++c) {
            m_color_classes[c].propagator.reset_to_zero();
        }
        for (VertexInfo& v : m_vertex_info) {
            v.open_classes.set();
            v.num_open_classes = num_classes;
            v.in_some_class = false;
        }
    }
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_complete_all_classes() {
    struct EmptyHandler {};
    EmptyHandler handler;
    ClassCompleter<EmptyHandler> completer{m_n_concrete, m_n_all, &handler};
    bool result = true;
    for (std::size_t cindex : range(m_color_classes.size())) {
        SharedDBPropagator& prop = m_color_classes[cindex].propagator;
        if (!completer.complete_class(prop)) {
            result = false;
        }
    }
    if (result)
        return true;
    p_reset();
    return false;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_reset() {
    for (VertexInfo& v : m_vertex_info) {
        v.open_classes.set();
        v.num_open_classes = m_color_classes.size();
        v.in_some_class = false;
    }
    p_reset_sat_model();
    m_empty_propagator.incorporate_or_throw();
    for (ColorClass& cc : m_color_classes) {
        cc.propagator.reset_or_throw();
    }
    m_currently_colored = m_best_clique;
    m_in_currently_colored.reset();
    m_zero_candidates.clear();
    m_one_candidates.clear();
    m_in_clique_uncolored.clear();
    for (std::size_t vi : m_currently_colored) {
        m_in_currently_colored[vi].set();
    }
    for (std::size_t ci : range(m_best_clique.size())) {
        std::size_t vi = m_best_clique[ci];
        m_color_classes[ci].add_vertex(this, vi);
    }
}

template <typename IncrementalSATSolver>
PartialSolution
CliqueSatDSaturSolver<IncrementalSATSolver>::get_partial_solution() const {
    auto transformer =
        [](const ColorClass& color) -> const SharedDBPropagator& {
        return color.propagator;
    };
    auto begin = boost::iterators::make_transform_iterator(
        m_color_classes.begin(), +transformer);
    auto end = boost::iterators::make_transform_iterator(m_color_classes.end(),
                                                         +transformer);
    return PartialSolution(m_n_all, m_inf_map, begin, end);
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_first_run_clique_model(
    const std::vector<std::size_t>& current_vertices) {
    p_setup_clique_model();
    for (const DynamicBitset& ca : m_covering_assignments) {
        p_clq_add_complete_constraint(ca);
    }
    for (const ColorClass& cc : m_color_classes) {
        p_clq_add_partial_constraint(cc.propagator);
    }
    for (std::size_t vertex : m_best_clique) {
        m_last_vertex_addition_source = "first_run_clique_model 1";
        p_clq_add_vertex(vertex);
    }
    auto added_vertices =
        sample_from_range(current_vertices, 1000, sammy::rng());
    for (std::size_t vertex : added_vertices) {
        m_last_vertex_addition_source = "first_run_clique_model 2";
        p_clq_add_vertex(vertex);
    }
    if (!p_solve_full_relaxation()) {
        return false;
    }
    return true;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_clq_add_complete_constraint(
    const DynamicBitset& complete_assignment) {
    m_buffer_expr.clear();
    for (std::size_t i = 0, ncv = m_existing_clique_vars.size(); i < ncv; ++i) {
        std::size_t vi = m_existing_clique_var_indices[i];
        Vertex v = m_all_vertices[vi];
        if (lit::is_true_in(v.first, complete_assignment) &&
            lit::is_true_in(v.second, complete_assignment))
        {
            m_buffer_expr += m_existing_clique_vars[i];
        }
    }
    m_clq_constr_complete_assignments.push_back(complete_assignment);
    m_clq_constr_complete_constraints.push_back(
        m_clique_model.addConstr(m_buffer_expr <= 1));
}

template <typename IncrementalSATSolver>
template <typename Propagator>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_clq_add_partial_constraint(
    Propagator&& partial_assignment, bool expr_in_buffer) {
    if (!expr_in_buffer) {
        m_buffer_expr.clear();
        for (std::size_t i = 0, ncv = m_existing_clique_vars.size(); i < ncv;
             ++i)
        {
            std::size_t vi = m_existing_clique_var_indices[i];
            Vertex v = m_all_vertices[vi];
            if (partial_assignment.is_true(v.first) &&
                partial_assignment.is_true(v.second))
            {
                m_buffer_expr += m_existing_clique_vars[i];
            }
        }
    }
    m_clq_constr_partial_assignments.emplace_back(
        std::forward<Propagator>(partial_assignment));
    m_clq_constr_partial_constraints.push_back(
        m_clique_model.addConstr(m_buffer_expr <= 1));
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_clq_add_vertex(
    std::size_t vertex) {
    if (m_present_in_clique_model[vertex])
        return;
    m_column_buffer.clear();
    Vertex v = m_all_vertices[vertex];
    std::size_t clique_model_index = m_existing_clique_vars.size();
    for (std::size_t complete_i :
         range(m_clq_constr_complete_assignments.size()))
    {
        const auto& assignment = m_clq_constr_complete_assignments[complete_i];
        if (lit::is_true_in(v.first, assignment) &&
            lit::is_true_in(v.second, assignment))
        {
            m_column_buffer.addTerm(
                1.0, m_clq_constr_complete_constraints[complete_i]);
        }
    }
    for (std::size_t partial_i : range(m_clq_constr_partial_assignments.size()))
    {
        const auto& assignment = m_clq_constr_partial_assignments[partial_i];
        if (assignment.is_true(v.first) && assignment.is_true(v.second)) {
            m_column_buffer.addTerm(
                1.0, m_clq_constr_partial_constraints[partial_i]);
        }
    }
    GRBVar var =
        m_clique_model.addVar(0.0, 1.0, 1.0, GRB_CONTINUOUS, m_column_buffer);
    m_present_in_clique_model[vertex].set();
    m_existing_clique_var_indices.push_back(vertex);
    m_existing_clique_vars.push_back(var);
    m_all_clique_vars[vertex] = var;
    m_clique_model_vertices_containing_concrete_literal[v.first].push_back(
        clique_model_index);
    m_clique_model_vertices_containing_concrete_literal[v.second].push_back(
        clique_model_index);
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_clq_extended_partial(
    std::size_t partial_index, std::size_t old_trail_size) {
    const auto& propagator = m_clq_constr_partial_assignments[partial_index];
    const auto& trail = propagator.get_trail();
    const Lit nclit = 2 * m_n_concrete;
    GRBConstr constr = m_clq_constr_partial_constraints[partial_index];
    for (Lit lnew : IteratorRange{trail.begin() + old_trail_size, trail.end()})
    {
        if (lnew >= nclit)
            continue;
        for (std::size_t clique_model_index :
             m_clique_model_vertices_containing_concrete_literal[lnew])
        {
            std::size_t vertex_index =
                m_existing_clique_var_indices[clique_model_index];
            Vertex v = m_all_vertices[vertex_index];
            if (propagator.is_true(v.first) && propagator.is_true(v.second)) {
                m_clique_model.chgCoeff(
                    constr, m_existing_clique_vars[clique_model_index], 1.0);
            }
        }
    }
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::
    p_clq_extend_partial_by_vertex(std::size_t partial_index,
                                   std::size_t vertex) {
    Vertex v = m_all_vertices[vertex];
    auto& propagator = m_clq_constr_partial_assignments[partial_index];
    const std::size_t old_trail_size = propagator.get_trail().size();
    int p = push_vertex(propagator, v);
    if (p < 0)
        return false;
    if (p == 0)
        return true;
    p_clq_extended_partial(partial_index, old_trail_size);
    return true;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<
    IncrementalSATSolver>::GurobiCallbackCheckAbort::callback() {
    if (that->m_aborted.load())
        this->abort();
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_setup_clique_model() {
    // thread count, sense and experimentally determined method/presolve level
    m_clique_model.set(GRB_IntParam_Threads, 1);
    m_clique_model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    m_clique_model.set(GRB_IntParam_Method, GRB_METHOD_DUAL);
    m_clique_model.set(GRB_IntParam_Presolve, 1);
    m_check_abort_cb = std::make_unique<GurobiCallbackCheckAbort>(this);
    m_clique_model.setCallback(m_check_abort_cb.get());
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_solve_relaxation() {
    bool repeat = true;
    while (repeat) {
        repeat = false;
        m_clique_model.optimize();
        int status = m_clique_model.get(GRB_IntAttr_Status);
        if (status == GRB_INTERRUPTED) {
            m_aborted.store(true);
            return false;
        }
        if (status != GRB_OPTIMAL) {
            throw std::logic_error(
                "Unexpected status after solving relaxation!");
        }
        m_existing_clique_var_value_count = m_existing_clique_vars.size();
        m_last_vertex_addition_source = "cleared after solving";
        m_existing_clique_var_values.reset(
            m_clique_model.get(GRB_DoubleAttr_X, m_existing_clique_vars.data(),
                               int(m_existing_clique_vars.size())));
        m_clq_constr_complete_dual_values.reset(m_clique_model.get(
            GRB_DoubleAttr_Pi, m_clq_constr_complete_constraints.data(),
            int(m_clq_constr_complete_constraints.size())));
        m_clq_constr_partial_dual_values.reset(m_clique_model.get(
            GRB_DoubleAttr_Pi, m_clq_constr_partial_constraints.data(),
            int(m_clq_constr_partial_constraints.size())));
        m_last_clique_objective = m_clique_model.get(GRB_DoubleAttr_ObjVal);
        m_last_rounded_objective = static_cast<std::size_t>(
            std::floor(m_last_clique_objective + 0.001));
        m_clq_ordered_solution.clear();
        for (std::size_t i = 0; i < m_existing_clique_vars.size(); ++i) {
            assert(i < m_existing_clique_var_value_count);
            double x = m_existing_clique_var_values[i];
            if (x >= 1.0e-4) {
                m_clq_ordered_solution.emplace_back(OrderedSolutionValue{i, x});
            }
        }
        std::sort(m_clq_ordered_solution.begin(), m_clq_ordered_solution.end());
        if (m_aborted.load())
            return false;
        if (p_round_greedy()) {
            repeat = true;
        }
    }
    return !m_aborted.load();
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_solve_full_relaxation() {
    for (;;) {
        if (!p_solve_relaxation())
            return false;
        if (m_last_rounded_objective <= m_best_clique.size())
            return true;
        if (!p_identify_violated_nonedges())
            return true;
        p_prohibit_violated_nonedges();
    }
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_identify_violated_nonedges() {
    m_clq_violated_nonedges.clear();
    for (OrderedSolutionIter i = m_clq_ordered_solution.begin(),
                             e = m_clq_ordered_solution.end();
         i != e; ++i)
    {
        double weight = i->value;
        double thresh = 1.01 - weight;
        if (weight < thresh)
            break;
        std::size_t clique_index = i->clique_model_index;
        std::size_t vertex_index = m_existing_clique_var_indices[clique_index];
        for (OrderedSolutionIter j = std::next(i); j != e; ++j) {
            double w_j = j->value;
            if (w_j < thresh)
                break;
            std::size_t vertex_index2 =
                m_existing_clique_var_indices[j->clique_model_index];
            if (p_is_nonedge(vertex_index, vertex_index2)) {
                m_clq_violated_nonedges.emplace_back(i, j);
            }
        }
    }
    return !m_clq_violated_nonedges.empty();
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_is_nonedge(
    std::size_t vertex1, std::size_t vertex2) {
    if (m_adjacency_matrix[vertex1][vertex2])
        return false;
    if (m_definitive_nonedges[vertex1][vertex2])
        return true;
    if (!push_vertex_pair(m_empty_propagator, m_all_vertices[vertex1],
                          m_all_vertices[vertex2]))
    {
        m_adjacency_matrix[vertex1][vertex2] = true;
        m_adjacency_matrix[vertex2][vertex1] = true;
        return false;
    } else {
        m_empty_propagator.reset_to_zero();
        m_definitive_nonedges[vertex1][vertex2] = true;
        m_definitive_nonedges[vertex2][vertex1] = true;
        return true;
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_prohibit_violated_nonedges() {
    m_covered_violated_nonedges.assign(m_clq_violated_nonedges.size(), false);
    auto last_first = m_clq_ordered_solution.cend();
    const auto nve = m_clq_violated_nonedges.size();
    std::size_t num_new_constraints = 0;
    for (std::size_t i = 0; i < nve; ++i) {
        auto [vit, wit] = m_clq_violated_nonedges[i];
        if (vit == last_first || m_covered_violated_nonedges[i])
            continue;
        if (p_prohibit_extend_existing(vit, wit)) {
            m_covered_violated_nonedges[i].set();
            continue;
        }
        last_first = vit;
        if (++num_new_constraints > CAP_PROHIBIT_VIOLATED_NUM_NEW_CONSTRAINTS)
            continue;
        std::size_t vind =
            m_existing_clique_var_indices[vit->clique_model_index];
        std::size_t wind =
            m_existing_clique_var_indices[wit->clique_model_index];
        if (!push_vertex_pair(m_empty_propagator, m_all_vertices[vind],
                              m_all_vertices[wind]))
        {
            throw std::logic_error(
                "Incorrectly identified edge as violated nonedge");
        }
        m_covered_violated_nonedges[i].set();
        p_extend_nonedge_cover_greedy(i + 1, nve);
        p_extend_nonedge_cover_greedy(0, i);
        p_clq_add_partial_constraint(m_empty_propagator);
        m_empty_propagator.reset_to_zero();
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_extend_nonedge_cover_greedy(
    std::size_t begin_index, std::size_t end_index) {
    for (std::size_t j = begin_index; j < end_index; ++j) {
        if (m_covered_violated_nonedges[j])
            continue;
        auto [vit, wit] = m_clq_violated_nonedges[j];
        std::size_t ci1 = vit->clique_model_index;
        std::size_t ci2 = wit->clique_model_index;
        Vertex v1 = m_all_vertices[m_existing_clique_var_indices[ci1]];
        Vertex v2 = m_all_vertices[m_existing_clique_var_indices[ci2]];
        if (push_vertex_pair(m_empty_propagator, v1, v2)) {
            m_covered_violated_nonedges[j].set();
        }
    }
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_prohibit_extend_existing(
    OrderedSolutionIter vit, OrderedSolutionIter wit) {
    std::size_t ci1 = vit->clique_model_index;
    std::size_t ci2 = wit->clique_model_index;
    Vertex v1 = m_all_vertices[m_existing_clique_var_indices[ci1]];
    Vertex v2 = m_all_vertices[m_existing_clique_var_indices[ci2]];
    for (std::size_t i = 0, ninc = m_clq_constr_partial_constraints.size();
         i < ninc; ++i)
    {
        SharedDBPropagator& prop = m_clq_constr_partial_assignments[i];
        if (prop.is_false(v1.first) || prop.is_false(v1.second) ||
            prop.is_false(v2.first) || prop.is_false(v2.second))
            continue;
        std::size_t trail_length = prop.get_trail().size();
        if (push_vertex_pair(prop, v1, v2)) {
            p_clq_extended_partial(i, trail_length);
            return true;
        }
    }
    return false;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_round_greedy() {
    // track the vertices in the clique
    m_vertex_set_buffer.clear();
    // track possible extensions of the clique
    m_vertex_bitset_buffer.assign(m_all_vertices.size(), true);
    // greedy rounding procedure - we use the basic G_1 definition
    // that has often been extended, at least on the
    // vertices in the ordered solution, by violated non-edge detection
    for (const auto& ov : m_clq_ordered_solution) {
        std::size_t ci = ov.clique_model_index;
        std::size_t vi = m_existing_clique_var_indices[ci];
        if (!m_vertex_bitset_buffer[vi])
            continue;
        assert(!m_adjacency_matrix[vi][vi]);
        m_vertex_bitset_buffer &= m_adjacency_matrix[vi];
        m_vertex_set_buffer.push_back(vi);
    }
    p_randomly_extend_clique();
    if (m_vertex_set_buffer.size() > m_best_clique.size()) {
        m_best_clique.swap(m_vertex_set_buffer);
        if (m_best_clique.size() > m_lower_bound) {
            m_lower_bound = m_best_clique.size();
            m_lower_bound_subgraph = m_best_clique;
            p_call_lb_callback();
        }
        p_ensure_best_clique_in_model();
        return true;
    }
    return false;
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_randomly_extend_clique() {
    auto& rng = sammy::rng();
    const std::size_t n = m_all_vertices.size();
    std::uniform_int_distribution<std::size_t> indices(0, n - 1);
    for (std::size_t trial_counter = 1; trial_counter <= 100; ++trial_counter) {
        std::size_t vi = indices(rng);
        if (!m_vertex_bitset_buffer[vi])
            continue;
        assert(!m_adjacency_matrix[vi][vi]);
        m_vertex_bitset_buffer &= m_adjacency_matrix[vi];
        m_vertex_set_buffer.push_back(vi);
        trial_counter = 0;
    }
    bool any = m_vertex_bitset_buffer.any();
    while (any) {
        any = false;
        for (std::size_t vi : m_vertex_bitset_buffer.ones()) {
            any = true;
            assert(!m_adjacency_matrix[vi][vi]);
            m_vertex_bitset_buffer &= m_adjacency_matrix[vi];
            m_vertex_set_buffer.push_back(vi);
            break;
        }
    }
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<
    IncrementalSATSolver>::p_ensure_best_clique_in_model() {
    for (std::size_t vi : m_best_clique) {
        m_last_vertex_addition_source = "ensure_best_clique_in_model";
        p_clq_add_vertex(vi);
    }
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_clique_model_cclass_cuts() {
    bool found = false;
    for (const auto& cc : m_color_classes) {
        const auto& prop = cc.propagator;
        double sum = 0.0;
        m_buffer_expr.clear();
        for (std::size_t i = 0, nc = m_existing_clique_var_indices.size();
             i != nc; ++i)
        {
            std::size_t vi = m_existing_clique_var_indices[i];
            Vertex v = m_all_vertices[vi];
            if (prop.is_false(v.first) || prop.is_false(v.second))
                continue;
            assert(i < m_existing_clique_var_value_count);
            sum += m_existing_clique_var_values[i];
            m_buffer_expr += m_existing_clique_vars[i];
        }
        if (sum >= 1.01) {
            p_clq_add_partial_constraint(prop, true);
            found = true;
        }
    }
    return found;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_use_clique_model(
    const std::vector<std::size_t>& current_vertices) {
    std::size_t old_clique_size = m_best_clique.size();
    for (std::size_t i = 0; i < 5; ++i) {
        if (m_last_rounded_objective > m_best_clique.size() &&
            !p_cheap_cut_rounds())
            return false;
        if (m_last_rounded_objective > m_best_clique.size() &&
            !p_expensive_cut_rounds())
            return false;
        if (!p_limited_price_vertices(current_vertices))
            return false;
    }
    if (!p_unlimited_price_vertices(current_vertices))
        return false;
    if (m_last_rounded_objective > m_best_clique.size() &&
        !p_cheap_cut_rounds())
        return false;
    if (m_last_rounded_objective > m_best_clique.size() &&
        !p_expensive_cut_rounds())
        return false;
    return !m_aborted.load() && m_best_clique.size() > old_clique_size;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_cheap_cut_rounds() {
    bool repeat = true;
    std::size_t total_it = 0;
    while (repeat && m_last_rounded_objective > m_best_clique.size()) {
        repeat = false;
        double gap_before = p_clique_gap();
        for (std::size_t it = 0; it < CHEAP_CUT_ROUNDS_PER_GAP_CHECK;
             ++it, ++total_it)
        {
            if (!p_cheap_cut_round()) {
                return !m_aborted.load();
            }
            if (m_last_rounded_objective <= m_best_clique.size())
                return true;
        }
        double gap_after = p_clique_gap();
        repeat = (gap_before - gap_after >= CHEAP_CUT_GAP_REDUCTION_REQUIRED);
        repeat &= (total_it < CHEAP_CUT_ROUNDS_HARD_CAP);
    }
    return true;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_cheap_cut_round() {
    if (p_greedy_add_to_cuts() || p_greedy_generate_cuts() ||
        p_clique_model_cclass_cuts())
    {
        return p_solve_full_relaxation();
    }
    return false;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_greedy_add_to_cuts() {
    bool found = false;
    for (std::size_t partial_i : range(m_clq_constr_partial_assignments.size()))
    {
        auto& prop = m_clq_constr_partial_assignments[partial_i];
        double total_value = 0.0;
        for (const OrderedSolutionValue& ov : m_clq_ordered_solution) {
            std::size_t ci = ov.clique_model_index;
            std::size_t vi = m_existing_clique_var_indices[ci];
            Vertex v = m_all_vertices[vi];
            if (!can_push(prop, v))
                continue;
            total_value += ov.value;
        }
        if (total_value < 1.01)
            continue;
        total_value = 0.0;
        std::size_t old_trail_length = prop.get_trail().size();
        for (const OrderedSolutionValue& ov : m_clq_ordered_solution) {
            std::size_t ci = ov.clique_model_index;
            std::size_t vi = m_existing_clique_var_indices[ci];
            Vertex v = m_all_vertices[vi];
            if (push_vertex(prop, v) >= 0)
                total_value += ov.value;
        }
        p_clq_extended_partial(partial_i, old_trail_length);
        if (total_value >= 1.01)
            found = true;
    }
    return found;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_greedy_generate_cuts() {
    m_greedy_cut_candidates.clear();
    for (OrderedSolutionIter b = m_clq_ordered_solution.begin(), i = b,
                             e = m_clq_ordered_solution.end();
         i != e; ++i)
    {
        m_vertex_set_buffer.clear();
        assert(m_empty_propagator.get_current_level() == 0);
        double v = p_extend_greedy_cut(i, e);
        v += p_extend_greedy_cut(b, i);
        if (v >= 1.01 && m_greedy_cut_candidates.would_push(v)) {
            m_greedy_cut_candidates.emplace(v, m_vertex_set_buffer);
        }
        m_empty_propagator.reset_to_zero();
    }
    for (const GreedyCutCandidate& candidate :
         m_greedy_cut_candidates.elements())
    {
        for (std::size_t vertex : candidate.vertex_indices) {
            Vertex v = m_all_vertices[vertex];
            if (push_vertex(m_empty_propagator, v) < 0) {
                throw std::logic_error(
                    "Incorrect greedy cut candidate: produced conflict!");
            }
        }
        p_clq_add_partial_constraint(m_empty_propagator);
        m_empty_propagator.reset_to_zero();
    }
    return !m_greedy_cut_candidates.elements().empty();
}

template <typename IncrementalSATSolver>
double CliqueSatDSaturSolver<IncrementalSATSolver>::p_extend_greedy_cut(
    OrderedSolutionIter begin, OrderedSolutionIter end) {
    double value = 0.0;
    for (OrderedSolutionIter i = begin; i != end; ++i) {
        std::size_t ci = i->clique_model_index;
        std::size_t vi = m_existing_clique_var_indices[ci];
        Vertex v = m_all_vertices[vi];
        if (push_vertex(m_empty_propagator, v) < 0)
            continue;
        value += i->value;
        m_vertex_set_buffer.push_back(vi);
    }
    return value;
}

template <typename IncrementalSATSolver>
double CliqueSatDSaturSolver<IncrementalSATSolver>::p_weight_in_dual(Vertex v) {
    double value = 0.0;
    for (std::size_t complete_i :
         range(m_clq_constr_complete_assignments.size()))
    {
        const auto& assignment = m_clq_constr_complete_assignments[complete_i];
        if (lit::is_true_in(v.first, assignment) &&
            lit::is_true_in(v.second, assignment))
        {
            value += m_clq_constr_complete_dual_values[complete_i];
            if (value >= 1.02)
                return 1.02;
        }
    }
    for (std::size_t partial_i : range(m_clq_constr_partial_assignments.size()))
    {
        const auto& assignment = m_clq_constr_partial_assignments[partial_i];
        if (assignment.is_true(v.first) && assignment.is_true(v.second)) {
            value += m_clq_constr_partial_dual_values[partial_i];
            if (value >= 1.02)
                return 1.02;
        }
    }
    return value;
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::
    p_price_collect_good_and_possible(
        const std::vector<std::size_t>& vertices) {
    m_pricing_good_vertices.clear();
    m_pricing_possible_vertices.clear();
    for (std::size_t vertex : vertices) {
        if (m_present_in_clique_model[vertex])
            continue;
        Vertex v = m_all_vertices[vertex];
        double w = p_weight_in_dual(v);
        if (w >= 1.01)
            continue;
        if (w <= 0.99) {
            m_pricing_good_vertices.emplace_back(vertex, w);
        } else {
            m_pricing_possible_vertices.emplace_back(vertex, w);
        }
    }
    return !m_pricing_good_vertices.empty() ||
           !m_pricing_possible_vertices.empty();
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_limited_price_vertices(
    const std::vector<std::size_t>& vertices) {
    std::size_t goal_vertices = (std::min)(
        vertices.size(),
        (std::max)(std::size_t(10),
                   std::size_t(0.2 * m_existing_clique_vars.size())));
    if (m_aborted.load())
        return false;
    if (!p_price_collect_good_and_possible(vertices))
        return false;
    if (m_aborted.load())
        return false;
    p_pricing_select_with_goal(goal_vertices);
    p_pricing_add_vertices();
    return p_solve_full_relaxation();
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_pricing_select_with_goal(
    std::size_t goal_vertices) {
    auto& good = m_pricing_good_vertices;
    auto& possible = m_pricing_possible_vertices;
    m_vertex_set_buffer.clear();
    if (good.empty()) {
        // no good vertices; take possible vertices
        if (goal_vertices > possible.size())
            goal_vertices = possible.size();
        std::transform(possible.begin(), possible.begin() + goal_vertices,
                       std::back_inserter(m_vertex_set_buffer),
                       [](const PricingEntry& e) { return e.vertex_index; });
        return;
    }
    if (good.size() <= goal_vertices) {
        // not enough good vertices; take all
        std::transform(good.begin(), good.end(),
                       std::back_inserter(m_vertex_set_buffer),
                       [](const PricingEntry& e) { return e.vertex_index; });
        return;
    }
    auto goal_iter = good.begin() + goal_vertices;
    std::nth_element(good.begin(), goal_iter, good.end(),
                     [](const PricingEntry& e1, const PricingEntry& e2) {
                         return e1.dual_weight > e2.dual_weight;
                     });
    std::transform(good.begin(), goal_iter,
                   std::back_inserter(m_vertex_set_buffer),
                   [](const PricingEntry& e) { return e.vertex_index; });
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_pricing_add_vertices() {
    if (m_vertex_set_buffer.empty())
        return;
    m_last_vertex_addition_source = "pricing_add_vertices";
    for (std::size_t v : m_vertex_set_buffer)
        p_clq_add_vertex(v);
}

template <typename IncrementalSATSolver>
bool CliqueSatDSaturSolver<IncrementalSATSolver>::p_unlimited_price_vertices(
    const std::vector<std::size_t>& vertices) {
    if (m_aborted.load())
        return false;
    if (!p_price_collect_good_and_possible(vertices))
        return false;
    if (m_aborted.load())
        return false;
    auto* source = &m_pricing_good_vertices;
    if (m_pricing_good_vertices.empty())
        source = &m_pricing_possible_vertices;
    m_vertex_set_buffer.clear();
    std::transform(source->begin(), source->end(),
                   std::back_inserter(m_vertex_set_buffer),
                   [](const PricingEntry& e) { return e.vertex_index; });
    p_pricing_add_vertices();
    return p_solve_full_relaxation();
}

template <typename IncrementalSATSolver>
void CliqueSatDSaturSolver<IncrementalSATSolver>::p_call_lb_callback() {
    if (!m_lower_bound_callback)
        return;
    std::vector<Vertex> buffer;
    buffer.reserve(m_lower_bound_subgraph.size());
    std::transform(m_lower_bound_subgraph.begin(), m_lower_bound_subgraph.end(),
                   std::back_inserter(buffer),
                   [this](std::size_t vi) { return m_all_vertices[vi]; });
    m_lower_bound_callback(m_lower_bound, buffer);
}

} // namespace sammy

#endif
==> ./dynamic_bitset.h <==
#ifndef HS_DYNAMIC_BITSET_H_
#define HS_DYNAMIC_BITSET_H_

#if (defined(WIN32) || defined(_WIN32)) && defined(_MSC_VER) &&                \
    !defined(__clang__)
#include <intrin.h>
#pragma intrinsic(_BitScanForward64)
#pragma intrinsic(__popcnt64)
#endif

#include "range.h"
#include <algorithm>
#include <cassert>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <vector>

namespace sammy {

// the template stuff is necessary because the call might be ambiguous
// otherwise.
template <typename T, std::enable_if_t<std::is_integral_v<T> &&
                                           sizeof(T) == sizeof(std::uint32_t),
                                       int> = 0>
static inline auto count_set_bits(std::uint32_t x) noexcept {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcount(x);
#elif defined(_MSC_VER)
    return __popcnt(x);
#else
    std::uint32_t count_2 =
        (x & UINT32_C(0x5555'5555)) + ((x & UINT32_C(0xAAAA'AAAA)) >> 1);
    std::uint32_t count_4 = (count_2 & UINT32_C(0x3333'3333)) +
                            ((count_2 & UINT32_C(0xCCCC'CCCC)) >> 2);
    std::uint32_t count_8 = (count_4 & UINT32_C(0x0F0F'0F0F)) +
                            ((count_4 & UINT32_C(0xF0F0'F0F0)) >> 4);
    std::uint32_t count_w = (count_8 & UINT32_C(0x00FF'00FF)) +
                            ((count_8 & UINT32_C(0xFF00'FF00)) >> 8);
    return (count_w & UINT32_C(0x0000'FFFF)) +
           ((count_w & UINT32_C(0xFFFF'0000)) >> 16);
#endif
}

// the template stuff is necessary because the call might be ambiguous
// otherwise.
template <typename T, std::enable_if_t<std::is_integral_v<T> &&
                                           sizeof(T) == sizeof(std::uint64_t),
                                       int> = 0>
static inline auto count_set_bits(T x) noexcept {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(x);
#elif defined(_MSC_VER)
    return __popcnt64(x);
#else
    std::uint64_t count_2 = (x & UINT64_C(0x5555'5555'5555'5555)) +
                            ((x & UINT64_C(0xAAAA'AAAA'AAAA'AAAA)) >> 1);
    std::uint64_t count_4 = (count_2 & UINT64_C(0x3333'3333'3333'3333)) +
                            ((count_2 & UINT64_C(0xCCCC'CCCC'CCCC'CCCC)) >> 2);
    std::uint64_t count_8 = (count_4 & UINT64_C(0x0F0F'0F0F'0F0F'0F0F)) +
                            ((count_4 & UINT64_C(0xF0F0'F0F0'F0F0'F0F0)) >> 4);
    std::uint64_t count_w = (count_8 & UINT64_C(0x00FF'00FF'00FF'00FF)) +
                            ((count_8 & UINT64_C(0xFF00'FF00'FF00'FF00)) >> 8);
    std::uint64_t count_d = (count_w & UINT64_C(0x0000'FFFF'0000'FFFF)) +
                            ((count_w & UINT64_C(0xFFFF'0000'FFFF'0000)) >> 16);
    return (count_d & UINT64_C(0x0000'0000'FFFF'FFFF)) +
           ((count_d & UINT64_C(0xFFFF'FFFF'0000'0000)) >> 32);
#endif
}

template <typename T, std::enable_if_t<std::is_integral_v<T> &&
                                           sizeof(T) == sizeof(std::uint32_t),
                                       int> = 0>
static inline auto count_trailing_zeros(T x) noexcept {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctz(x);
#elif defined(_MSC_VER)
    unsigned long result;
    _BitScanForward(&result, x);
    return result;
#else
    std::size_t result = 0;
    for (std::uint32_t mask = 1; result < 32; ++result, mask <<= 1) {
        if (x & mask)
            break;
    }
    return result;
#endif
}

template <typename T, std::enable_if_t<std::is_integral_v<T> &&
                                           sizeof(T) == sizeof(std::uint64_t),
                                       int> = 0>
static inline auto count_trailing_zeros(T x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(x);
#elif defined(_MSC_VER)
    unsigned long result;
    _BitScanForward64(&result, x);
    return result;
#else
    std::size_t result = 0;
    for (std::uint64_t mask = 1; result < 64; ++result, mask <<= 1) {
        if (x & mask)
            break;
    }
    return result;
#endif
}

/**
 * @brief Since std::vector<bool> lacks efficient counting/iteration,
 *        and boost::dynamic_bitset does not do doubling growth
 *        unless std::vector::resize does, here goes another bitset.
 *        Also, we are hopefully more auto-vectorization-friendly
 *        than boost::dynamic_bitset, which struggles to get
 *        auto-vectorized (at least by GCC 11) for operations
 *        like -= and such.
 */
class DynamicBitset {
  public:
    using Block = std::size_t;
    static constexpr std::size_t bits_per_word = sizeof(Block) * CHAR_BIT;
    using const_reference = bool;
    using size_type = std::size_t;

    static constexpr std::size_t
    num_blocks_required(std::size_t n_bits) noexcept {
        return (n_bits / bits_per_word) + bool(n_bits % bits_per_word);
    }

    static constexpr Block broadcast(bool v) noexcept {
        // usually, the compiler is clever enough for this to
        // get compiled to sensible branch-free code
        return v ? ~Block(0) : Block(0);
    }

    // Constructors
    DynamicBitset() noexcept : m_num_bits(0) {}

    explicit DynamicBitset(std::size_t s, bool value = false)
        : m_blocks(num_blocks_required(s), broadcast(value)), m_num_bits(s) {
        p_zero_unused();
    }

    explicit DynamicBitset(const std::vector<bool>& vbool)
        : m_blocks(num_blocks_required(vbool.size()), broadcast(false)),
          m_num_bits(vbool.size()) {
        for (std::size_t i = 0, s = vbool.size(); i != s; ++i) {
            if (vbool[i])
                (*this)[i].set();
        }
    }

    bool operator==(const DynamicBitset& o) const noexcept {
        return size() == o.size() &&
               std::equal(m_blocks.begin(), m_blocks.end(), o.m_blocks.begin());
    }

    bool operator!=(const DynamicBitset& o) const noexcept {
        return !(*this == o);
    }

    DynamicBitset& operator|=(const DynamicBitset& o) noexcept {
        assert(m_num_bits == o.m_num_bits);
        std::transform(m_blocks.begin(), m_blocks.end(), o.m_blocks.begin(),
                       m_blocks.begin(),
                       [](Block b1, Block b2) { return b1 | b2; });
        return *this;
    }

    void binary_or(const DynamicBitset& o, std::size_t begin_offset) {
        assert(begin_offset < m_num_bits);
        assert(m_num_bits == o.m_num_bits);
        std::size_t bblock(begin_offset / bits_per_word);
        std::uint8_t boffs(begin_offset % bits_per_word);
        Block o_first = o.m_blocks[bblock];
        o_first &= ~Block((Block(1) << boffs) - 1);
        m_blocks[bblock] |= o_first;
        ++bblock;
        std::transform(m_blocks.begin() + bblock, m_blocks.end(),
                       o.m_blocks.begin() + bblock, m_blocks.begin() + bblock,
                       [](Block b1, Block b2) { return b1 | b2; });
    }

    void binary_subtract(const DynamicBitset& o, std::size_t begin_offset) {
        assert(begin_offset < m_num_bits);
        assert(m_num_bits <= o.m_num_bits);
        std::size_t bblock(begin_offset / bits_per_word);
        std::size_t boffs(begin_offset % bits_per_word);
        Block o_first = o.m_blocks[bblock];
        o_first &= ~Block((Block(1) << boffs) - 1);
        m_blocks[bblock] &= ~o_first;
        ++bblock;
        std::transform(m_blocks.begin() + bblock, m_blocks.end(),
                       o.m_blocks.begin() + bblock, m_blocks.begin() + bblock,
                       [](Block b1, Block b2) { return b1 & ~b2; });
    }

    void binary_subtract(const std::uint32_t* bits_begin) {
        constexpr std::size_t blocks_per_u32 =
            sizeof(Block) / sizeof(std::uint32_t);
        const std::size_t nblocks = m_blocks.size();
        for (std::size_t i = 0; i < nblocks; ++i) {
            Block blk = p_read_block(bits_begin + i * blocks_per_u32);
            m_blocks[i] &= ~blk;
        }
    }

    void binary_or(const std::uint32_t* bits_begin) {
        constexpr std::size_t blocks_per_u32 =
            sizeof(Block) / sizeof(std::uint32_t);
        const std::size_t nblocks = m_blocks.size();
        for (std::size_t i = 0; i < nblocks; ++i) {
            Block blk = p_read_block(bits_begin + i * blocks_per_u32);
            m_blocks[i] |= blk;
        }
    }

    void binary_and(const std::uint32_t* bits_begin) {
        constexpr std::size_t blocks_per_u32 =
            sizeof(Block) / sizeof(std::uint32_t);
        const std::size_t nblocks = m_blocks.size();
        for (std::size_t i = 0; i < nblocks; ++i) {
            Block blk = p_read_block(bits_begin + i * blocks_per_u32);
            m_blocks[i] &= blk;
        }
    }

    void binary_and(const DynamicBitset& o, std::size_t begin_offset) {
        assert(begin_offset < m_num_bits);
        assert(m_num_bits <= o.m_num_bits);
        std::size_t bblock(begin_offset / bits_per_word);
        std::uint8_t boffs(begin_offset % bits_per_word);
        Block o_first = o.m_blocks[bblock];
        o_first |= Block((Block(1) << boffs) - 1);
        m_blocks[bblock] &= o_first;
        ++bblock;
        std::transform(m_blocks.begin() + bblock, m_blocks.end(),
                       o.m_blocks.begin() + bblock, m_blocks.begin() + bblock,
                       [](Block b1, Block b2) { return b1 & b2; });
    }

    DynamicBitset& operator&=(const DynamicBitset& o) noexcept {
        assert(m_num_bits <= o.m_num_bits);
        std::transform(m_blocks.begin(), m_blocks.end(), o.m_blocks.begin(),
                       m_blocks.begin(),
                       [](Block b1, Block b2) { return b1 & b2; });
        return *this;
    }

    DynamicBitset& operator^=(const DynamicBitset& o) noexcept {
        assert(m_num_bits == o.m_num_bits);
        std::transform(m_blocks.begin(), m_blocks.end(), o.m_blocks.begin(),
                       m_blocks.begin(),
                       [](Block b1, Block b2) { return b1 ^ b2; });
        return *this;
    }

    DynamicBitset& operator-=(const DynamicBitset& o) noexcept {
        assert(m_num_bits == o.m_num_bits);
        std::transform(m_blocks.begin(), m_blocks.end(), o.m_blocks.begin(),
                       m_blocks.begin(),
                       [](Block b1, Block b2) { return b1 & ~b2; });
        return *this;
    }

    DynamicBitset operator~() const {
        DynamicBitset result(*this);
        result.flip();
        return result;
    }

    // Pseudo-reference
    class reference {
      public:
        reference(const reference&) noexcept = default;

        /* implicit */ operator bool() const noexcept { return *m_ptr & m_msk; }

        bool operator~() const noexcept { return !(*m_ptr & m_msk); }

        bool operator!() const noexcept { return ~*this; }

        reference& flip() noexcept {
            *m_ptr ^= m_msk;
            return *this;
        }

        reference& set() noexcept {
            *m_ptr |= m_msk;
            return *this;
        }

        reference& set(bool v) noexcept { return *this = v; }

        reference& reset() noexcept {
            *m_ptr &= ~m_msk;
            return *this;
        }

        reference& operator=(bool v) noexcept {
            if (v)
                *m_ptr |= m_msk;
            else
                *m_ptr &= ~m_msk;
            return *this;
        }

        reference& operator=(const reference& o) noexcept {
            return *this = static_cast<bool>(o);
        }

        reference& operator|=(bool x) noexcept {
            if (x)
                *this = true;
            return *this;
        }

        reference& operator&=(bool x) noexcept {
            if (!x)
                *this = false;
            return *this;
        }

        reference& operator^=(bool x) noexcept {
            if (x)
                this->flip();
            return *this;
        }

        reference& operator-=(bool x) noexcept {
            if (x)
                this->reset();
            return *this;
        }

      private:
        void operator&() = delete;

        reference(Block* ptr, Block msk) noexcept : m_ptr(ptr), m_msk(msk) {}

        Block* m_ptr;
        Block m_msk;

        friend class DynamicBitset;
    };

    reference operator[](std::size_t idx) noexcept {
        std::size_t blk_i = idx / bits_per_word;
        std::uint8_t sub_i(idx % bits_per_word);
        return reference{m_blocks.data() + blk_i, Block(1) << sub_i};
    }

    bool operator[](std::size_t idx) const noexcept {
        std::size_t blk_i = idx / bits_per_word;
        std::uint8_t sub_i(idx % bits_per_word);
        return m_blocks[blk_i] & (Block(1) << sub_i);
    }

    std::size_t size() const noexcept { return m_num_bits; }

    void reserve(size_type bits) {
        m_blocks.reserve(num_blocks_required(bits));
    }

    void resize(size_type num_bits, bool value = false) {
        auto nreq = num_blocks_required(num_bits);
        if (num_bits < m_num_bits) {
            m_blocks.resize(nreq);
        } else if (num_bits > m_num_bits) {
            m_blocks.reserve(nreq);
            if (value)
                p_one_unused();
            m_blocks.resize(nreq, broadcast(value));
        }
        m_num_bits = num_bits;
        p_zero_unused();
    }

    void assign(size_type num_bits, bool value) {
        if (num_bits != size())
            resize(num_bits, value);
        set(value);
    }

    void clear() noexcept {
        m_num_bits = 0;
        m_blocks.clear();
    }

    void push_back(bool bit) {
        std::size_t new_bit_subind = m_num_bits % bits_per_word;
        if (!new_bit_subind) {
            m_blocks.push_back(bit);
        } else {
            Block msk = bit;
            msk <<= new_bit_subind;
            m_blocks.back() |= msk;
        }
        ++m_num_bits;
    }

    void pop_back() noexcept {
        std::uint8_t bits_used_in_last(--m_num_bits % bits_per_word);
        if (bits_used_in_last == 0) {
            m_blocks.pop_back();
        } else {
            m_blocks.back() &= (Block(1) << bits_used_in_last) - 1;
        }
    }

    void set(bool value = true) noexcept {
        std::fill(m_blocks.begin(), m_blocks.end(), broadcast(value));
        p_zero_unused();
    }

    void reset() noexcept { set(false); }

    void flip() noexcept {
        std::transform(m_blocks.begin(), m_blocks.end(), m_blocks.begin(),
                       [](Block b) { return ~b; });
        p_zero_unused();
    }

    bool any() const noexcept {
        return std::any_of(m_blocks.begin(), m_blocks.end(),
                           [](Block b) { return b != 0; });
    }

    bool none() const noexcept { return !any(); }

    bool all() const noexcept {
        auto l = std::find_if(m_blocks.begin(), m_blocks.end(),
                              [](Block b) { return b != ~Block(0); });
        if (l == m_blocks.end())
            return true;
        if (l != m_blocks.end() - 1)
            return false;
        Block last = m_blocks.back();
        std::uint8_t bits_used_in_last = m_num_bits % bits_per_word;
        if (!bits_used_in_last)
            return false;
        return (last | ~Block((Block(1) << bits_used_in_last) - 1)) ==
               ~Block(0);
    }

    std::size_t count() const noexcept {
        return std::transform_reduce(m_blocks.begin(), m_blocks.end(),
                                     std::size_t(0), std::plus<>{},
                                     [](Block b) { return count_set_bits(b); });
    }

    std::size_t count_from(std::size_t begin_index) const noexcept {
        std::size_t word_idx = begin_index / bits_per_word;
        std::size_t sub_idx = begin_index % bits_per_word;
        Block partial = m_blocks[word_idx];
        partial &= ~Block((Block(1) << sub_idx) - 1);
        std::size_t initial = count_set_bits(partial);
        return std::transform_reduce(m_blocks.begin() + word_idx + 1,
                                     m_blocks.end(), initial, std::plus<>{},
                                     [](Block b) { return count_set_bits(b); });
    }

    class OnesIterator {
      public:
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using reference = std::size_t;
        using pointer = const std::size_t*;
        using value_type = std::size_t;

        bool operator==(const OnesIterator& o) const noexcept {
            return m_pos == o.m_pos && m_cnt == o.m_cnt;
        }

        bool operator!=(const OnesIterator& o) const noexcept {
            return m_pos != o.m_pos || m_cnt != o.m_cnt;
        }

        OnesIterator operator++(int) const noexcept {
            OnesIterator result(*this);
            ++result;
            return result;
        }

        OnesIterator& operator++() noexcept {
            if (--m_cnt == 0) {
                do {
                    m_offs += bits_per_word;
                    if (++m_pos == m_end) {
                        m_cnt = 0;
                        return *this;
                    }
                    m_buf = *m_pos;
                } while (!m_buf);
                m_cnt = count_set_bits(m_buf);
            } else {
                m_buf &= m_buf - 1;
            }
            return *this;
        }

        std::size_t operator*() const noexcept {
            return count_trailing_zeros(m_buf) + m_offs;
        }

        OnesIterator() noexcept = default;

      private:
        friend class DynamicBitset;

        explicit OnesIterator(const Block* end) noexcept
            : m_pos(end), m_end(end), m_buf(0), m_cnt(0), m_offs(0) {}

        explicit OnesIterator(const Block* begin, const Block* end) noexcept
            : m_pos(begin), m_end(end), m_buf(0), m_cnt(0), m_offs(0) {
            p_scroll();
        }

        explicit OnesIterator(const Block* begin, const Block* end,
                              std::uint8_t index_in_word,
                              std::size_t offs) noexcept
            : m_pos(begin), m_end(end), m_buf(0), m_cnt(0), m_offs(offs) {
            if (m_pos != m_end) {
                m_buf = *m_pos;
                m_buf &= ~Block((Block(1) << index_in_word) - 1);
                if (m_buf) {
                    m_cnt = count_set_bits(m_buf);
                    return;
                }
                ++m_pos;
                m_offs += bits_per_word;
                p_scroll();
            }
        }

        void p_scroll() {
            while (m_pos != m_end) {
                m_buf = *m_pos;
                if (m_buf) {
                    m_cnt = count_set_bits(m_buf);
                    return;
                }
                ++m_pos;
                m_offs += bits_per_word;
            }
        }

        const Block* m_pos;
        const Block* m_end;
        Block m_buf;
        Block m_cnt;
        std::size_t m_offs;
    };

    OnesIterator ones_begin() const noexcept {
        return OnesIterator(m_blocks.data(), m_blocks.data() + m_blocks.size());
    }

    OnesIterator ones_end() const noexcept {
        return OnesIterator(m_blocks.data() + m_blocks.size());
    }

    IteratorRange<OnesIterator> ones() const noexcept {
        return {ones_begin(), ones_end()};
    }

    OnesIterator ones_from_begin(std::size_t begin_index) const noexcept {
        std::size_t word_idx = begin_index / bits_per_word;
        std::size_t sub_idx = begin_index % bits_per_word;
        std::size_t offs = word_idx * bits_per_word;
        return OnesIterator(m_blocks.data() + word_idx,
                            m_blocks.data() + m_blocks.size(), sub_idx, offs);
    }

    IteratorRange<OnesIterator>
    ones_from(std::size_t begin_index) const noexcept {
        return {ones_from_begin(begin_index), ones_end()};
    }

    std::size_t bytes_used() const noexcept {
        return m_blocks.capacity() * sizeof(Block) * CHAR_BIT / 8;
    }

    DynamicBitset(const DynamicBitset&) = default;
    DynamicBitset& operator=(const DynamicBitset&) = default;
    ~DynamicBitset() = default;

    DynamicBitset(DynamicBitset&& o) noexcept
        : m_blocks(std::move(o.m_blocks)), m_num_bits(o.m_num_bits) {
        o.m_blocks.clear();
        o.m_num_bits = 0;
    }

    DynamicBitset& operator=(DynamicBitset&& o) noexcept {
        std::swap(m_num_bits, o.m_num_bits);
        m_blocks.swap(o.m_blocks);
        return *this;
    }

    explicit operator std::vector<bool>() const {
        const std::size_t n = size();
        std::vector<bool> result(n, false);
        for (std::size_t i = 0; i < n; ++i) {
            result[i] = bool((*this)[i]);
        }
        return result;
    }

    const std::vector<Block>& blocks() const noexcept { return m_blocks; }

  private:
    template <typename T, std::enable_if_t<sizeof(T) == sizeof(Block), int> = 0>
    Block p_read_block(const T* from) {
        return Block(*from);
    }

    template <typename T,
              std::enable_if_t<sizeof(T) * 2 == sizeof(Block) &&
                                   sizeof(Block) == sizeof(std::uint64_t),
                               int> = 0>
    Block p_read_block(const T* from) {
        Block b1(*from);
        Block b2(*(from + 1));
        return (b2 << 32) | b1;
    }

    void p_zero_unused() {
        std::uint8_t bits_used_in_last = m_num_bits % bits_per_word;
        if (!bits_used_in_last)
            return;
        m_blocks.back() &= (Block(1) << bits_used_in_last) - 1;
    }

    void p_one_unused() {
        std::uint8_t bits_used_in_last = m_num_bits % bits_per_word;
        if (!bits_used_in_last)
            return;
        m_blocks.back() |= ~Block((Block(1) << bits_used_in_last) - 1);
    }

    std::vector<Block> m_blocks;
    std::size_t m_num_bits;
};

using Bitset = DynamicBitset;

} // namespace sammy

#endif
==> ./pair_infeasibility_map.h <==
#ifndef SAMMY_PAIR_INFEASIBILITY_MAP_H_
#define SAMMY_PAIR_INFEASIBILITY_MAP_H_

#include "clause_db.h"
#include "cuda_iteration.h"
#include "dynamic_bitset.h"
#include "literals.h"
#include "rng.h"
#include "shared_db_propagator.h"
#include "thread_group.h"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <memory>
#include <new>

namespace sammy {

using Index = Lit;

inline Bitset make_bitset(std::size_t length, bool value = false) {
    return Bitset(length, value);
}

inline Bitset make_larger_bitset(const Bitset& from) {
    Bitset bs;
    bs.reserve(from.size() * 2);
    bs = from;
    return bs;
}

#ifdef SAMMY_CUDA_SUPPORTED
static inline void cuda_extract_feasibilities(
    Var n_concrete, std::vector<DynamicBitset>& definitely_feasible_matrix,
    const std::vector<DynamicBitset>& literals_in_class) {
    if (literals_in_class.empty())
        return;

    const Lit nclit = 2 * n_concrete;
    if (definitely_feasible_matrix.size() != nclit ||
        definitely_feasible_matrix[0].size() != nclit)
    {
        throw std::invalid_argument(
            "definitely_feasible_matrix has wrong size");
    }
    std::size_t u32_per_bitset = literals_in_class[0].blocks().size() *
                                 sizeof(DynamicBitset::Block) /
                                 sizeof(std::uint32_t);
    std::vector<std::uint32_t> host_prepare_buffer;
    host_prepare_buffer.reserve(literals_in_class.size() * u32_per_bitset);
    for (const DynamicBitset& bset : literals_in_class) {
        for (DynamicBitset::Block b : bset.blocks()) {
            detail::to_prepare_buffer(b, host_prepare_buffer);
        }
    }
    detail::CUDADevicePointer<const std::uint32_t> device_bit_data(
        host_prepare_buffer);
    detail::CUDADevicePointer<std::uint32_t> device_output_buffer(
        nclit * u32_per_bitset);
    Var v = 0;
    while (v < n_concrete) {
        Var num_rows = n_concrete - v;
        if (num_rows > 2 * detail::GOAL_ROWS_PER_CALL()) {
            num_rows = 2 * detail::GOAL_ROWS_PER_CALL();
        }
        detail::call_cuda_extract_kernel(device_bit_data.get(), u32_per_bitset,
                                         device_output_buffer.get(), nclit, v,
                                         num_rows, literals_in_class.size());
        v += num_rows;
    }
    std::vector<std::uint32_t> host_output =
        device_output_buffer.get_host_copy();
    for (Lit l = 0; l < nclit; ++l) {
        definitely_feasible_matrix[l].binary_or(
            &host_output[l * u32_per_bitset]);
    }
}
#endif

/**
 * Matrix that stores pairs that are
 * known-feasible or known-infeasible.
 */
class PairInfeasibilityMap {
  public:
    using Row = Bitset;
    using Matrix = std::vector<Row>;

  private:
    std::size_t num_vars;
    Matrix m_matrix;
    Matrix m_def_feasible;
    std::vector<Lit> m_incorporate_buffer;

    explicit PairInfeasibilityMap(std::size_t n_concrete,
                                  const std::vector<std::vector<bool>>& bits)
        : num_vars(n_concrete), m_matrix(), m_def_feasible() {
        m_matrix.reserve(bits.size());
        m_def_feasible.reserve(bits.size());
        for (std::size_t i = 0; i < bits.size(); ++i) {
            m_matrix.emplace_back(bits[i]);
            m_def_feasible.emplace_back(m_matrix.back());
            m_def_feasible.back().flip();
        }
    }

  public:
    std::size_t get_n_concrete() const noexcept { return num_vars; }

    explicit PairInfeasibilityMap(std::size_t n_concrete)
        : num_vars(n_concrete),
          m_matrix(2 * num_vars, Bitset(2 * num_vars, false)),
          m_def_feasible(2 * num_vars, Bitset(2 * num_vars, false)) {
        Lit l = 0;
        for (auto& r : m_matrix) {
            r[lit::negate(l)] = true;
            ++l;
        }
    }

    Row& operator[](Lit l) noexcept { return m_matrix[l]; }

    const Row& operator[](Lit l) const noexcept { return m_matrix[l]; }

    void literal_infeasible(Lit l) noexcept {
        const auto nl = 2 * num_vars;
        m_matrix[l].set();
        for (Lit i = 0; i < nl; ++i) {
            m_matrix[i][l] = true;
        }
    }

    void literal_pair_infeasible(Lit l1, Lit l2) noexcept {
        m_matrix[l1][l2] = true;
        m_matrix[l2][l1] = true;
    }

    bool operator()(Lit lmin, Lit lmax) const noexcept {
        return m_matrix[lmin][lmax];
    }

    bool is_definitely_feasible(Lit lmin, Lit lmax) const noexcept {
        return m_def_feasible[lmin][lmax];
    }

    void set_definitely_feasible(Lit lmin, Lit lmax) noexcept {
        m_def_feasible[lmin][lmax] = true;
        m_def_feasible[lmax][lmin] = true;
    }

    std::vector<std::vector<bool>> export_bits() const {
        std::vector<std::vector<bool>> result;
        result.reserve(m_matrix.size());
        for (const auto& row : m_matrix) {
            result.emplace_back(static_cast<std::vector<bool>>(row));
        }
        return result;
    }

    static PairInfeasibilityMap
    import_bits(const std::vector<std::vector<bool>>& bits) {
        if (bits.size() == 0)
            throw std::runtime_error("Empty matrix in import_bits!");
        if (bits.size() % 2 == 1)
            throw std::runtime_error("Odd matrix in import_bits!");
        if (bits[0].size() != bits.size())
            throw std::runtime_error("Non-square matrix in import_bits!");
        std::size_t num_concrete = bits.size() / 2;
        return PairInfeasibilityMap{num_concrete, bits};
    }

    void incorporate_complete_class(const SharedDBPropagator& propagator) {
        const Lit nc_lit = 2 * num_vars;
        m_incorporate_buffer.clear();
        std::copy_if(propagator.get_trail().begin(),
                     propagator.get_trail().end(),
                     std::back_inserter(m_incorporate_buffer),
                     [&](Lit l) { return l < nc_lit; });
        std::sort(m_incorporate_buffer.begin(), m_incorporate_buffer.end());
        for (auto i = m_incorporate_buffer.begin(),
                  e = m_incorporate_buffer.end();
             i != e; ++i)
        {
            auto& row = m_def_feasible[*i];
            for (auto j = m_incorporate_buffer.begin(); j != i; ++j) {
                row[*j] = true;
            }
            for (auto j = std::next(i); j != e; ++j) {
                row[*j] = true;
            }
        }
    }

    void incorporate_complete_classes(
        const std::vector<DynamicBitset>& literals_in_class,
        ThreadGroup<void>& tpool) {
        auto handle_var = [&](Var v) {
            Lit plit = lit::positive_lit(v);
            DynamicBitset& positive_row = m_def_feasible[plit];
            DynamicBitset& negative_row = m_def_feasible[lit::negative_lit(v)];
            for (const DynamicBitset& bset : literals_in_class) {
                DynamicBitset& out_row =
                    (bset[plit] ? positive_row : negative_row);
                out_row |= bset;
            }
        };

        tpool.parallel_foreach_iterator(Var(0), Var(num_vars), handle_var);
    }

#ifdef SAMMY_CUDA_SUPPORTED
    void cuda_incorporate_complete_classes(
        const std::vector<DynamicBitset>& literals_in_class,
        ThreadGroup<void>& tpool) {
        try {
            cuda_extract_feasibilities(num_vars, m_def_feasible,
                                       literals_in_class);
            return;
        } catch (const CUDAError& err) {
            std::cerr << "CUDA error in cuda_incorporate_complete_classes: "
                      << err.what() << std::endl;
            had_cuda_error(err);
        }
        incorporate_complete_classes(literals_in_class, tpool);
    }
#endif

    /*
        template<typename PropIterator>
        void incorporate_complete_classes(PropIterator begin, PropIterator end,
       ThreadGroup<void>& tpool) { std::vector<Lit> prop_trails; std::size_t
       nclasses = std::size_t(end - begin); const Var ncvars = num_vars; const
       Lit nclits = 2 * ncvars; prop_trails.reserve(nclasses * num_vars);
            std::for_each(begin, end, [&] (const SharedDBPropagator& prop) {
                const auto& t = prop.get_trail();
                std::copy_if(t.begin(), t.end(),
       std::back_inserter(prop_trails),
                             [&] (Lit l) { return l < nclits; });
            });
            auto trail_range = [&] (std::size_t class_index) {
                auto tbegin = prop_trails.begin() + (class_index * ncvars);
                return std::make_pair(tbegin, tbegin + ncvars);
            };
            auto handle_var = [&] (Var v) {
                const Lit p = lit::positive_lit(v), n = lit::negative_lit(v);
                for(std::size_t iclass = 0; iclass < nclasses; ++iclass) {
                    Row& out = begin[iclass].is_true(p) ? m_def_feasible[p] :
       m_def_feasible[n]; auto [tbegin, tend] = trail_range(iclass);
                    std::for_each(tbegin, tend, [&] (Lit l) {out[l].set();});
                }
            };
            tpool.parallel_foreach_iterator(Var(0), ncvars, handle_var);
        } */

    std::size_t total_memory_usage() const noexcept {
        return 2 * (m_matrix.capacity() * sizeof(Row) +
                    m_matrix[0].bytes_used() * m_matrix.size()) +
               sizeof(PairInfeasibilityMap);
    }

    std::size_t count_vertices() const noexcept {
        auto res = std::transform_reduce(
            m_def_feasible.begin(), m_def_feasible.end(), std::size_t(0),
            std::plus<>{}, [](const Row& r) { return r.count(); });
        res /= 2;
        return res;
    }

    std::vector<Vertex>
    collect_vertices(std::size_t reserve_size = 0) const noexcept {
        std::vector<Vertex> result;
        if (reserve_size)
            result.reserve(reserve_size);
        for (Lit i = 0, nclit = 2 * num_vars; i < nclit - 2; ++i) {
            std::transform(
                m_def_feasible[i].ones_from_begin(i + 1),
                m_def_feasible[i].ones_end(), std::back_inserter(result),
                [&](std::size_t lmax) { return Vertex{i, Lit(lmax)}; });
        }
        return result;
    }

    std::vector<Vertex> sample_vertices(std::size_t target_size,
                                        std::size_t total_size = 0) {
        if (!total_size) {
            total_size = count_vertices();
        }
        double sample_prob = double(target_size) / total_size;
        sample_prob = std::min(1.0, sample_prob);
        std::geometric_distribution<std::size_t> dist{sample_prob};
        auto& rng = sammy::rng();
        std::size_t skip_count = dist(rng);
        std::vector<Vertex> result;
        result.reserve(std::size_t(1.1 * target_size));
        for (Lit i = 0, nclit = 2 * num_vars; i < nclit - 2; ++i) {
            for (std::size_t j : m_def_feasible[i].ones_from(i + 1)) {
                if (!skip_count) {
                    skip_count = dist(rng);
                    result.emplace_back(i, Lit(j));
                } else {
                    --skip_count;
                }
            }
        }
        return result;
    }
};

} // namespace sammy

#endif
==> ./problem_input.h <==
==> ./output.h <==
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
#include <numeric>
#include <optional>
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
==> ./fast_clique.h <==
#ifndef SAMMY_FAST_CLIQUE_H_INCLUDED_
#define SAMMY_FAST_CLIQUE_H_INCLUDED_

#include "literals.h"
#include "rng.h"
#include "shared_db_propagator.h"
#include "thread_group.h"

namespace sammy {

/**
 * @brief Quickly compute some maximal clique (i.e. a set of mutually exclusive
 * literal pairs) on a subgraph of vertices (i.e., literal pairs). The vertices
 * are assumed to be feasible if they can be pushed to a propagator; normally,
 * with vertices extracted from a sample, this should be a non-issue. Running
 * this heuristic should be feasible even if storing/computing the complete
 * subgraph is too expensive. Mainly useful to quickly get lower bounds.
 */
class FastCliqueBuilder {
  public:
    explicit FastCliqueBuilder(SharedDBPropagator prop)
        : m_prop(std::move(prop)) {}

    std::vector<Vertex> compute_clique(Vertex start_vertex,
                                       const std::vector<Vertex>& subgraph) {
        p_check_prop();
        m_subg = subgraph;
        if (!p_can_push(start_vertex))
            return {};
        p_filter_invalid();
        return p_compute_clique(start_vertex);
    }

    /**
     * @brief Compute a clique starting from a set of vertices. All vertices
     *        (starting vertices and subgraph) are assumed to be feasible
     * interactions.
     */
    template <typename StartVerticesIterator>
    std::vector<Vertex>
    compute_clique_known_valid(StartVerticesIterator begin,
                               StartVerticesIterator end,
                               const std::vector<Vertex>& subgraph) {
        std::vector<Vertex> result;
        p_check_prop();
        m_subg = subgraph;
        for (auto v : IteratorRange{begin, end}) {
            result.push_back(v);
            p_filter_nonneighbors(v);
        }
        while (!m_subg.empty()) {
            Vertex v = p_select_next_vertex();
            result.push_back(v);
            p_filter_nonneighbors(v);
        }
        return result;
    }

    std::vector<Vertex>
    random_multistart_best_clique(std::size_t iterations,
                                  const std::vector<Vertex>& subgraph) {
        p_check_prop();
        std::vector<Vertex> result;
        m_subg = subgraph;
        p_filter_invalid();
        std::vector<Vertex> filtered_whole = m_subg;
        for (std::size_t i = 0; i < iterations; ++i) {
            m_subg = filtered_whole;
            Vertex vbegin = p_select_next_vertex();
            std::vector<Vertex> current = p_compute_clique(vbegin);
            if (current.size() > result.size())
                result = current;
        }
        return result;
    }

    std::vector<Vertex> random_multistart_best_clique_known_valid(
        std::size_t iterations, const std::vector<Vertex>& subgraph) {
        p_check_prop();
        std::vector<Vertex> result;
        for (std::size_t i = 0; i < iterations; ++i) {
            m_subg = subgraph;
            Vertex vbegin = p_select_next_vertex();
            std::vector<Vertex> current = p_compute_clique(vbegin);
            if (current.size() > result.size())
                result = current;
        }
        return result;
    }

    // compute degree sequence (deg(vertices[i]) == result[i]);
    // much slower than random_multistart_best_clique!
    std::vector<std::size_t>
    subgraph_degrees(const std::vector<Vertex>& vertices) {
        p_check_prop();
        std::vector<std::size_t> result;
        for (auto i = vertices.begin(), e = vertices.end(); i != e; ++i) {
            Vertex v_i = *i;
            if (!p_push(v_i)) {
                result.push_back(0);
            } else {
                auto count = static_cast<std::size_t>(std::count_if(
                    vertices.begin(), vertices.end(), [&](Vertex v_j) {
                        return v_i != v_j && !p_can_push(v_j);
                    }));
                result.push_back(count);
            }
            while (m_prop.get_current_level() > 0)
                m_prop.pop_level();
        }
        return result;
    }

    std::vector<double>
    subgraph_degree_estimates(const std::vector<Vertex>& vertices,
                              std::size_t num_edge_samples) {
        p_check_prop();
        if (vertices.empty())
            return {};
        if (vertices.size() == 1)
            return std::initializer_list<double>{0.0};
        std::size_t total_edges = (vertices.size() * (vertices.size() - 1)) / 2;
        double sample_prob = double(num_edge_samples) / total_edges;
        std::geometric_distribution<std::size_t> skip_dist(sample_prob);
        double inv_prob = 1.0 / sample_prob;
        auto& rng = sammy::rng();
        std::vector<double> estimates(vertices.size(), 0.0);
        auto b = vertices.begin(), i = b + 1, j = b, e = vertices.end();
        p_push(*i);
        std::size_t skip = skip_dist(rng);
        for (;;) {
            if (p_advance_pair(i, j, b, e, skip)) {
                while (m_prop.get_current_level() > 0)
                    m_prop.pop_level();
                if (i == e)
                    return estimates;
                p_push(*i);
            }
            if (!p_can_push(*j)) {
                estimates[i - b] += inv_prob;
                estimates[j - b] += inv_prob;
            }
            skip = 1 + skip_dist(rng);
        }
    }

    SharedDBPropagator& propagator() noexcept { return m_prop; }

  private:
    void p_check_prop() {
        if (m_prop.get_current_level() != 0)
            throw std::logic_error("Propagator not at level 0");
        m_prop.incorporate_or_throw();
    }

    template <typename Iterator>
    bool p_advance_pair(Iterator& i, Iterator& j, Iterator begin, Iterator end,
                        std::size_t dist) {
        bool result = false;
        while (i != end && dist > 0) {
            std::size_t row_remaining(i - j);
            if (row_remaining > dist) {
                j += dist;
                break;
            } else {
                dist -= row_remaining;
                ++i;
                j = begin;
                result = true;
            }
        }
        return result;
    }

    std::vector<Vertex> p_compute_clique(Vertex start_vertex) {
        std::vector<Vertex> result;
        p_filter_nonneighbors(start_vertex);
        while (!m_subg.empty()) {
            Vertex v = p_select_next_vertex();
            p_filter_nonneighbors(v);
            result.push_back(v);
        }
        return result;
    }

    Vertex p_select_next_vertex() const {
        std::uniform_int_distribution<std::size_t> dist(0, m_subg.size() - 1);
        return m_subg[dist(sammy::rng())];
    }

    bool p_push(Vertex v) {
        if (m_prop.is_false(v.first) || m_prop.is_false(v.second))
            return false;
        bool pushed = false;
        if (m_prop.is_open(v.first)) {
            if (!m_prop.push_level(v.first) || m_prop.is_false(v.second)) {
                m_prop.pop_level();
                return false;
            }
            pushed = true;
        }
        if (m_prop.is_open(v.second)) {
            if (!m_prop.push_level(v.second)) {
                m_prop.pop_level();
                if (pushed)
                    m_prop.pop_level();
                return false;
            }
        }
        return true;
    }

    bool p_can_push(Vertex v) {
        auto old_level = m_prop.get_current_level();
        if (!p_push(v))
            return false;
        while (m_prop.get_current_level() > old_level)
            m_prop.pop_level();
        return true;
    }

    void p_filter_invalid() {
        m_subg.erase(std::remove_if(m_subg.begin(), m_subg.end(),
                                    [&](Vertex v) { return !p_can_push(v); }),
                     m_subg.end());
    }

    void p_filter_nonneighbors(Vertex vnew) {
        p_push(vnew);
        auto nend = std::remove_if(m_subg.begin(), m_subg.end(), [&](Vertex v) {
            return v == vnew || p_can_push(v);
        });
        m_subg.erase(nend, m_subg.end());
        while (m_prop.get_current_level() > 0)
            m_prop.pop_level();
    }

    SharedDBPropagator m_prop;
    std::vector<Vertex> m_subg;
};

class ParallelFastCliqueBuilder {
  public:
    explicit ParallelFastCliqueBuilder(SharedDBPropagator base_prop,
                                       ThreadGroup<void>* thread_pool)
        : m_base_prop(std::move(base_prop)), m_thread_pool(thread_pool) {}

    std::vector<Vertex>
    random_multistart_best_clique(std::size_t iterations_per_thread,
                                  const std::vector<Vertex>& subgraph) {
        m_base_prop.incorporate_or_throw();
        std::mutex m_out_lock;
        std::vector<Vertex> result;
        m_thread_pool->run_n_copies(m_thread_pool->num_threads() + 1, [&]() {
            FastCliqueBuilder builder{m_base_prop};
            auto rclique = builder.random_multistart_best_clique(
                iterations_per_thread, subgraph);
            {
                std::unique_lock<std::mutex> l{m_out_lock};
                if (rclique.size() > result.size()) {
                    result = std::move(rclique);
                }
            }
        });
        return result;
    }

  private:
    SharedDBPropagator m_base_prop;
    ThreadGroup<void>* m_thread_pool;
};

} // namespace sammy

#endif
==> ./barrage_worker_exact.h <==
#ifndef SAMMY_BARRAGE_WORKER_EXACT_H_INCLUDED_
#define SAMMY_BARRAGE_WORKER_EXACT_H_INCLUDED_

#include "barrage.h"
#include "clique_sat_dsatur.h"
#include "implied_vertices.h"

#include <sammy/cadical_solver.h>

namespace sammy {

class CliqueSatDSaturExactSolverCore {
  public:
    static std::unique_ptr<CliqueSatDSaturExactSolverCore> factory(
        PortfolioSolver* solver,
        PortfolioElementWithCore<CliqueSatDSaturExactSolverCore>* our_element) {
        return std::make_unique<CliqueSatDSaturExactSolverCore>(
            solver, our_element->get_mutable_recorder());
    }

    using IncrementalSolver = CadicalSolver;
    using CSDSSolver = CliqueSatDSaturSolver<IncrementalSolver>;

    explicit CliqueSatDSaturExactSolverCore(PortfolioSolver* solver,
                                            EventRecorder* local_recorder)
        : m_solver(solver), m_clauses(solver->get_clauses()),
          m_local_recorder(*local_recorder), m_icache(p_check_cache()),
          m_exact(m_icache->get_reduced_universe(),
                  &solver->get_infeasibility_map(), solver->get_clauses(),
                  solver->get_best_mes(), solver->get_best_mes().size(),
                  solver->get_best_lower_bound(),
                  solver->get_best_solution().assignments()) {
        m_exact.set_clique_candidate_callback([this]() -> std::vector<Vertex> {
            return m_solver->get_best_mes();
        });
        m_exact.set_lower_bound_callback(
            [this](std::size_t lb, const std::vector<Vertex>& subgraph) {
                if (lb == subgraph.size()) {
                    m_solver->report_mes(subgraph,
                                         "C&P/SATDSatur Exact Solver");
                } else {
                    m_solver->report_lower_bound(lb, subgraph,
                                                 "C&P/SATDSatur Exact Solver");
                }
            });
        m_exact.set_event_recorder(&m_local_recorder);
    }

    /**
     * Aside from handling the termination flag,
     * we do not need to check for other interrupt reasons.
     */
    void interrupt_if_necessary(const InterruptionCheckInfo& /*info*/) {}

    /**
     * Called by the PortfolioElementWithCore
     * if it finds the termination flag to be set.
     */
    void termination_flag_set() { m_exact.abort(); }

    void main() {
        auto status = m_exact.solve();
        if (status == CSDSSolver::SolveResult::ABORTED) {
            m_local_recorder.store_event("EXACT_ABORTED");
            return;
        }
        // the lower bound is reported via the callback;
        // the solution is reported here
        if (status == CSDSSolver::SolveResult::IMPROVED_SOLUTION) {
            auto solution = m_exact.get_partial_solution();
            m_solver->report_solution(solution, "C&P/SATDSatur Exact Solver");
        }
        m_local_recorder.store_event("EXACT_DONE");
    }

  private:
    const ImpliedVertexCache* p_check_cache() {
        const ImpliedVertexCache& icache = m_solver->implied_cache();
        if (!icache.have_reduced_universe()) {
            throw std::logic_error(
                "Did not reduce universe in exact portfolio element!");
        }
        return &icache;
    }

    PortfolioSolver* m_solver;
    ClauseDB& m_clauses;
    EventRecorder& m_local_recorder;
    const ImpliedVertexCache* m_icache;
    CSDSSolver m_exact;
};

} // namespace sammy

#endif
==> ./run_initial.h <==
#ifndef SAMMY_RUN_INITIAL_H_INCLUDED_
#define SAMMY_RUN_INITIAL_H_INCLUDED_

#include "fast_clique.h"
#include "initial_coloring_heuristic.h"
#include "input.h"
#include "learn_infeasibilities.h"
#include "output.h"
#include "simplification.h"
#include <limits>

namespace sammy {

struct RunInitialConfig {
    double min_time = 60.0, max_time = std::numeric_limits<double>::infinity();
    std::size_t random_clique_restarts_per_iteration = 40;
    std::size_t goal_iterations = 20;
    bool simplify = true;
    bool quiet = true;
};

class InitialHeuristicData {
  public:
    InitialHeuristicData(const ClauseDB& simplified, Var n_concrete,
                         EventRecorder* recorder,
                         const RunInitialConfig& config,
                         ThreadGroup<void>* thread_pool)
        : simplified(simplified), inf_map(n_concrete),
          col_solver(&this->simplified, n_concrete, thread_pool, &inf_map),
          clq_solver(SharedDBPropagator(&this->simplified), thread_pool),
          recorder(recorder), config(config) {
        col_solver.set_quiet(config.quiet);
    }

    void first_solution() {
        recorder->store_event("FIRST_SOLVE_BEGINS");
        col_solver.initialize_feasibilities();
        col_solver.color_lazy();
        recorder->store_event("FIRST_COLORING_DONE");
        col_solver.extract_feasibilities();
        recorder->store_event("FEASIBILITIES_EXTRACTED");
        learn_infeasibilities(simplified, &inf_map);
        recorder->store_event("INFEASIBILITIES_LEARNED");
        update_best_ub(0, col_solver.class_spawners());
        recorder->store_event(
            "FIRST_SOLVE_DONE",
            OutputObject{{"num_configs", col_solver.all_classes().size()}});
        all_spawners = col_solver.class_spawners();
        recorder->store_event("FIRST_LB_BEGINS");
        auto r = config.random_clique_restarts_per_iteration;
        update_best_lb(
            0, clq_solver.random_multistart_best_clique(r, all_spawners));
        recorder->store_event("FIRST_LB_DONE");
    }

    void follow_up_iteration(std::size_t iteration) {
        col_solver.reset_coloring();
        col_solver.color_lazy(best_lb);
        update_best_ub(iteration, col_solver.class_spawners());
        auto new_spawners = col_solver.class_spawners();
        add_spawners(new_spawners);
        auto r = config.random_clique_restarts_per_iteration;
        update_best_lb(iteration, clq_solver.random_multistart_best_clique(
                                      r / 2, all_spawners));
        update_best_lb(iteration, clq_solver.random_multistart_best_clique(
                                      r / 2, new_spawners));
    }

    void rerun_primal(const std::vector<Vertex>& initial) {
        auto iteration = std::numeric_limits<std::size_t>::max();
        col_solver.reset_coloring();
        col_solver.color_lazy(initial);
        update_best_ub(iteration, col_solver.class_spawners());
        add_spawners(col_solver.class_spawners());
    }

    void rerun_primal_preset(const std::vector<Vertex>& initial) {
        auto iteration = std::numeric_limits<std::size_t>::max();
        col_solver.color_lazy(initial);
        update_best_ub(iteration, col_solver.class_spawners());
        add_spawners(col_solver.class_spawners());
    }

    void add_spawners(const std::vector<Vertex>& new_spawners) {
        all_spawners.insert(all_spawners.end(), new_spawners.begin(),
                            new_spawners.end());
        std::sort(all_spawners.begin(), all_spawners.end());
        all_spawners.erase(
            std::unique(all_spawners.begin(), all_spawners.end()),
            all_spawners.end());
    }

    void update_best_ub(std::size_t iteration,
                        const std::vector<Vertex>& spawners) {
        update_best_ub(iteration, col_solver.all_classes(), spawners);
    }

    void update_best_ub(std::size_t iteration,
                        const std::vector<SharedDBPropagator>& propagators,
                        const std::vector<Vertex>& spawners) {
        if (best_ub.empty() || propagators.size() < best_ub.size()) {
            best_ub.clear();
            std::transform(propagators.begin(), propagators.end(),
                           std::back_inserter(best_ub),
                           [](const SharedDBPropagator& prop) {
                               return prop.extract_assignment();
                           });
            best_spawners = spawners;
            recorder->store_event(
                "UB_IMPROVED",
                OutputObject{{"num_configs", propagators.size()},
                             {"iteration", iteration + 1}},
                "num_configs", "iteration");
        }
    }

    void update_best_lb(std::size_t iteration, std::vector<Vertex> clique) {
        if (clique.size() > best_lb.size()) {
            recorder->store_event(
                "LB_IMPROVED",
                OutputObject{{"num_interactions", clique.size()},
                             {"iteration", iteration + 1}},
                "num_interactions", "iteration");
            best_lb = std::move(clique);
        }
    }

    std::vector<std::vector<bool>>
    reconstruct_solution(SimplifyDatastructure& simplify_ds,
                         SimplifiedInstance& simplified_inst) {
        return simplify_ds.reconstruct_sample(simplified_inst, best_ub);
    }

    std::vector<Vertex> reconstruct_lb(SimplifyDatastructure& simplify_ds,
                                       SimplifiedInstance& simplified_inst) {
        return simplify_ds.reconstruct_lb(simplified_inst, best_lb);
    }

    const std::vector<std::vector<bool>>& get_solution() const noexcept {
        return best_ub;
    }

    const std::vector<Vertex>& get_bound() const noexcept { return best_lb; }

    std::size_t best_lb_value() const noexcept {
        return (std::max)(std::size_t(1), best_lb.size());
    }

    std::size_t best_ub_value() const noexcept { return best_ub.size(); }

    ColoringHeuristicSolver& primal_solver() noexcept { return col_solver; }

    ParallelFastCliqueBuilder& dual_solver() noexcept { return clq_solver; }

    const ColoringHeuristicSolver& primal_solver() const noexcept {
        return col_solver;
    }

    const ParallelFastCliqueBuilder& dual_solver() const noexcept {
        return clq_solver;
    }

    PairInfeasibilityMap& infeasibility_map() noexcept { return inf_map; }

    const PairInfeasibilityMap& infeasibility_map() const noexcept {
        return inf_map;
    }

    const std::vector<Vertex>& get_all_spawners() const noexcept {
        return all_spawners;
    }

    ClauseDB& formula() noexcept { return simplified; }
    const ClauseDB& formula() const noexcept { return simplified; }
    Var num_concrete() const noexcept { return inf_map.get_n_concrete(); }

    const std::vector<Vertex>& get_best_spawners() const noexcept {
        return best_spawners;
    }

  private:
    ClauseDB simplified;
    PairInfeasibilityMap inf_map;
    ColoringHeuristicSolver col_solver;
    ParallelFastCliqueBuilder clq_solver;
    EventRecorder* recorder;
    std::vector<std::vector<bool>> best_ub;
    std::vector<Vertex> best_spawners;
    std::vector<Vertex> best_lb;
    std::vector<Vertex> all_spawners;
    RunInitialConfig config;
};

inline std::pair<std::size_t, std::size_t> run_initial_heuristic(
    const InputData& input, OutputObject& output, EventRecorder& recorder,
    RunInitialConfig config, ThreadGroup<void>& tgroup,
    std::optional<InitialHeuristicData>& initial_data,
    std::optional<SimplifyDatastructure>& simplifier,
    std::optional<SimplifiedInstance>& simplified,
    std::function<std::vector<std::vector<bool>>()>& get_best_solution_fn,
    std::function<std::vector<Vertex>()>& get_best_mut_ex_set) {
    const ClauseDB* formula = &input.formula;
    Var num_concrete = input.num_concrete;
    if (config.simplify) {
        SimplificationStats simplification_stats;
        simplification_stats.capture_before(input.formula, input.num_concrete);
        simplifier.emplace(input.formula, input.num_concrete);
        simplified.emplace(
            sammy::run_simplifier(*simplifier, simplification_stats));
        num_concrete = simplified->num_concrete;
        formula = &simplified->formula;
        simplification_stats.capture_after(*formula, num_concrete);
        recorder.store_event("SIMPLIFICATION_DONE");
        export_simplification(output, *simplifier, *simplified,
                              &simplification_stats);
    }
    auto before_initial = Clock::now();
    initial_data.emplace(*formula, num_concrete, &recorder, config, &tgroup);
    auto should_stop = [&](std::size_t iteration) -> bool {
        double current_time = seconds_between(before_initial, Clock::now());
        if (iteration >= config.goal_iterations &&
            current_time >= config.min_time)
            return true;
        if (initial_data->best_lb_value() >= initial_data->best_ub_value())
            return true;
        return current_time >= config.max_time;
    };
    if (config.simplify) {
        get_best_solution_fn = [&initial_data, &simplifier, &simplified]() {
            return initial_data->reconstruct_solution(*simplifier, *simplified);
        };
        get_best_mut_ex_set = [&initial_data, &simplifier, &simplified]() {
            return initial_data->reconstruct_lb(*simplifier, *simplified);
        };
    } else {
        get_best_solution_fn = [&initial_data]() {
            return initial_data->get_solution();
        };
        get_best_mut_ex_set = [&initial_data]() {
            return initial_data->get_bound();
        };
    }
    initial_data->first_solution();
    std::size_t actual_iterations;
    for (actual_iterations = 1; !should_stop(actual_iterations);
         ++actual_iterations)
    {
        initial_data->follow_up_iteration(actual_iterations);
    }
    export_solution(output, get_best_solution_fn(), "initial_heuristic",
                    OutputObject{{"iterations", actual_iterations}});
    export_bound(output, get_best_mut_ex_set(), "initial_heuristic",
                 OutputObject{{"random_clique_restarts",
                               config.random_clique_restarts_per_iteration}});
    return {initial_data->best_lb_value(), initial_data->best_ub_value()};
}

inline std::pair<std::size_t, std::size_t>
run_initial_heuristic(const InputData& input, OutputObject& output,
                      EventRecorder& recorder, RunInitialConfig config,
                      ThreadGroup<void>* thread_pool = nullptr) {
    std::optional<ThreadGroup<void>> tg;
    if (!thread_pool) {
        tg.emplace();
        thread_pool = &*tg;
    }
    std::optional<InitialHeuristicData> initial_data;
    std::optional<SimplifyDatastructure> simplifier;
    std::optional<SimplifiedInstance> simplified;
    std::function<std::vector<std::vector<bool>>()> get_best_solution;
    std::function<std::vector<Vertex>()> get_best_mut_ex_set;
    return run_initial_heuristic(input, output, recorder, config, *thread_pool,
                                 initial_data, simplifier, simplified,
                                 get_best_solution, get_best_mut_ex_set);
}

} // namespace sammy

#endif
==> ./cuda_iteration.h <==
#ifndef SAMMY_CUDA_ITERATION_H_INCLUDED_
#define SAMMY_CUDA_ITERATION_H_INCLUDED_

#include "basic_cuda_iteration.h"
#include "dynamic_bitset.h"
#include "literals.h"

namespace sammy {
namespace detail {

template <typename Container,
          std::enable_if_t<sizeof(DynamicBitset::Block) ==
                               0 * sizeof(Container) + sizeof(std::uint64_t),
                           int> = 0>
static void to_prepare_buffer(DynamicBitset::Block b, Container& buffer) {
    buffer.push_back(
        std::uint32_t(b & std::numeric_limits<std::uint32_t>::max()));
    buffer.push_back(std::uint32_t(b >> 32));
}

template <typename Container,
          std::enable_if_t<sizeof(DynamicBitset::Block) ==
                               0 * sizeof(Container) + sizeof(std::uint32_t),
                           int> = 0>
static void to_prepare_buffer(DynamicBitset::Block b, Container& buffer) {
    buffer.push_back(std::uint32_t(b));
}

template <typename T> class CUDADevicePointer {
    static_assert(std::is_pod_v<T>, "T must be POD");

    std::remove_cv_t<T>* m_device_ptr{nullptr};
    std::size_t m_device_size{0};

  public:
    /**
     * Create an empty device pointer.
     */
    CUDADevicePointer() = default;

    /**
     * Create a device pointer pointing to a buffer
     * of size objects of type T.
     */
    CUDADevicePointer(std::size_t size) {
        std::size_t bytes = size * sizeof(T) * (CHAR_BIT / 8);
        cuda_malloc_or_throw((void**)(&m_device_ptr), bytes);
        m_device_size = size;
    }

    /**
     * Create a device pointer pointing to a buffer
     * of size objects of type T, and copy the data
     * from the given vector.
     */
    CUDADevicePointer(const std::vector<std::remove_cv_t<T>>& host_data) {
        std::size_t bytes = host_data.size() * sizeof(T) * (CHAR_BIT / 8);
        cuda_malloc_or_throw((void**)(&m_device_ptr), bytes);
        try {
            cuda_memcpy_htd_or_throw(static_cast<void*>(m_device_ptr),
                                     static_cast<const void*>(host_data.data()),
                                     bytes);
        } catch (...) {
            cuda_free(static_cast<void*>(m_device_ptr));
            m_device_ptr = nullptr;
            throw;
        }
        m_device_size = host_data.size();
    }

    /**
     * Destroy the object (freeing the device memory).
     */
    ~CUDADevicePointer() {
        if (m_device_ptr) {
            cuda_free(static_cast<void*>(m_device_ptr));
        }
    }

    /**
     * Get the pointer to the device memory.
     */
    T* get() const noexcept { return m_device_ptr; }

    /**
     * Copy the device memory to the given host buffer.
     */
    void copy_to_host(std::vector<std::remove_cv_t<T>>& host_buffer) const {
        host_buffer.resize(m_device_size);
        cuda_memcpy_dth_or_throw(host_buffer.data(), m_device_ptr,
                                 m_device_size * sizeof(T) * (CHAR_BIT / 8));
    }

    /**
     * Copy the device memory to a vector and return it.
     */
    std::vector<std::remove_cv_t<T>> get_host_copy() const {
        std::vector<std::remove_cv_t<T>> result;
        copy_to_host(result);
        return result;
    }
};

} // namespace detail

} // namespace sammy

#endif
==> ./lazy_g2_adjacency_matrix.h <==
#ifndef SAMMY_LAZY_G2_ADJACENCY_MATRIX_H_INCLUDED_
#define SAMMY_LAZY_G2_ADJACENCY_MATRIX_H_INCLUDED_

#include "dynamic_bitset.h"
#include "literals.h"
#include "pair_infeasibility_map.h"
#include "shared_db_propagator.h"
#include "thread_interrupt.h"
#include "vertex_operations.h"

namespace sammy {

/**
 * Adjacency matrix of a subgraph of G2 on a
 * fixed subset of all possible vertices; this
 * may either be a reduced vertex set of the
 * full universe or a real subset, e.g., during LNS.
 */
class LazyG2AdjacencyMatrix {
  public:
    using Vertices = std::vector<Vertex>;
    using Indices = std::vector<std::size_t>;
    using VertexMap = VertexMapTo<std::size_t>;

    LazyG2AdjacencyMatrix(std::vector<Vertex> considered_vertices,
                          ClauseDB& clause_db, std::size_t n_concrete);

    /**
     * Get a list of vertex indices containing the given concrete literal l.
     */
    const Indices& vertices_with_concrete_literal(Lit l) const noexcept {
        assert(l < m_vertices_containing_concrete_literal.size());
        return m_vertices_containing_concrete_literal[l];
    }

    /**
     * A vector-based map of:
     *  - concrete literal l -> vertex indices implying l by UP.
     */
    const std::vector<DynamicBitset>&
    vertices_implying_literal() const noexcept {
        return m_definitive_nonedges;
    }

    /**
     * Get a bitset of vertices that imply the given literal l.
     */
    const DynamicBitset& vertices_implying_literal(Lit l) const noexcept {
        assert(l < m_vertices_implying_literals.size());
        return m_vertices_implying_literals[l];
    }

    /**
     * Get a vector-based map of:
     * concrete literal l -> vertex indices containing l.
     */
    const std::vector<Indices>&
    vertices_with_concrete_literal() const noexcept {
        return m_vertices_containing_concrete_literal;
    }

    /**
     * Get a reference to the list of all vertices.
     * Vertex indices point into this list.
     */
    const Vertices& all_vertices() const noexcept { return m_all_vertices; }

    /**
     * Get the number of vertices.
     */
    std::size_t n() const noexcept { return m_all_vertices.size(); }

    /**
     * Get a vertex by its index.
     */
    Vertex vertex(std::size_t index) const noexcept {
        return m_all_vertices[index];
    }

    /**
     * Get a temporary propagator reference.
     * Clear after each use.
     */
    SharedDBPropagator& temp_propagator() noexcept {
        return m_empty_propagator;
    }

    /**
     * Get a reference to the given row of the adjacency matrix.
     */
    DynamicBitset& row(std::size_t index) noexcept {
        return m_adjacency_matrix[index];
    }

    /**
     * Get a const reference to the given row of the adjacency matrix.
     */
    const DynamicBitset& row(std::size_t index) const noexcept {
        return m_adjacency_matrix[index];
    }

    /**
     * Check if there definitely is no edge between
     * vertices with index i and j; if necessary, does
     * the pairwise propagation to confirm this,
     * and cache the result.
     */
    bool is_definitive_nonedge(std::size_t i, std::size_t j) noexcept {
        if (m_adjacency_matrix[i][j])
            return false;
        if (m_definitive_nonedges[i][j])
            return true;
        if (push_vertex_pair(m_empty_propagator, m_all_vertices[i],
                             m_all_vertices[j]))
        {
            m_empty_propagator.reset_to_zero();
            m_definitive_nonedges[i][j].set();
            m_definitive_nonedges[j][i].set();
            return true;
        } else {
            m_adjacency_matrix[i][j].set();
            m_adjacency_matrix[j][i].set();
            return false;
        }
    }

    std::size_t num_vars() const noexcept { return m_n_all; }

    std::size_t num_concrete() const noexcept { return m_n_concrete; }

    void nonedge_to_edge(std::size_t v, std::size_t w) noexcept {
        m_adjacency_matrix[v][w].set();
        m_adjacency_matrix[w][v].set();
        m_definitive_nonedges[v][w].reset();
        m_definitive_nonedges[w][v].reset();
    }

    std::size_t index_of(Vertex v) const noexcept {
        return m_vertex_indices.at(v);
    }

    template <typename Range>
    std::vector<std::size_t> indices_of(Range&& r) const {
        std::vector<std::size_t> result;
        for (Vertex v : r) {
            result.push_back(index_of(v));
        }
        return result;
    }

    template <typename Range> std::vector<Vertex> vertices_of(Range&& r) const {
        std::vector<Vertex> result;
        for (std::size_t i : r) {
            result.push_back(vertex(i));
        }
        return result;
    }

  private:
    void p_initialize_vertices_implying_literal();
    void p_initialize_vertices_containing_concrete_literal();
    void p_initialize_matrix_from_implied_literals();

    /**
     * The number of all literals.
     */
    std::size_t m_n_all;

    /**
     * The number of concrete literals.
     */
    std::size_t m_n_concrete;

    /**
     * A propagator we can reuse.
     * After use, it must be emptied again.
     */
    SharedDBPropagator m_empty_propagator;

    /**
     * The set of all vertices; in other places, they
     * are referenced by their index in this vector.
     */
    std::vector<Vertex> m_all_vertices;

    /**
     * For each concrete literal, the list of vertices containing that literal.
     */
    std::vector<std::vector<std::size_t>>
        m_vertices_containing_concrete_literal;

    /**
     * For each (concrete or non-concrete) literal, the set of vertices
     * _implying_ that literal by propagation.
     */
    std::vector<DynamicBitset> m_vertices_implying_literals;

    /**
     * Adjacency matrix that contains ones for definite edges;
     * zeros might turn out to be edges in G2, but not in the
     * current subgraph.
     */
    std::vector<DynamicBitset> m_adjacency_matrix;

    /**
     * A 'true' means that there is definitely no edge,
     * i.e., propagation of interactions m_all_vertices[i]
     * and m_all_vertices[j] at the same time found no conflict.
     */
    std::vector<DynamicBitset> m_definitive_nonedges;

    /**
     * Get the index of a given vertex by hash-table lookup.
     */
    VertexMap m_vertex_indices;
};

LazyG2AdjacencyMatrix::LazyG2AdjacencyMatrix(
    std::vector<Vertex> considered_vertices, ClauseDB& clause_db,
    std::size_t n_concrete)
    : m_n_all(clause_db.num_vars()), m_n_concrete(n_concrete),
      m_empty_propagator(&clause_db),
      m_all_vertices(std::move(considered_vertices)) {
    m_vertex_indices.reserve(m_all_vertices.size());
    for (std::size_t i = 0, n = m_all_vertices.size(); i < n; ++i) {
        m_vertex_indices[m_all_vertices[i]] = i;
    }
    p_initialize_vertices_containing_concrete_literal();
    p_initialize_vertices_implying_literal();
    p_initialize_matrix_from_implied_literals();
}

/**
 * Initialize the list of vertices containing each concrete literal.
 */
void LazyG2AdjacencyMatrix::
    p_initialize_vertices_containing_concrete_literal() {
    m_vertices_containing_concrete_literal.assign(2 * m_n_concrete, Indices{});
    for (std::size_t vi = 0, vn = m_all_vertices.size(); vi < vn; ++vi) {
        Vertex v = m_all_vertices[vi];
        assert(lit::var(v.first) < m_n_concrete);
        assert(lit::var(v.second) < m_n_concrete);
        m_vertices_containing_concrete_literal[v.first].push_back(vi);
        m_vertices_containing_concrete_literal[v.second].push_back(vi);
    }
}

/**
 * Initialize the vertices_implying_literals member.
 */
void LazyG2AdjacencyMatrix::p_initialize_vertices_implying_literal() {
    m_vertices_implying_literals.reserve(2 * m_n_all);
    std::generate_n(
        std::back_inserter(m_vertices_implying_literals), 2 * m_n_all,
        [&]() { return DynamicBitset(m_all_vertices.size(), false); });
    for (std::size_t vi = 0, n = m_all_vertices.size(); vi < n; ++vi) {
        assert(m_empty_propagator.get_current_level() == 0);
        Vertex v = m_all_vertices[vi];
        if (push_vertex(m_empty_propagator, v) < 0)
            throw std::logic_error("Infeasible interaction in m_all_vertices!");
        for (Lit l : m_empty_propagator.get_trail()) {
            m_vertices_implying_literals[l][vi].set();
        }
        m_empty_propagator.reset_to_zero();
    }
}

void LazyG2AdjacencyMatrix::p_initialize_matrix_from_implied_literals() {
    m_adjacency_matrix.assign(m_all_vertices.size(),
                              DynamicBitset(m_all_vertices.size(), false));
    m_definitive_nonedges.assign(m_all_vertices.size(),
                                 DynamicBitset(m_all_vertices.size(), false));
    assert(m_adjacency_matrix.size() == m_all_vertices.size());
    assert(m_adjacency_matrix.empty() ||
           m_adjacency_matrix[0].size() == m_all_vertices.size());
    for (std::size_t vi = 0, n = m_all_vertices.size(); vi < n; ++vi) {
        assert(m_empty_propagator.get_current_level() == 0);
        Vertex v = m_all_vertices[vi];
        DynamicBitset& row = m_adjacency_matrix[vi];
        push_vertex(m_empty_propagator, v);
        for (Lit lpos : m_empty_propagator.get_trail()) {
            row |= m_vertices_implying_literals[lit::negate(lpos)];
        }
        m_empty_propagator.reset_to_zero();
        if (vi % 16 == 15)
            throw_if_interrupted();
    }
}

} // namespace sammy

#endif
==> ./shared_db_propagator.h <==
#ifndef SAMMY_SHARED_DB_PROPAGATOR_H_INCLUDED_
#define SAMMY_SHARED_DB_PROPAGATOR_H_INCLUDED_

#include "clause_db.h"
#include "error.h"
#include "literals.h"
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
    std::vector<std::pair<std::int32_t, Lit>> supporting_decision_buffer;

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
        while (lbpos < learn_buffer.size()) {
            Lit next = learn_buffer[lbpos++];
            std::size_t tindex = get_trail_index(next);
            if (trail_reasons[tindex].reason_length == 0) {
                supporting_decision_buffer.emplace_back(
                    get_decision_level(next), next);
            } else {
                Reason reason = trail_reasons[tindex];
                for (Lit lr : reason.lits(db())) {
                    if (lr != next) {
                        Var v = lit::var(lr);
                        if (variables[v].get_stamp() != current_stamp) {
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
    template <typename FullAssignment>
    void incorporate_assignment(const FullAssignment& assignment) {
        for (Var v = 0, nv = db().num_vars(); v < nv; ++v) {
            Lit l = assignment[v] ? lit::positive_lit(v) : lit::negative_lit(v);
            if (is_false(l) || (is_open(l) && !push_level(l))) {
                throw std::logic_error("Assignment is infeasible!");
            }
        }
    }

    /**
     * Compute a list of [Level, Literal] pairs of decisions
     * that ultimately led to including l in the trail.
     */
    const std::vector<std::pair<std::int32_t, Lit>>&
    decisions_leading_to(Lit l) {
        if (conflicting)
            throw std::logic_error(
                "decisions_leading_to called on propagator with conflict!");
        if (is_open(l))
            throw std::logic_error(
                "decisions_leading_to called with open literal!");
        supporting_decision_buffer.clear();
        std::size_t tindex = get_trail_index(l);
        if (trail_reasons[tindex].reason_length == 0) {
            supporting_decision_buffer.emplace_back(get_decision_level(l), l);
            return supporting_decision_buffer;
        }

        auto current = p_increase_stamp();
        Reason reason = trail_reasons[tindex];
        for (Lit lr : reason.lits(db())) {
            if (lr != l) {
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
    const std::vector<std::pair<std::int32_t, Lit>>&
    decisions_leading_to_conflict() {
        if (!conflicting)
            throw std::logic_error("decisions_leading_to_conflict called on "
                                   "non-conflicting propagator!");

        supporting_decision_buffer.clear();
        auto current = p_increase_stamp();
        for (Lit lr : conflict_reason.lits(db())) {
            if (lr != conflict_lit) {
                variables[lit::var(lr)].stamp_with(current);
                learn_buffer.push_back(lit::negate(lr));
            }
        }
        variables[lit::var(conflict_lit)].stamp_with(current);
        Lit lc = lit::negate(conflict_lit);
        for (Lit lr : get_reason(lc).lits(db())) {
            if (variables[lit::var(lr)].get_stamp() != current) {
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

    std::vector<Reason>::const_iterator
    current_level_reasons_begin() const noexcept {
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
     * @brief Get the list of reasons for the literals that are currently
     * assigned to true.
     */
    const std::vector<Reason>& get_reasons() const noexcept {
        return trail_reasons;
    }

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
     * @throws UNSATError if this results in a conflict at level 0 (the formula
     * is UNSAT).
     */
    void reset_or_throw() {
        if (is_conflicting())
            resolve_or_throw();
        while (get_current_level() > 0) {
            pop_level();
        }
        incorporate_or_throw();
    }

    /**
     * @brief Resolve conflicts without handler.
     * @throws UNSATError if the conflict cannot be resolved (the formula is
     * UNSAT).
     */
    void resolve_or_throw() {
        if (!resolve_conflicts())
            throw UNSATError();
    }

    /**
     * Reset propagator to level 0.
     */
    void reset_to_zero() noexcept {
        while (get_current_level() > 0) {
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
        if (!incorporate_new_clauses_at_level_0())
            throw UNSATError();
    }

    /**
     * @brief Extract an assignment as bit-vector, where result[i] == true means
     *        that variable i (internal 0-based indexing) is set to true.
     */
    std::vector<bool> extract_assignment() const {
        const Var nv = db().num_vars();
        if (get_trail().size() != nv) {
            throw std::logic_error("Trail incomplete in extract_assignment!");
        }
        std::vector<bool> result(nv, false);
        for (Lit l : get_trail()) {
            if (!lit::negative(l)) {
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

} // namespace sammy

#endif
==> ./subproblem_solver_with_mes.h <==
#ifndef SAMMY_SUBPROBLEM_SOLVER_WITH_MES_H_INCLUDED_
#define SAMMY_SUBPROBLEM_SOLVER_WITH_MES_H_INCLUDED_

#include "global_stat_info.h"
#include "gurobi_clique_solver_g2.h"
#include "incremental_sat_lns.h"
#include <boost/circular_buffer.hpp>

namespace sammy {

constexpr std::size_t MAX_CUT_ROUNDS_WITHOUT_PRICING = 5;
constexpr std::size_t SUB_MES_SOLVER_FULL_GRAPH_THRESHOLD = 200;
constexpr std::size_t PRICING_SAMPLE_SIZE_THRESHOLD = 10'000;

/**
 * Solver for finding good MESs in a subproblem
 * before running another solver.
 */
class SubproblemMESSolver {
  public:
    SubproblemMESSolver(PortfolioSolver* portfolio, LNSSubproblem&& subproblem,
                        SharedDBPropagator prop, EventRecorder* recorder,
                        std::size_t worker_id,
                        std::size_t iterations_without_improvement = 30)
        : m_portfolio(portfolio), m_recorder(recorder), m_worker_id(worker_id),
          m_initial_best_global_mes_size(m_portfolio->get_best_mes_size()),
          m_single_thread(0), m_clauses(&portfolio->get_clauses()),
          m_subproblem(std::move(subproblem)),
          // make & extend subgraph; on interrupt, move subproblem back
          m_subgraph(p_make_universe_subgraph(subproblem)),
          m_base_solver(&m_subgraph, m_clauses->num_vars(), m_recorder,
                        m_subproblem.removed_configurations,
                        m_subproblem.mutually_exclusive_set, p_get_config()),
          m_iterations_without_improvement_limit(
              iterations_without_improvement) {
        m_base_solver.make_abortable();
    }

    void abort() { m_base_solver.abort(); }

    const std::vector<Vertex>& mes_vertices() const {
        return m_base_solver.get_best_mes();
    }

    void return_subproblem(LNSSubproblem&& subproblem) {
        m_subproblem = std::move(subproblem);
    }

    LNSSubproblem move_out_subproblem() { return std::move(m_subproblem); }

    const std::vector<DynamicBitset>& get_solution() const {
        throw std::logic_error(
            "SubproblemMESSolver asked for solution; this should not happen");
    }

    std::size_t num_lp_variables() const {
        return m_base_solver.num_lp_variables();
    }

    std::size_t num_lp_constraints() const {
        return m_base_solver.num_lp_constraints();
    }

    std::size_t num_lp_nonzeros() const {
        return m_base_solver.num_lp_nonzeros();
    }

    std::size_t num_lp_solve_calls() const {
        return m_base_solver.num_lp_solve_calls();
    }

    std::size_t num_graph_vertices() const noexcept {
        return m_base_solver.num_graph_vertices();
    }

    template<typename IterationCallback/*(SubproblemMESSolver&, bool was_interrupted, bool proved_mes_optimality, bool failed)*/>
    std::optional<bool> extended_solve(IterationCallback&& callback) {
        try {
            for (;;) {
                if (get_and_clear_interrupt_flag()) {
                    return std::nullopt;
                }
                auto iteration_result = p_solve_iteration();
                if (!iteration_result) {
                    callback(*this, true, false, false);
                    return std::nullopt;
                }
                if (iteration_result->improved) {
                    if (!p_handle_improvement_should_continue()) {
                        callback(*this, false, true, false);
                        return false;
                    }
                }
                if (iteration_result->optimal_on_subgraph) {
                    if (!p_optimal_on_subgraph_can_continue()) {
                        callback(*this, false, true, false);
                        return true;
                    }
                } else {
                    if (!p_non_optimal_on_subgraph_can_continue()) {
                        callback(*this, false, false, true);
                        return true;
                    }
                }
                if (!callback(*this, false, false, false)) {
                    return true;
                }
            }
        } catch (InterruptError&) {
            return std::nullopt;
        }
    }

    std::optional<bool> solve() {
        try {
            for (;;) {
                if (get_and_clear_interrupt_flag()) {
                    return std::nullopt;
                }
                auto iteration_result = p_solve_iteration();
                if (!iteration_result) {
                    // interrupted
                    return std::nullopt;
                }
                if (iteration_result->improved) {
                    if (!p_handle_improvement_should_continue()) {
                        return false;
                    }
                }
                if (iteration_result->optimal_on_subgraph) {
                    if (!p_optimal_on_subgraph_can_continue()) {
                        return true;
                    }
                } else {
                    if (!p_non_optimal_on_subgraph_can_continue()) {
                        return true;
                    }
                }
                if (!p_should_continue()) {
                    return true;
                }
            }
        } catch (InterruptError&) {
            return std::nullopt;
        }
    }

    /**
     * Whether the last change was a pricing round,
     * i.e., addition of variables.
     */
    bool last_change_was_pricing() const { return m_last_was_pricing; }

    /**
     * Whether the last change was a cutting round,
     * i.e., addition or strengthening of constraints.
     */
    bool last_change_was_cutting() const { return m_last_was_cutting; }

    /**
     * Get the last relaxation value from the clique solver.
     */
    double get_relaxation_value() const noexcept {
        return m_base_solver.get_last_value();
    }

    /**
     * Get the integral upper bound on the current subgraph.
     */
    std::size_t get_subgraph_ub() const noexcept {
        return m_base_solver.get_subgraph_ub();
    }

    /**
     * Get the number of non-trivial vertex values in the current relaxation.
     */
    std::size_t fractional_mes_support_size() const noexcept {
        return m_base_solver.get_fractional_mes_support().size();
    }

  private:
    /**
     * Check if we should continue the LNS
     * subproblem solver MES computation.
     */
    bool p_should_continue() {
        if (m_iterations_without_improvement >=
            m_iterations_without_improvement_limit)
        {
            m_recorder->store_event(
                "LNS_MES_STOPPING_DUE_TO_MISSING_IMPROVEMENT",
                {{"worker_id", m_worker_id}}, "worker_id");
            return false;
        }
        return true;
    }

    /**
     * Price a sample of vertices from the universe.
     */
    bool p_price_sample() {
        if (m_subproblem.uncovered_universe.size() <=
            PRICING_SAMPLE_SIZE_THRESHOLD)
        {
            return p_price_all();
        }
        std::size_t sample_goal_size = 1000;
        sample_goal_size = (std::max)(
            sample_goal_size,
            std::size_t(0.1 * m_subproblem.uncovered_universe.size()));
        if (sample_goal_size >= m_subproblem.uncovered_universe.size()) {
            return p_price_all();
        }
        auto sample = sample_from_range(m_subproblem.uncovered_universe,
                                        sample_goal_size, sammy::rng());
        return m_base_solver.price_vertices(sample.begin(), sample.end()) > 0;
    }

    /**
     * Price all vertices of the universe.
     */
    bool p_price_all() {
        return m_base_solver.price_vertices(
                   m_subproblem.uncovered_universe.begin(),
                   m_subproblem.uncovered_universe.end()) > 0;
    }

    /**
     * Price a sample or, if necessary, all vertices of the universe.
     */
    bool p_price_sample_or_all() {
        if (m_subproblem.uncovered_universe.size() <=
            PRICING_SAMPLE_SIZE_THRESHOLD)
        {
            return p_price_all();
        }
        return p_price_sample() || p_price_all();
    }

    /**
     * Handle the case when we know that we have
     * the optimum MES on the current subgraph.
     * Here, we have to begin by pricing.
     */
    bool p_optimal_on_subgraph_can_continue() {
        m_cut_rounds_without_pricing = 0;
        if (p_price_sample_or_all()) {
            m_last_was_cutting = false;
            m_last_was_pricing = true;
            return true;
        }
        return false;
    }

    /**
     * Handle the case when we don't know if we have the
     * optimum MES on the current subgraph.
     * Here, we start by cuts; if separation fails or
     * every few rounds, we also perform pricing.
     */
    bool p_non_optimal_on_subgraph_can_continue() {
        if (++m_cut_rounds_without_pricing == 6) {
            m_cut_rounds_without_pricing = 0;
            if (p_price_sample()) {
                m_last_was_pricing = true;
                m_last_was_cutting = false;
                return true;
            }
        }
        if (!p_cheap_cut_round()) {
            return p_cheap_cuts_failed();
        }
        m_last_was_pricing = false;
        m_last_was_cutting = true;
        return true;
    }

    /**
     * Run a cheap cutting plane round;
     * return true if a cut was found or
     * an existing constraint strengthened
     * to cut of the current relaxed solution.
     */
    bool p_cheap_cut_round() {
        return m_base_solver.greedy_add_to_cutting_planes() ||
               m_base_solver.greedy_generate_cutting_planes();
    }

    bool p_cheap_cuts_failed() {
        // TODO expensive cuts?
        m_cut_rounds_without_pricing = 0;
        if (p_price_sample_or_all()) {
            m_last_was_pricing = true;
            m_last_was_cutting = false;
            return true;
        }
        return false;
    }

    bool p_handle_improvement_should_continue() {
        std::size_t new_mes_size = m_base_solver.get_best_mes().size();
        m_recorder->store_event(
            "LNS_MES_IMPROVED_BY_CNP",
            {{"mes_size", new_mes_size}, {"worker_id", m_worker_id}},
            "mes_size", "worker_id");
        if (new_mes_size > m_initial_best_global_mes_size) {
            m_portfolio->report_mes(m_base_solver.get_best_mes(),
                                    "Cut & Price on LNS subproblem");
            m_initial_best_global_mes_size = new_mes_size;
        }
        if (new_mes_size >= m_subproblem.removed_configurations.size()) {
            m_recorder->store_event("LNS_MES_BY_CNP_PRECLUDES_IMPROVEMENT",
                                    {{"worker_id", m_worker_id}}, "worker_id");
            return false;
        }
        return true;
    }

    struct IterationResult {
        bool improved;
        bool optimal_on_subgraph;
    };

    std::optional<IterationResult> p_solve_iteration() {
        std::size_t old_mes_size = m_base_solver.get_best_mes().size();
        IterationResult result{false, false};
        ++m_iterations_without_improvement;
        auto before_solve_full = std::chrono::steady_clock::now();
        auto sfr_result = m_base_solver.solve_full_relaxation();
        auto after_solve_full = std::chrono::steady_clock::now();
        get_global_stats().double_stat_add(
            "MESImprovementSolveFullRelaxationTime",
            seconds_between(before_solve_full, after_solve_full));
        switch (sfr_result) {
        default:
        case SolverState::TIMEOUT_IMPROVEMENT:
        case SolverState::TIMEOUT_NO_IMPROVEMENT:
            return std::nullopt;

        case SolverState::OPTIMUM_ON_SUBGRAPH:
            result.optimal_on_subgraph = true;
            if (old_mes_size < m_base_solver.get_best_mes().size()) {
                result.improved = true;
                m_iterations_without_improvement = 0;
            }
            return result;

        case SolverState::IMPROVEMENT_FOUND:
            result.improved = true;
            m_iterations_without_improvement = 0;
            return result;

        case SolverState::NO_IMPROVEMENT_FOUND:
            return result;
        }
    }

    std::vector<Vertex> p_get_initial_vertices() {
        if (m_subproblem.uncovered_universe.size() <=
            SUB_MES_SOLVER_FULL_GRAPH_THRESHOLD)
        {
            return m_subproblem.uncovered_universe;
        }
        return p_get_initial_vertices_sample(
            SUB_MES_SOLVER_FULL_GRAPH_THRESHOLD);
    }

    std::vector<Vertex> p_get_initial_vertices_sample(std::size_t sample_size) {
        return sample_from_range(m_subproblem.uncovered_universe, sample_size,
                                 sammy::rng());
    }

    static LowerBoundMIPConfig p_get_config() {
        LowerBoundMIPConfig config;
        config.quiet_gurobi = true;
        return config;
    }

    UniverseSubgraph
    p_make_universe_subgraph(LNSSubproblem& restore_subproblem_to) {
        try {
            std::vector<Vertex> initial_vertices = p_get_initial_vertices();
            UniverseSubgraph subgraph{m_clauses, &m_single_thread,
                                      &m_portfolio->get_infeasibility_map(),
                                      std::move(initial_vertices)};
            subgraph.extend_matrix_by_propagation(true);
            return subgraph;
        } catch (InterruptError&) {
            restore_subproblem_to = std::move(m_subproblem);
            throw;
        }
    }

    PortfolioSolver* m_portfolio;
    EventRecorder* m_recorder;
    std::size_t m_worker_id;
    std::size_t m_initial_best_global_mes_size;
    ThreadGroup<void> m_single_thread;
    ClauseDB* m_clauses;
    LNSSubproblem m_subproblem;
    UniverseSubgraph m_subgraph;
    GurobiCliqueSolverG2 m_base_solver;
    std::size_t m_cut_rounds_without_pricing = 0;
    std::size_t m_iterations_without_improvement = 0;
    std::size_t m_iterations_without_improvement_limit;
    bool m_last_was_cutting{false};
    bool m_last_was_pricing{false};
};

template <typename WrappedSubproblemSolver> class SubproblemSolverWithMES {
  public:
    static std::string name() {
        return std::string("MESImprovement<") +
               WrappedSubproblemSolver::name() + ">";
    }

    SubproblemSolverWithMES(PortfolioSolver* portfolio,
                            LNSSubproblem&& subproblem, SharedDBPropagator prop,
                            EventRecorder* recorder, std::size_t worker_id)
        : m_portfolio(portfolio), m_recorder(recorder), m_worker_id(worker_id),
          m_time_before_construction(std::chrono::steady_clock::now()),
          m_mes_solver(std::make_unique<SubproblemMESSolver>(
              portfolio, std::move(subproblem), std::move(prop), recorder,
              worker_id)),
          m_sub_solver() {
        get_global_stats().double_stat_add(
            "MESImprovementConstructionTime",
            seconds_between(m_time_before_construction,
                            std::chrono::steady_clock::now()));
    }

    void abort() {
        // we only need the lock when:
        //  - checking the pointers from control thread
        //  - changing the pointers in the worker thread
        std::unique_lock l{m_abort_mutex};
        if (m_mes_solver) {
            m_mes_solver->abort();
        } else {
            m_sub_solver->abort();
        }
    }

    std::optional<bool> solve() {
        auto time_before = std::chrono::steady_clock::now();
        std::optional<bool> mes_result = m_mes_solver->solve();
        auto time_after = std::chrono::steady_clock::now();
        get_global_stats().double_stat_add(
            "MESImprovementTime", seconds_between(time_before, time_after));
        if (!mes_result || !*mes_result) {
            get_global_stats().int_stat_add("MESImprovementDidResolveLNS", 1);
            return mes_result;
        } else {
            get_global_stats().int_stat_add("MESImprovementDidNotResolveLNS",
                                            1);
        }
        LNSSubproblem subproblem = m_mes_solver->move_out_subproblem();
        m_original_mes = std::move(subproblem.mutually_exclusive_set);
        subproblem.mutually_exclusive_set = m_mes_solver->mes_vertices();
        try {
            p_switch_phase(std::move(subproblem));
        } catch (InterruptError&) {
            subproblem.mutually_exclusive_set = std::move(m_original_mes);
            m_mes_solver->return_subproblem(std::move(subproblem));
            return std::nullopt;
        }
        auto time_after_phase_switch = std::chrono::steady_clock::now();
        get_global_stats().double_stat_add(
            "CoreSolverConstructionTime",
            seconds_between(time_after, time_after_phase_switch));
        auto solver_result = m_sub_solver->solve();
        get_global_stats().double_stat_add(
            "CoreSolverSolveTime",
            seconds_between(time_after_phase_switch,
                            std::chrono::steady_clock::now()));
        return solver_result;
    }

    LNSSubproblem move_out_subproblem() noexcept {
        if (m_mes_solver) {
            return m_mes_solver->move_out_subproblem();
        } else {
            LNSSubproblem subproblem = m_sub_solver->move_out_subproblem();
            subproblem.mutually_exclusive_set = std::move(m_original_mes);
            return subproblem;
        }
    }

    const std::vector<Vertex>& mes_vertices() const {
        if (m_mes_solver) {
            return m_mes_solver->mes_vertices();
        } else {
            return m_sub_solver->mes_vertices();
        }
    }

    const std::vector<DynamicBitset>& get_solution() const {
        if (m_mes_solver) {
            return m_mes_solver->get_solution();
        } else {
            return m_sub_solver->get_solution();
        }
    }

  private:
    void p_switch_phase(LNSSubproblem&& subproblem) {
        assert(m_mes_solver);
        assert(!m_sub_solver);
        auto new_sub_solver = std::make_unique<WrappedSubproblemSolver>(
            m_portfolio, std::move(subproblem),
            SharedDBPropagator(&m_portfolio->get_clauses()), m_recorder,
            m_worker_id);
        std::unique_ptr<SubproblemMESSolver>
            tmp; // let destructor run outside of locked region
        {
            std::unique_lock l{m_abort_mutex};
            m_sub_solver = std::move(new_sub_solver);
            tmp = std::move(m_mes_solver);
        }
        assert(m_sub_solver);
        assert(!m_mes_solver);
    }

    PortfolioSolver* m_portfolio;
    EventRecorder* m_recorder;
    std::size_t m_worker_id;
    std::chrono::steady_clock::time_point m_time_before_construction;

    std::mutex m_abort_mutex;
    std::unique_ptr<SubproblemMESSolver> m_mes_solver;
    std::vector<Vertex> m_original_mes;
    std::unique_ptr<WrappedSubproblemSolver> m_sub_solver;
};

} // namespace sammy

#endif
==> ./learn_infeasibilities.h <==
#ifndef SAMMY_LEARN_INFEASIBILITIES_H_INCLUDED_
#define SAMMY_LEARN_INFEASIBILITIES_H_INCLUDED_

#include "literals.h"
#include "pair_infeasibility_map.h"
#include "shared_db_propagator.h"

namespace sammy {

namespace detail {

static inline void
learn_unary_infeasibilities(SharedDBPropagator& propagator,
                            const PairInfeasibilityMap* infeasibilities) {
    const auto& inf = *infeasibilities;
    const Lit nclit = 2 * inf.get_n_concrete();
    for (Lit l = 0; l < nclit; ++l) {
        if (inf[l][l]) {
            if (!propagator.is_false(l)) {
                Lit new_unary[1] = {lit::negate(l)};
                propagator.db().add_clause(+new_unary, new_unary + 1);
                propagator.incorporate_or_throw();
            }
        }
    }
}

} // namespace detail

/**
 * @brief Learn clauses so that all infeasible pairs
 *        are recognized easily by propagation, no matter
 *        which literal of the pair is pushed.
 *        Only learns clauses as they're needed, does
 *        not blindly create a new binary clause for each
 *        literal pair. Similarly handles missing unary clauses.
 *
 * @param clauses
 * @param infeasibilities
 */
inline void learn_infeasibilities(ClauseDB& clauses,
                                  const PairInfeasibilityMap* infeasibilities) {
    const Lit nclit = 2 * infeasibilities->get_n_concrete();
    const auto& inf = *infeasibilities;
    SharedDBPropagator propagator(&clauses);
    detail::learn_unary_infeasibilities(propagator, infeasibilities);
    for (Lit l = 0; l < nclit; ++l) {
        // the literal is false at level 0
        if (propagator.is_false(l)) {
            if (!inf[l][l]) {
                throw std::logic_error(
                    "Non-infeasible literal false at level 0!");
            }
            continue;
        }
        // the literal is true or open at level 0
        bool pushed_outer = false;
        if (propagator.is_open(l)) {
            if (!propagator.push_level(l)) {
                throw std::logic_error(
                    "Non-infeasible literal conflicts at level 1!");
            }
            pushed_outer = true;
        }
        // we are either at level 0 or 1, and l is now set to true
        for (Lit l2 : inf[l].ones()) {
            if (!propagator.is_false(l2)) {
                Lit new_binary[2] = {lit::negate(l), lit::negate(l2)};
                clauses.add_clause(+new_binary, new_binary + 2);
                if (pushed_outer)
                    propagator.pop_level();
                propagator.incorporate_or_throw();
                pushed_outer = false;
                if (propagator.is_open(l)) {
                    if (!propagator.push_level(l)) {
                        throw std::logic_error(
                            "Non-infeasible literal conflicts at level 1!");
                    }
                    pushed_outer = true;
                } else if (propagator.is_false(l)) {
                    throw std::logic_error(
                        "Non-infeasible literal false at level 0!");
                }
                assert(propagator.is_false(l2));
            }
        }
        if (pushed_outer)
            propagator.pop_level();
        assert(propagator.get_current_level() == 0);
    }
}

} // namespace sammy

#endif
==> ./partial_solution.h <==
#ifndef SAMMY_PARTIAL_SOLUTION_H_INCLUDED_
#define SAMMY_PARTIAL_SOLUTION_H_INCLUDED_

#include "dynamic_bitset.h"
#include "literals.h"
#include "pair_infeasibility_map.h"
#include "thread_interrupt.h"

namespace sammy {

class PartialSolution {
  public:
    explicit PartialSolution(Var num_vars, const PairInfeasibilityMap* inf_map)
        : PartialSolution(num_vars, inf_map,
                          static_cast<const DynamicBitset*>(nullptr),
                          static_cast<const DynamicBitset*>(nullptr)) {}

    template <typename ClassIterator>
    explicit PartialSolution(Var num_vars, const PairInfeasibilityMap* inf_map,
                             ClassIterator classes_begin,
                             ClassIterator classes_end)
        : m_assignments_with_literal(2 * inf_map->get_n_concrete()),
          m_inf_map(inf_map), m_num_vars(num_vars),
          m_row_uncolored(2 * inf_map->get_n_concrete(), true) {
        using ClassType = std::remove_reference_t<std::remove_cv_t<
            typename std::iterator_traits<ClassIterator>::value_type>>;

        m_assignments.reserve(std::distance(classes_begin, classes_end));
        m_true_concrete_literals.reserve(
            std::distance(classes_begin, classes_end));
        if constexpr (std::is_same_v<ClassType, SharedDBPropagator>) {
            p_from_props(classes_begin, classes_end);
        } else {
            p_from_bitsets(classes_begin, classes_end);
        }
    }

    template <typename Callable /*(Lit lmin, Lit lmax)*/>
    void iterate_all_uncovered(Callable&& callable,
                               bool interruptible = false) {
        Lit nclit = 2 * m_inf_map->get_n_concrete();
        std::size_t count = 0;
        for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
            m_row_uncolored.set();
            m_row_uncolored ^= (*m_inf_map)[lmin];
            const auto& with_lit = m_assignments_with_literal[lmin];
            for (Index i : with_lit) {
                m_row_uncolored.binary_subtract(m_true_concrete_literals[i],
                                                lmin + 1);
            }
            for (Lit lmax : m_row_uncolored.ones_from(lmin + 1)) {
                (std::forward<Callable>(callable))(lmin, lmax);
            }
            if (interruptible && ++count == 512) {
                count = 0;
                throw_if_interrupted();
            }
        }
    }

    std::size_t get_n_concrete() const noexcept {
        return m_inf_map->get_n_concrete();
    }

    std::size_t count_uncovered() noexcept {
        std::size_t result = 0;
        Lit nclit = 2 * m_inf_map->get_n_concrete();
        for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
            m_row_uncolored.set();
            m_row_uncolored ^= (*m_inf_map)[lmin];
            const auto& with_lit = m_assignments_with_literal[lmin];
            for (Index i : with_lit) {
                m_row_uncolored.binary_subtract(m_true_concrete_literals[i],
                                                lmin + 1);
            }
            result += m_row_uncolored.count_from(lmin + 1);
        }
        return result;
    }

    template <typename Callable /*(Lit lmin, Lit lmax)*/>
    void iterate_all_uniquely_covered(Callable&& callable,
                                      bool cover_is_incomplete = false) {
        if (cover_is_incomplete) {
            p_iterate_uniquely_covered_incomplete(
                std::forward<Callable>(callable));
        } else {
            p_iterate_uniquely_covered_complete(
                std::forward<Callable>(callable));
        }
    }

    template <typename Callable /*(Index class_index, Lit lmin, Lit lmax)*/>
    void iterate_all_uniquely_covered_with_class(Callable&& callable) {
        Lit nclit = 2 * m_inf_map->get_n_concrete();
        // track everything that is multiply-covered
        DynamicBitset multi_covered = make_bitset(nclit, false);
        // track everything that is covered or infeasible
        DynamicBitset covered = make_bitset(nclit, false);
        // temporary buffer
        DynamicBitset tmp = make_bitset(nclit, false);
        for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
            const auto& with_lit = m_assignments_with_literal[lmin];
            covered = (*m_inf_map)[lmin];
            multi_covered = covered;
            for (Index cci : with_lit) {
                const auto& ccbits = m_true_concrete_literals[cci];
                tmp = covered;
                tmp &= ccbits; // tmp = set of previously covered that are by
                               // config cci
                multi_covered |= tmp;
                covered |= ccbits;
            }
            multi_covered.flip(); // multi_covered = everything that's not
                                  // covered multiple times
            multi_covered &=
                covered; // multi_covered = everything that's not covered
                         // multiple times and not uncovered
            for (Index cci : with_lit) {
                tmp = m_true_concrete_literals[cci];
                tmp &= multi_covered;
                for (Lit lmax : tmp.ones_from(lmin + 1)) {
                    (std::forward<Callable>(callable))(cci, lmin, lmax);
                }
            }
        }
    }

    Index add_class(const SharedDBPropagator& prop) {
        Index i = m_assignments.size();
        p_add_class_from_prop(i, prop);
        return i;
    }

    template <typename BitsetType>
    Index add_assignment(BitsetType&& assignment) {
        Index i = m_assignments.size();
        p_add_class_from_bitset(i, std::forward<BitsetType>(assignment));
        return i;
    }

    std::vector<DynamicBitset>
    remove_assignments(const std::vector<Index>& indices) {
        std::vector<DynamicBitset> result;
        result.reserve(indices.size());
        DynamicBitset is_deleted(m_assignments.size(), false);
        for (Index i : indices)
            is_deleted[i].set();
        std::vector<Index> old_to_new(m_assignments.size(), NIL);
        for (std::size_t i = 0, o = 0, s = m_assignments.size(); i < s; ++i) {
            if (!is_deleted[i]) {
                old_to_new[i] = o++;
            } else {
                result.emplace_back(std::move(m_assignments[i]));
                m_true_concrete_literals[i].clear();
            }
        }
        m_assignments.erase(std::remove_if(m_assignments.begin(),
                                           m_assignments.end(),
                                           [](const DynamicBitset& b) {
                                               return b.size() == 0;
                                           }),
                            m_assignments.end());
        m_true_concrete_literals.erase(
            std::remove_if(
                m_true_concrete_literals.begin(),
                m_true_concrete_literals.end(),
                [](const DynamicBitset& b) { return b.size() == 0; }),
            m_true_concrete_literals.end());
        std::for_each(
            m_assignments_with_literal.begin(),
            m_assignments_with_literal.end(), [&](std::vector<Index>& sinds) {
                std::transform(sinds.begin(), sinds.end(), sinds.begin(),
                               [&](Index x) { return old_to_new[x]; });
                sinds.erase(std::remove(sinds.begin(), sinds.end(), NIL),
                            sinds.end());
            });
        return result;
    }

    std::size_t size() const noexcept { return m_assignments.size(); }

    bool empty() const noexcept { return m_assignments.empty(); }

    const DynamicBitset&
    get_assignment(std::size_t class_index) const noexcept {
        return m_assignments[class_index];
    }

    const DynamicBitset&
    get_true_concrete_literals(std::size_t class_index) const noexcept {
        return m_true_concrete_literals[class_index];
    }

    const std::vector<DynamicBitset>& assignments() const noexcept {
        return m_assignments;
    }

    template <typename BitsetType>
    std::vector<BitsetType> assignments_as() const {
        std::vector<BitsetType> result;
        for (const DynamicBitset& assignment : m_assignments) {
            BitsetType a(m_num_vars, false);
            for (Var v : range(m_num_vars)) {
                if (assignment[v])
                    a[v] = true;
            }
            result.emplace_back(std::move(a));
        }
        return result;
    }

    Index find_covering_class(Lit lmin, Lit lmax) const {
        if ((*m_inf_map)(lmin, lmax))
            return NIL;
        const auto& l1 = m_assignments_with_literal[lmin];
        const auto& l2 = m_assignments_with_literal[lmax];
        const auto& list = l1.size() < l2.size() ? l1 : l2;
        const Lit find_lit = l1.size() < l2.size() ? lmax : lmin;
        for (Index cls_index : list) {
            if (m_true_concrete_literals[cls_index][find_lit])
                return cls_index;
        }
        throw std::out_of_range(
            "Vertex (" + std::to_string(lmin) + ", " + std::to_string(lmax) +
            ") is not covered and not marked as infeasible!");
    }

    template <typename Callable /*(Index)*/>
    void find_covering_classes(Lit lmin, Lit lmax, Callable&& callable) const {
        if ((*m_inf_map)(lmin, lmax))
            return;
        const auto& l1 = m_assignments_with_literal[lmin];
        const auto& l2 = m_assignments_with_literal[lmax];
        const auto& list = l1.size() < l2.size() ? l1 : l2;
        const Lit find_lit = l1.size() < l2.size() ? lmax : lmin;
        for (Index cls_index : list) {
            if (m_true_concrete_literals[cls_index][find_lit]) {
                std::forward<Callable>(callable)(cls_index);
            }
        }
    }

    std::size_t bytes_used() const noexcept {
        std::size_t result = sizeof(PartialSolution);
        result += m_assignments.capacity() * sizeof(DynamicBitset);
        if (!m_assignments.empty()) {
            result += m_assignments[0].bytes_used() * m_assignments.size();
        }
        result += m_true_concrete_literals.capacity() * sizeof(DynamicBitset);
        if (!m_true_concrete_literals.empty()) {
            result += m_true_concrete_literals[0].bytes_used() *
                      m_true_concrete_literals.size();
        }
        result +=
            sizeof(std::vector<Index>) * m_assignments_with_literal.capacity();
        for (const auto& v : m_assignments_with_literal) {
            result += sizeof(Index) * v.capacity();
        }
        return result * (CHAR_BIT / 8);
    }

  private:
    template <typename Callable>
    void p_iterate_uniquely_covered_complete(Callable&& callable) {
        Lit nclit = 2 * m_inf_map->get_n_concrete();
        DynamicBitset multi_covered = make_bitset(nclit, false);
        DynamicBitset prev_covered = make_bitset(nclit, false);
        DynamicBitset tmp = make_bitset(nclit, false);
        for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
            multi_covered.reset();
            prev_covered.reset();
            const auto& with_lit = m_assignments_with_literal[lmin];
            for (Index cci : with_lit) {
                const auto& ccbits = m_true_concrete_literals[cci];
                tmp = prev_covered;
                tmp &= ccbits;
                multi_covered |= tmp;
                prev_covered |= ccbits;
            }
            tmp = multi_covered;
            tmp.flip();
            tmp -= (*m_inf_map)[lmin];
            for (Lit lmax : tmp.ones_from(lmin + 1)) {
                (std::forward<Callable>(callable))(lmin, lmax);
            }
        }
    }

    template <typename Callable>
    void p_iterate_uniquely_covered_incomplete(Callable&& callable) {
        Lit nclit = 2 * m_inf_map->get_n_concrete();
        // track everything that is multiply-covered
        DynamicBitset multi_covered = make_bitset(nclit, false);
        // track everything that is covered or infeasible
        DynamicBitset covered = make_bitset(nclit, false);
        // temporary buffer
        DynamicBitset tmp = make_bitset(nclit, false);
        for (Lit lmin = 0; lmin < nclit - 2; ++lmin) {
            const auto& with_lit = m_assignments_with_literal[lmin];
            covered = (*m_inf_map)[lmin];
            multi_covered = covered;
            for (Index cci : with_lit) {
                const auto& ccbits = m_true_concrete_literals[cci];
                tmp = covered;
                tmp &= ccbits; // tmp = set of previously covered that are by
                               // config cci
                multi_covered |= tmp;
                covered |= ccbits;
            }
            multi_covered.flip(); // multi_covered = everything that's not
                                  // covered multiple times
            multi_covered &=
                covered; // multi_covered = everything that's not covered
                         // multiple times and not uncovered
            for (Lit lmax : multi_covered.ones_from(lmin + 1)) {
                (std::forward<Callable>(callable))(lmin, lmax);
            }
        }
    }

    template <typename Iterator>
    void p_from_props(Iterator begin, Iterator end) {
        for (Index iclass = 0; begin != end; ++begin, ++iclass) {
            const auto& propagator = *begin;
            p_add_class_from_prop(iclass, propagator);
        }
    }

    void p_add_class_from_prop(Index class_index,
                               const SharedDBPropagator& prop) {
        Var n_conc = m_inf_map->get_n_concrete();
        DynamicBitset true_concrete(2 * n_conc, false);
        DynamicBitset assignment(m_num_vars, false);
        for (Var v : range(n_conc)) {
            bool vtrue = prop.is_true(lit::positive_lit(v));
            Lit t = vtrue ? lit::positive_lit(v) : lit::negative_lit(v);
            true_concrete[t].set();
            m_assignments_with_literal[t].push_back(class_index);
            assignment[v] = vtrue;
        }
        for (Var v : range(n_conc, m_num_vars)) {
            bool vtrue = prop.is_true(lit::positive_lit(v));
            assignment[v] = vtrue;
        }
        m_assignments.emplace_back(std::move(assignment));
        m_true_concrete_literals.emplace_back(std::move(true_concrete));
    }

    template <typename BitsetType>
    void p_add_class_from_bitset(Index class_index, BitsetType&& assignment) {
        Var n_conc = m_inf_map->get_n_concrete();
        DynamicBitset true_concrete(2 * n_conc, false);
        for (Var v : range(n_conc)) {
            Lit t = assignment[v] ? lit::positive_lit(v) : lit::negative_lit(v);
            true_concrete[t].set();
            m_assignments_with_literal[t].push_back(class_index);
        }
        m_assignments.emplace_back(std::forward<BitsetType>(assignment));
        m_true_concrete_literals.emplace_back(std::move(true_concrete));
    }

    template <typename Iterator>
    void p_from_bitsets(Iterator begin, Iterator end) {
        std::size_t index = 0;
        for (const auto& c : IteratorRange{begin, end}) {
            p_add_class_from_bitset(index++, c);
        }
    }

    std::vector<DynamicBitset> m_assignments;
    std::vector<DynamicBitset> m_true_concrete_literals;
    std::vector<std::vector<Index>> m_assignments_with_literal;
    const PairInfeasibilityMap* m_inf_map;
    Var m_num_vars;
    DynamicBitset m_row_uncolored;
};

} // namespace sammy

#endif
==> ./barrage_worker_lns.h <==
#ifndef SAMMY_BARRAGE_WORKER_LNS_H_INCLUDED_
#define SAMMY_BARRAGE_WORKER_LNS_H_INCLUDED_

#include "algorithm_ex.h"
#include "barrage.h"
#include "barrage_worker_with_core.h"
#include "clique_sat_dsatur.h"
#include "gurobi_clique_solver_g2.h"
#include "implied_vertices.h"
#include "lns_destroy.h"
#include "satdsatur_colorer.h"
#include "thread_interrupt.h"

namespace sammy {

static constexpr double RANDOM_DESTRUCTION_PROB = 0.3;
static constexpr double RANDOM_WITH_TABLE_PROB = 0.2;
static constexpr double RANDOM_AND_LEAST_UNIQUE_DESTRUCTION_PROB = 0.5;

template <typename SubproblemSolverType> class SubproblemLNSSolverCore {
  public:
    using SubproblemSolver = SubproblemSolverType;

    static std::unique_ptr<SubproblemLNSSolverCore>
    factory(PortfolioSolver* solver,
            PortfolioElementWithCore<SubproblemLNSSolverCore>* element) {
        return std::make_unique<SubproblemLNSSolverCore>(solver, element);
    }

    SubproblemLNSSolverCore(
        PortfolioSolver* solver,
        PortfolioElementWithCore<SubproblemLNSSolverCore>* element)
        : m_portfolio(solver), m_element_mutex(&element->get_mutex()),
          m_termination_flag(get_interrupt_flag_ptr()), m_interrupt_flag(false),
          m_local_recorder(element->get_mutable_recorder()),
          m_source("LNSCore<" + SubproblemSolverType::name() + ">"),
          m_clauses(&solver->get_clauses()), m_propagator(m_clauses),
          m_destroy(m_local_recorder, m_worker_counter++,
                    &solver->implied_cache(), solver->get_best_solution(),
                    solver->get_best_mes(), m_propagator),
          m_destroy_seen_solution(m_destroy.total_num_configs()) {
        set_interrupt_flag_ptr(&m_interrupt_flag);
        if (m_termination_flag->load()) {
            throw InterruptError();
        }
    }

    /**
     * Called by the PortfolioElementWithCore
     * if it finds the termination flag to be set.
     */
    void termination_flag_set() {
        assert(m_termination_flag->load());
        p_interrupt();
    }

    void interrupt_if_necessary(const InterruptionCheckInfo& info) {
        if (m_termination_flag->load()) {
            p_interrupt();
        } else if (info.best_upper_bound < m_destroy_seen_solution) {
            m_destroy_seen_solution = info.best_upper_bound;
            p_interrupt();
        }
    }

    void main() {
        while (!m_termination_flag->load()) {
            auto phase_result = p_destroy_phase();
            if (phase_result == DestroyPhaseResult::GLOBAL_OPT) {
                break;
            } else if (phase_result == DestroyPhaseResult::INTERRUPTED) {
                continue;
            }
            LNSSubproblem subproblem = m_destroy.move_out_subproblem();
            if (m_portfolio->subproblem_reporting_enabled()) {
                m_portfolio->report_subproblem(
                    subproblem, m_destroy.get_partial(),
                    m_destroy.best_global_mes().size(),
                    m_portfolio->get_best_lower_bound(), m_source.c_str());
            }
            auto before = std::chrono::steady_clock::now();
            try {
                p_create_subsolver(std::move(subproblem));
            } catch (InterruptError&) {
                m_destroy.return_subproblem_on_abort(std::move(subproblem));
                continue;
            }
            auto res = m_solver->solve();
            double time_taken =
                seconds_between(before, std::chrono::steady_clock::now());
            if (!res) {
                m_destroy.return_subproblem_on_abort(
                    m_solver->move_out_subproblem());
                continue;
            }
            LNSSubproblem tmp = m_solver->move_out_subproblem();
            if (!*res) {
                auto old_lb = m_portfolio->get_best_lower_bound();
                if (old_lb < tmp.removed_configurations.size()) {
                    m_portfolio->report_lower_bound(
                        tmp.removed_configurations.size(),
                        tmp.uncovered_universe, m_source.c_str());
                }
                m_portfolio->lns_report_failure(
                    tmp.removed_configurations.size(), time_taken);
                m_destroy.improvement_impossible(std::move(tmp),
                                                 m_solver->mes_vertices());
            } else {
                m_portfolio->lns_report_success(
                    tmp.removed_configurations.size(), time_taken);
                const auto& improvement = m_solver->get_solution();
                PartialSolution improved =
                    m_destroy.improve_destroyed(improvement);
                if (m_portfolio->report_solution(improved, m_source.c_str())) {
                    m_destroy.return_subproblem_on_success(std::move(tmp),
                                                           std::move(improved));
                    std::size_t new_size = m_destroy.total_num_configs();
                    {
                        std::unique_lock l{*m_element_mutex};
                        if (new_size < m_destroy_seen_solution) {
                            m_destroy_seen_solution = new_size;
                        }
                    }
                } else {
                    m_destroy.return_subproblem_on_abort(std::move(tmp));
                }
            }
        }
        m_local_recorder->store_event("LNS_WORKER_EXITING",
                                      {{"worker_id", m_destroy.worker_index()}},
                                      "worker_id");
    }

  private:
    enum class DestroyPhaseResult { GLOBAL_OPT, INTERRUPTED, SUBPROBLEM_FOUND };

    DestroyPhaseResult p_destroy_phase() {
        {
            std::unique_lock l{*m_element_mutex};
            if (m_solver) {
                m_solver.reset();
            }
        }
        try {
            p_update_destroy_if_needed();
            if (get_and_clear_interrupt_flag() || m_termination_flag->load()) {
                return DestroyPhaseResult::INTERRUPTED;
            }
            std::size_t num_removed = m_destroy.destroy(
                m_portfolio->lns_select_removed_class_count());
            if (num_removed == 0) {
                m_portfolio->report_lower_bound(
                    m_destroy.total_num_configs(),
                    m_portfolio->get_infeasibility_map().collect_vertices(),
                    "LNS whole solution destruction");
                return DestroyPhaseResult::GLOBAL_OPT;
            }
            return DestroyPhaseResult::SUBPROBLEM_FOUND;
        } catch (InterruptError&) {
            return DestroyPhaseResult::INTERRUPTED;
        }
    }

    void p_update_destroy_if_needed() {
        if (p_mes_update_looks_necessary()) {
            m_destroy.update_global_mes_if_better(m_portfolio->get_best_mes());
        }
        if (!p_update_looks_necessary())
            return;
        PartialSolution new_full = m_portfolio->get_best_solution();
        {
            std::unique_lock l{*m_element_mutex};
            m_destroy_seen_solution = new_full.size();
        }
        m_destroy.update_full_solution_if_better(std::move(new_full));
    }

    bool p_mes_update_looks_necessary() const {
        return m_portfolio->get_best_mes_size() >
               m_destroy.best_global_mes().size();
    }

    bool p_update_looks_necessary() const {
        std::unique_lock l{*m_element_mutex};
        return m_portfolio->get_best_solution_size() <
               m_destroy.total_num_configs();
    }

    void p_interrupt() {
        m_interrupt_flag.store(true);
        if (m_solver) {
            m_solver->abort();
        }
    }

    void p_create_subsolver(LNSSubproblem&& subproblem) {
        SharedDBPropagator prop{m_clauses};
        // relies on the subsolver guarantee that,
        // on interrupt exceptions, subproblem remains usable
        auto subsolver = std::make_unique<SubproblemSolver>(
            m_portfolio, std::move(subproblem), std::move(prop),
            m_local_recorder, m_destroy.worker_index());
        std::unique_lock l{*m_element_mutex};
        m_solver = std::move(subsolver);
    }

    static std::atomic<std::size_t> m_worker_counter;
    PortfolioSolver* m_portfolio;
    std::mutex* m_element_mutex;
    std::atomic<bool>* m_termination_flag;
    std::atomic<bool> m_interrupt_flag;
    EventRecorder* m_local_recorder;
    std::string m_source;
    ClauseDB* m_clauses;
    SharedDBPropagator m_propagator;
    LNSDestroy m_destroy;
    std::size_t m_destroy_seen_solution;
    std::unique_ptr<SubproblemSolver> m_solver;
};

template <typename S>
std::atomic<std::size_t> SubproblemLNSSolverCore<S>::m_worker_counter{0};

/**
 * Portfolio worker that does LNS using
 * the CliqueSatDSaturSolver class template.
 * At any point in time, this is in one of
 * two stages: either we are looking for a
 * subproblem that we can hopefully solve
 * to improve the current solution (m_solver == std::nullopt),
 * or we are currently solving such a subproblem
 * (m_solver != std::nullopt).
 */
template <typename IncrementalSolverType>
class CliqueSatDSaturLNSElement : public PortfolioElement {
  public:
    using IncrementalSolver = IncrementalSolverType;
    using CSDSSolver = CliqueSatDSaturSolver<IncrementalSolver>;

    CliqueSatDSaturLNSElement(PortfolioSolver* solver, std::size_t upper_bound,
                              std::size_t worker_index)
        : PortfolioElement(solver), m_worker_index(worker_index),
          m_clauses(nullptr), m_last_seen_upper_bound(upper_bound),
          m_our_old_solution(upper_bound),
          m_current_subproblem(solver->get_best_solution()),
          m_should_interrupt(false) {}

    void interrupt_if_necessary(const InterruptionCheckInfo& info) override {
        if (should_terminate.load()) {
            p_interrupt();
        }
        if (events & static_cast<EventMask>(PortfolioEvent::BETTER_UPPER_BOUND))
        {
            if (info.best_upper_bound < m_last_seen_upper_bound) {
                m_last_seen_upper_bound = info.best_upper_bound;
                if (m_last_seen_upper_bound < m_our_old_solution) {
                    p_interrupt();
                }
            }
        }
        events = 0;
    }

    void main() override {
        p_init();
        set_interrupt_flag_ptr(&m_should_interrupt);
        m_last_start_time = Clock::now();
        bool begun_selection = false;
        while (!should_terminate.load()) {
            if (!begun_selection)
                m_local_recorder.store_event("BEGIN_SUBPROBLEM_SELECTION");
            else
                m_local_recorder.store_event("RETRY_SUBPROBLEM_SELECTION");
            begun_selection = true;
            p_switch_to_subproblem_selection_stage();
            if (!p_select_next_subproblem())
                continue;
            if (!p_create_solver())
                continue;
            p_switch_to_solving_stage();
            p_solve_cover();
            begun_selection = false;
        }
        m_local_recorder.store_event("LNS_WORKER_EXITING",
                                     {{"description", get_description()}},
                                     "description");
    }

    const EventRecorder* get_recorder() const override {
        return &m_local_recorder;
    }

    void synchronize_recorder(const EventRecorder& other) override {
        m_local_recorder.synchronize_with(other);
    }

    void set_recorder_quiet(bool quiet) override {
        m_local_recorder.set_print_events(!quiet);
    }

    std::string get_description() const override {
        return "C&P/SATDSatur LNS Worker #" + std::to_string(m_worker_index);
    }

  private:
    void p_init() {
        m_clauses = &solver->get_clauses();
        m_empty_propagator.emplace(m_clauses);
        m_fast_clique_builder.emplace(*m_empty_propagator);
    }

    /**
     * Select the next subproblem we try to solve.
     * If we found a subproblem, we also have computed
     * the filtered vertex set and a heuristic clique.
     *
     * @return true if we found a subproblem, false if
     *         we were told to terminate/interrupt or
     *         did not find a suitable subproblem.
     */
    bool p_select_next_subproblem() {
        m_current_subproblem = solver->get_best_solution();
        if (!p_updated_current_subproblem())
            return false;
        if (!p_destroy())
            return false;
        m_local_recorder.store_event("SELECTED_NEXT_SUBPROBLEM");
        return true;
    }

    /**
     * Destroy the current solution to find the next subproblem
     * we try to solve to improve the solution.
     * Must set the filtered vertex set and heuristic clique.
     *
     * @return true if we should continue with this subproblem,
     *         false if we should retry.
     */
    bool p_destroy() {
        std::size_t goal_removed = solver->lns_select_removed_class_count();
        if (goal_removed < 3) {
            goal_removed = 3;
        }
        if (goal_removed > m_current_subproblem.size()) {
            m_local_recorder.store_event("DESTROY_WHOLE_SOLUTION");
            goal_removed = m_current_subproblem.size();
            return p_destroy_random(goal_removed);
        }
        double destroy_method_sel =
            std::uniform_real_distribution<double>{0.0, 1.0}(sammy::rng());
        if (destroy_method_sel < RANDOM_DESTRUCTION_PROB) {
            m_local_recorder.store_event("DESTROY_RANDOM",
                                         {{"goal_removed", goal_removed}},
                                         "goal_removed");
            return p_destroy_random(goal_removed);
        } else if (destroy_method_sel <
                   RANDOM_DESTRUCTION_PROB + RANDOM_WITH_TABLE_PROB)
        {
            p_update_cached_mes_and_table();
            m_local_recorder.store_event("DESTROY_RANDOM_USING_BEST_MES",
                                         {{"goal_removed", goal_removed}},
                                         "goal_removed");
            return p_destroy_random_with_table(goal_removed);
        } else {
            m_local_recorder.store_event("DESTROY_RANDOM_AND_LEAST_UNIQUE",
                                         {{"goal_removed", goal_removed}},
                                         "goal_removed");
            return p_destroy_random_and_least_unique(goal_removed);
        }
    }

    /**
     * Update the cached MES and the table
     * of classes covering MES vertices.
     */
    void p_update_cached_mes_and_table() {
        m_cached_best_mes = solver->get_best_mes();
        m_class_index_to_mes_index.clear();
        m_class_index_to_mes_index.resize(m_current_subproblem.size());
        m_mes_vertices_nonuniquely_covered.assign(m_cached_best_mes.size(),
                                                  false);
        for (auto& v : m_mes_index_to_class_indices)
            v.clear();
        m_mes_index_to_class_indices.resize(m_cached_best_mes.size());
        std::size_t idx = 0;
        for (Vertex v : m_cached_best_mes) {
            m_current_subproblem.find_covering_classes(
                v.first, v.second, [&](Index cls) {
                    m_class_index_to_mes_index[cls] = idx;
                    m_mes_index_to_class_indices[idx].push_back(cls);
                });
            if (m_mes_index_to_class_indices[idx].size() != 1) {
                m_mes_vertices_nonuniquely_covered[idx].set();
            }
            ++idx;
        }
    }

    /**
     * Check whether there is some hope for finding an improvement
     * by removing the classes with the given indices, or if
     * each of the indices covers a unique vertex in the MES.
     */
    bool p_index_set_could_improve(const std::vector<Index>& indices) {
        DynamicBitset mes_vertices_covered = m_mes_vertices_nonuniquely_covered;
        m_initial_clique.clear();
        for (Index i : indices) {
            auto v = m_class_index_to_mes_index[i];
            if (v) {
                if (!mes_vertices_covered[*v]) {
                    mes_vertices_covered[*v].set();
                    m_initial_clique.push_back(m_cached_best_mes[*v]);
                }
            }
        }
        if (m_initial_clique.size() == indices.size())
            return false;
        return true;
    }

    /**
     * Remove goal_removed randomly selected classes.
     * Repeat the removal process up to max_iterations times
     * if the removal process as long as the removal process
     * results in a 'hopeless' subproblem according to the best MES.
     */
    bool p_remove_random_if_improvement_possible(std::size_t goal_removed,
                                                 std::size_t max_iterations) {
        p_update_cached_mes_and_table();
        if (p_interrupt_check())
            return false;
        auto& rng = sammy::rng();
        for (std::size_t iter = 0; iter < max_iterations; ++iter) {
            std::vector<Index> indices =
                vector(range(Index(m_current_subproblem.size())));
            std::shuffle(indices.begin(), indices.end(), rng);
            indices.resize(goal_removed);
            std::sort(indices.begin(), indices.end());
            if (p_index_set_could_improve(indices)) {
                m_currently_removed =
                    m_current_subproblem.remove_assignments(indices);
                return true;
            }
        }
        return false;
    }

    /**
     * Destroy the current solution, making use of the
     * table of classes covering MES vertices to find
     * a deletion that may cause an improvement.
     */
    bool p_destroy_random_with_table(std::size_t goal_removed) {
        std::vector<Index> classes_without_mes_vertex;
        std::vector<std::size_t> multiply_covered_mes_vertices;
        for (Index i : range(Index(m_current_subproblem.size()))) {
            if (!m_class_index_to_mes_index[i]) {
                classes_without_mes_vertex.push_back(i);
            }
        }
        for (std::size_t vi : m_mes_vertices_nonuniquely_covered.ones()) {
            multiply_covered_mes_vertices.push_back(vi);
        }
        if (classes_without_mes_vertex.empty() &&
            multiply_covered_mes_vertices.empty())
        {
            // we actually have optimality here; may reach this depending on
            // race condition
            should_terminate.store(true);
            return false;
        }
        auto& rng = sammy::rng();
        std::uniform_real_distribution<double> sel(0.0, 1.0);
        std::vector<Index> selected_classes;
        if (!classes_without_mes_vertex.empty() &&
            (multiply_covered_mes_vertices.empty() || sel(rng) < 0.5))
        {
            std::size_t select = (std::min)(goal_removed / 2 + 1,
                                            classes_without_mes_vertex.size());
            std::shuffle(classes_without_mes_vertex.begin(),
                         classes_without_mes_vertex.end(), rng);
            selected_classes.assign(classes_without_mes_vertex.begin(),
                                    classes_without_mes_vertex.begin() +
                                        select);
        } else {
            std::size_t srnd = std::uniform_int_distribution<std::size_t>(
                0, multiply_covered_mes_vertices.size() - 1)(rng);
            std::size_t mvi = multiply_covered_mes_vertices[srnd];
            std::vector<Index> classes = m_mes_index_to_class_indices[mvi];
            selected_classes = m_mes_index_to_class_indices[mvi];
            if (selected_classes.size() > goal_removed) {
                std::shuffle(selected_classes.begin(), selected_classes.end(),
                             rng);
                selected_classes.resize(goal_removed);
            }
        }
        m_initial_clique.clear();
        std::sort(selected_classes.begin(), selected_classes.end());
        auto is_in = [&](Index idx) {
            return std::binary_search(selected_classes.begin(),
                                      selected_classes.end(), idx);
        };
        std::vector<Index> additional =
            vector(range(Index(m_current_subproblem.size())));
        additional.erase(
            std::remove_if(additional.begin(), additional.end(), is_in),
            additional.end());
        std::shuffle(additional.begin(), additional.end(), rng);
        additional.resize(goal_removed - selected_classes.size());
        auto inserted = selected_classes.insert(
            selected_classes.end(), additional.begin(), additional.end());
        std::sort(inserted, selected_classes.end());
        std::inplace_merge(selected_classes.begin(), inserted,
                           selected_classes.end());
        m_currently_removed =
            m_current_subproblem.remove_assignments(selected_classes);
        return p_post_destroy();
    }

    /**
     * Destroy m_current_solution by removing
     * goal_removed randomly selected elements.
     */
    bool p_destroy_random(std::size_t goal_removed) {
        bool res = p_remove_random_if_improvement_possible(goal_removed, 10);
        if (!res) {
            m_local_recorder.store_event("DESTROY_RANDOM_FAILED",
                                         {{"goal_removed", goal_removed}},
                                         "goal_removed");
            if (p_interrupt_check())
                return false;
            return p_destroy_random_with_table(goal_removed);
        }
        return p_post_destroy();
    }

    /**
     * Destroy m_current_solution by removing
     * roughly 1/3 * goal_removed randomly selected configurations
     * followed by roughly 2/3 * goal_removed configurations
     * that cover the least amount of uniquely covered elements.
     */
    bool p_destroy_random_and_least_unique(std::size_t goal_removed) {
        std::size_t goal_removed_random =
            (std::min)(goal_removed, goal_removed / 3 + 1);
        if (!p_remove_random_if_improvement_possible(goal_removed_random, 10)) {
            if (p_interrupt_check())
                return false;
            return p_destroy_random_with_table(goal_removed);
        }
        std::size_t goal_removed_least_unique =
            goal_removed - goal_removed_random;
        std::vector<std::size_t> unique_per_class(m_current_subproblem.size(),
                                                  0);
        m_current_subproblem.iterate_all_uniquely_covered_with_class(
            [&](Index i, Lit, Lit) { unique_per_class[i] += 1; });
        std::vector<Index> indices =
            vector(range(Index(m_current_subproblem.size())));
        std::nth_element(indices.begin(),
                         indices.begin() + goal_removed_least_unique,
                         indices.end(), [&](Index a, Index b) {
                             return unique_per_class[a] < unique_per_class[b];
                         });
        indices.resize(goal_removed_least_unique);
        std::sort(indices.begin(), indices.end());
        auto further_removed = m_current_subproblem.remove_assignments(indices);
        for (auto& r : further_removed)
            m_currently_removed.emplace_back(std::move(r));
        return p_post_destroy();
    }

    /**
     * After destruction of m_current_solution,
     * compute the filtered set of uncovered vertices
     * and an initial clique/MES.
     */
    bool p_post_destroy() {
        if (p_interrupt_check())
            return false;
        m_local_recorder.store_event("BEGIN_FILTER_UNCOVERED");
        m_filtered_uncovered.clear();
        m_current_subproblem.iterate_all_uncovered([&](Lit lmin, Lit lmax) {
            m_filtered_uncovered.emplace_back(lmin, lmax);
        });
        solver->implied_cache().remove_implied(m_filtered_uncovered,
                                               *m_empty_propagator);
        m_local_recorder.store_event("DONE_FILTER_UNCOVERED");
        if (m_filtered_uncovered.empty()) {
            solver->report_solution(m_current_subproblem,
                                    "empty uncovered set");
            return false;
        }
        m_local_recorder.store_event("BEGIN_INITIAL_MES");
        m_initial_clique = p_find_initial_clique();
        m_local_recorder.store_event(
            "DONE_INITIAL_MES", {{"size", m_initial_clique.size()}}, "size");
        if (m_initial_clique.size() >= m_currently_removed.size()) {
            m_local_recorder.store_event("INITIAL_MES_PRECLUDES_IMPROVEMENT");
            solver->report_mes(m_initial_clique, "LNS initial MES heuristic");
            return false;
        }
        return !p_interrupt_check();
    }

    /**
     * Notify that we updated the current subproblem from
     * the best solution from the PortfolioSolver.
     * @return true if we should continue with this subproblem,
     *         false if we should retry or terminate.
     */
    bool p_updated_current_subproblem() {
        std::size_t new_size = m_current_subproblem.size();
        std::unique_lock l{mutex};
        m_our_old_solution = new_size;
        if (m_last_seen_upper_bound < m_our_old_solution) {
            return false;
        }
        if (m_our_old_solution < m_last_seen_upper_bound) {
            m_last_seen_upper_bound = m_our_old_solution;
        }
        return !should_terminate.load();
    }

    /**
     * Find an initial clique on m_filtered_uncovered
     * by heuristic means.
     */
    std::vector<Vertex> p_find_initial_clique() {
        assert(m_empty_propagator->get_current_level() == 0);
        std::vector<Vertex> best_clique =
            m_fast_clique_builder->random_multistart_best_clique_known_valid(
                10, m_filtered_uncovered);
        if (m_initial_clique.size() > best_clique.size()) {
            std::vector<Vertex> missing = m_initial_clique;
            auto rem =
                std::remove_if(missing.begin(), missing.end(), [&](Vertex v) {
                    return std::binary_search(m_filtered_uncovered.begin(),
                                              m_filtered_uncovered.end(), v);
                });
            missing.erase(rem, missing.end());
            if (!missing.empty()) {
                std::sort(missing.begin(), missing.end());
                auto ins = m_filtered_uncovered.insert(
                    m_filtered_uncovered.end(), missing.begin(), missing.end());
                std::inplace_merge(m_filtered_uncovered.begin(), ins,
                                   m_filtered_uncovered.end());
            }
            best_clique = m_initial_clique;
            best_clique = m_fast_clique_builder->compute_clique_known_valid(
                best_clique.begin(), best_clique.end(), m_filtered_uncovered);
        }
        if (p_walk_cached_cliques(best_clique)) {
            std::size_t old_size = best_clique.size();
            best_clique = m_fast_clique_builder->compute_clique_known_valid(
                best_clique.begin(), best_clique.end(), m_filtered_uncovered);
            if (best_clique.size() > old_size) {
                solver->add_clique(best_clique);
            }
        }
        return best_clique;
    }

    bool p_walk_cached_cliques(std::vector<Vertex>& best_clique) {
        auto cview = solver->clique_cache_view();
        std::vector<Vertex> current;
        auto used = cview.end();
        for (auto i = cview.begin(), e = cview.end(); i != e; ++i) {
            auto clique = *i;
            if (clique.size() <= best_clique.size())
                continue;
            current.clear();
            std::copy_if(clique.begin(), clique.end(),
                         std::back_inserter(current), [&](Vertex v) {
                             return std::binary_search(
                                 m_filtered_uncovered.begin(),
                                 m_filtered_uncovered.end(), v);
                         });
            if (current.size() > best_clique.size()) {
                best_clique = current;
                used = i;
            }
        }
        if (used != cview.end()) {
            solver->clique_was_used(cview, used);
            return true;
        } else {
            solver->add_clique(best_clique);
            return false;
        }
    }

    /**
     * Create the solver.
     * @return true if we should continue, false if we were interrupted.
     */
    bool p_create_solver() {
        m_local_recorder.store_event("BEGIN_CREATE_SUBPROBLEM_SOLVER");
        try {
            std::size_t best_mes_size = solver->get_best_mes_size();
            std::size_t best_lower_bound = solver->get_best_lower_bound();
            solver->report_subproblem(m_filtered_uncovered, m_initial_clique,
                                      m_current_subproblem, m_currently_removed,
                                      best_mes_size, best_lower_bound,
                                      "satdsatur");
            m_solver.emplace(m_filtered_uncovered,
                             &solver->get_infeasibility_map(), *m_clauses,
                             m_initial_clique, best_mes_size, best_lower_bound,
                             m_currently_removed);
            m_solver->set_event_recorder(&m_local_recorder);
        } catch (const InterruptError&) {
            m_local_recorder.store_event("ABORTED_CREATE_SUBPROBLEM_SOLVER");
            return false;
        }
        m_local_recorder.store_event("DONE_CREATE_SUBPROBLEM_SOLVER");
        if (get_and_clear_interrupt_flag())
            return false;
        return !p_interrupt_check();
    }

    /**
     * Start the solver and wait for it to finish or be aborted.
     * Automatically reports improved solutions to the
     * PortfolioSolver instance controlling this worker.
     */
    void p_solve_cover() {
        auto state = m_solver->solve();
        if (state == CSDSSolver::SolveResult::ABORTED) {
            m_local_recorder.store_event("ABORTED_SUBPROBLEM_SOLVE");
            return;
        }
        m_local_recorder.store_event("DONE_SUBPROBLEM_SOLVE");
        if (m_solver->improved_global_bound()) {
            solver->report_lower_bound(m_solver->get_best_bound(),
                                       m_solver->get_best_bound_subgraph(),
                                       "C&P/SATDSatur LNS");
        }
        if (m_solver->improved_mes()) {
            solver->report_mes(m_solver->get_best_mes(), "C&P/SATDSatur LNS");
        }
        double time_taken = seconds_between(m_last_start_time, Clock::now());
        if (state == CSDSSolver::SolveResult::SOLUTION_WAS_OPTIMAL) {
            solver->lns_report_failure(m_currently_removed.size(), time_taken);
        } else {
            PartialSolution subgraph_solution =
                m_solver->get_partial_solution();
            for (const auto& assignment : subgraph_solution.assignments()) {
                m_current_subproblem.add_assignment(assignment);
            }
            solver->report_solution(m_current_subproblem, "C&P/SATDSatur LNS");
            solver->lns_report_success(m_currently_removed.size(), time_taken);
        }
    }

    /**
     * Check whether we should interrupt what we are doing
     * (called mostly during subproblem selection stage).
     * This can be due to a better solution, optimality,
     * or some termination criterion.
     *
     * @return true if we should interrupt, false otherwise.
     */
    bool p_interrupt_check() {
        std::unique_lock l{mutex};
        if (should_terminate.load() ||
            m_our_old_solution > m_last_seen_upper_bound)
        {
            return true;
        }
        return false;
    }

    /**
     * Cause the solver to terminate
     * as soon as possible.
     * Called with lock held.
     */
    void p_interrupt() {
        if (m_have_solver) {
            m_solver->abort();
        } else {
            m_should_interrupt.store(true);
        }
    }

    /**
     * Switch to the subproblem selection stage.
     * This means taking the lock, resetting m_have_solver,
     * and destroying m_solver.
     */
    void p_switch_to_subproblem_selection_stage() {
        std::unique_lock l{mutex};
        m_have_solver = false;
        m_solver.reset();
    }

    /**
     * Switch to the subproblem solving stage.
     * This means taking the lock and setting m_have_solver.
     * m_solver MUST BE non-empty before calling this.
     */
    void p_switch_to_solving_stage() {
        std::unique_lock l{mutex};
        m_local_recorder.store_event(
            "BEGIN_SUBPROBLEM_SOLVE",
            {{"removed_configurations", m_currently_removed.size()},
             {"old_solution_value", m_our_old_solution},
             {"reduced_uncovered", m_solver->get_graph_size()}},
            "removed_configurations", "reduced_uncovered");
        assert(!!m_solver);
        m_have_solver = true;
        m_last_start_time = Clock::now();
    }

    /**
     * The stage indicator variable.
     * Invariant: m_solver is != std::nullopt if m_have_solver is true.
     * The worker thread will always populate
     * m_solver before setting it to true.
     * This variable is only changed with the mutex held.
     */
    bool m_have_solver{false};

    /**
     * The solver, if we currently have it.
     */
    std::optional<CSDSSolver> m_solver;

    /**
     * Index of this worker (mostly for debugging/logging/time measurement
     * purposes).
     */
    std::size_t m_worker_index;

    /**
     * Pointer to our clauses.
     */
    ClauseDB* m_clauses;

    /**
     * Empty propagator we can reuse.
     */
    std::optional<SharedDBPropagator> m_empty_propagator;

    /**
     * Clique heuristic.
     */
    std::optional<FastCliqueBuilder> m_fast_clique_builder;

    /**
     * Our local recorder.
     */
    EventRecorder m_local_recorder;

    /**
     * The last seen upper bound.
     */
    std::size_t m_last_seen_upper_bound;

    /**
     * The value of the solution we destroyed to obtain
     * our current repair subproblem.
     */
    std::size_t m_our_old_solution;

    /**
     * The last time we started looking for an improvement.
     */
    Clock::time_point m_last_start_time;

    /**
     * The partial solution we currently try to extend.
     */
    PartialSolution m_current_subproblem;

    /**
     * The set of configurations we removed to
     * obtain the current subproblem.
     */
    std::vector<DynamicBitset> m_currently_removed;

    /**
     * A filtered set of uncovered vertices
     * that guarantees all uncovered vertices are covered.
     */
    std::vector<Vertex> m_filtered_uncovered;

    /**
     * Initial clique found on m_filtered_uncovered.
     */
    std::vector<Vertex> m_initial_clique;

    /**
     * The best global MES at the time of destruction of
     * the solution we currently work on.
     */
    std::vector<Vertex> m_cached_best_mes;

    /**
     * Index of the vertex in m_cached_best_mes
     * that is covered by class with index i, if any.
     */
    std::vector<std::optional<std::size_t>> m_class_index_to_mes_index;

    /**
     * Classes covering vertex i in m_cached_best_mes.
     */
    std::vector<std::vector<Index>> m_mes_index_to_class_indices;

    /**
     * A bitset of MES vertices that are not uniquely covered.
     */
    DynamicBitset m_mes_vertices_nonuniquely_covered;

    /**
     * Interruption flag for the worker thread.
     */
    std::atomic<bool> m_should_interrupt;
};

} // namespace sammy

#endif
==> ./basic_cuda_iteration.h <==
#ifndef SAMMY_BASIC_CUDA_ITERATION_H_INCLUDED_
#define SAMMY_BASIC_CUDA_ITERATION_H_INCLUDED_

#include <cstddef>
#include <cstdint>
#include <exception>
#include <stdexcept>

/**
 * This file exists as a workaround to bugs with nvcc (or gcc/clang/stdlibc++,
 * depends on who you ask) which causes compilation problems when nvcc
 * is exposed to some standard C++ headers.
 */

/**
 * Conditional compilation (if no CUDA support,
 * the methods in this file cannot be used and
 * will lead to linker errors, but the file can be included).
 */
#ifndef SAMMY_CUDA_SUPPORTED
using cudaError_t = int;
#define SAMMY_HD
#else
#include <cuda_runtime.h>
#define SAMMY_HD __host__ __device__
#endif
#include <atomic>

namespace sammy {

/**
 * Exception thrown when CUDA errors occur.
 */
class CUDAError : public std::runtime_error {
  public:
    explicit CUDAError(cudaError_t err);
    explicit CUDAError(cudaError_t err, const std::string& info);

    /**
     * Get the CUDA error code.
     */
    cudaError_t error_code() const noexcept;

    /**
     * Get the CUDA error string.
     */
    const char* error_string() const noexcept;

    /**
     * Get the CUDA error name.
     */
    const char* error_name() const noexcept;

  private:
    cudaError_t m_error_code;
};

/**
 * CUDA use mode: if available, use CUDA; on error, fallback to CPU.
 */
static constexpr int CUDA_USAGE_IF_AVAILABLE = 0;

/**
 * CUDA use mode: use CUDA; on error, fail.
 */
static constexpr int CUDA_USAGE_FORCED = 1;

/**
 * CUDA use mode: do not use CUDA.
 */
static constexpr int CUDA_USAGE_DISABLED = 2;

/**
 * CUDA use mode: had an error previously, don't use CUDA.
 */
static constexpr int CUDA_USAGE_HAD_ERROR = 3;

#ifdef SAMMY_CUDA_SUPPORTED

/**
 * Usage info flag for CUDA.
 */
inline std::atomic<int>& cuda_usage_info() {
    static std::atomic<int> usage_info{CUDA_USAGE_IF_AVAILABLE};
    return usage_info;
}

inline bool should_use_cuda() {
    switch (cuda_usage_info()) {
    default:
    case CUDA_USAGE_IF_AVAILABLE:
    case CUDA_USAGE_FORCED:
        return true;

    case CUDA_USAGE_HAD_ERROR:
    case CUDA_USAGE_DISABLED:
        return false;
    }
}

inline void had_cuda_error_(const CUDAError& error) {
    auto& info = cuda_usage_info();
    int loaded = info.load();
    if (loaded == CUDA_USAGE_FORCED) {
        throw error;
    }
    info.store(CUDA_USAGE_HAD_ERROR);
}

inline void had_cuda_error(const CUDAError& error) noexcept {
    had_cuda_error_(error);
}

inline void set_cuda_mode(int mode) { cuda_usage_info().store(mode); }

#else

inline constexpr bool should_use_cuda() { return false; }

inline void had_cuda_error(const CUDAError& error) noexcept {}

inline void set_cuda_mode(int mode) {
    if (mode == CUDA_USAGE_FORCED) {
        throw std::runtime_error(
            "Requested CUDA mode 'force' without CUDA support compiled in!");
    }
}

#endif

namespace detail {

void cuda_malloc_or_throw(void** ptr, std::size_t bytes);
void cuda_memcpy_htd_or_throw(void* dst, const void* src, std::size_t bytes);
void cuda_memcpy_dth_or_throw(void* dst, const void* src, std::size_t bytes);
void cuda_free(void* ptr);

/**
 * We produce rows of output.
 * Each row is a bitset of size num_literals.
 * Each bitset consists of 32-bit integers.
 * THREADS_PER_ROW threads produce one output bitset together.
 */
constexpr static SAMMY_HD std::size_t BLOCK_SIZE() { return 256; }
constexpr static SAMMY_HD std::size_t THREADS_PER_ROW() { return 16; }
constexpr static SAMMY_HD std::size_t ROWS_PER_BLOCK() {
    return BLOCK_SIZE() / THREADS_PER_ROW();
}
constexpr static SAMMY_HD std::size_t GOAL_ROWS_PER_CALL() {
    return 256 * ROWS_PER_BLOCK();
}

void cuda_call_bit_filter_kernel(
    const std::uint32_t* bit_data, std::size_t u32_per_bitset,
    const std::uint32_t* classes_with_literal, std::size_t num_literals,
    const std::uint32_t* classes_with_literal_offsets,
    std::uint32_t* output_buffer, std::size_t begin_row, std::size_t num_rows);

void call_cuda_extract_kernel(const std::uint32_t* device_bit_data,
                              std::size_t u32_per_bitset,
                              std::uint32_t* output_buffer,
                              std::size_t num_literals, std::size_t begin_var,
                              std::size_t num_vars, std::size_t num_classes);

} // namespace detail

} // namespace sammy

#endif
==> ./satdsatur_colorer.h <==
#ifndef SAMMY_SATDSATUR_COLORER_H_INCLUDED_
#define SAMMY_SATDSATUR_COLORER_H_INCLUDED_

#include "algorithm_ex.h"
#include "best_k.h"
#include "class_completer.h"
#include "partial_solution.h"
#include "shared_db_propagator.h"
#include "universe_subgraph.h"
#include "vertex_operations.h"
#include <atomic>
#include <variant>

namespace sammy {

template <typename IncrementalSolver> class SATDSaturCoverer {
    // Options for keeping open_classes updated?
    // - use only the graph datastructure
    //   * iterate neighbors on addition
    //   * no deletion/conflict learning
    //   * un 'unexpected' conflict, remove class
    //   * resets of classes on SAT coverage
    // - use graph datastructure and propagated literals
    //   * iterate neighbors on addition
    //   * have/iterate list of vertices with literal
    //   * no deletion/conflict learning
    //   * on unexpected conflict, remove class
    //   * the only way to remove vertices from a class
    //     is a successful SAT return (recoloring)
    //   * how to track covering classes?
    struct VertexInfo {
        DynamicBitset
            open_classes; //< a superset of the open classes for this vertex
        std::size_t num_open_classes; //< number of set bits in open_classes
        bool in_some_class; //< number of classes containing this vertex
        std::size_t degree; //< degree of the vertex in the graph

        /**
         * Mark the given class as unavailable for this vertex.
         */
        bool remove_from(std::size_t class_index) noexcept {
            if (open_classes[class_index]) {
                open_classes[class_index].reset();
                --num_open_classes;
                return true;
            }
            return false;
        }
    };

    struct CompareByInfo {
        bool operator()(std::size_t v1, std::size_t v2) const noexcept {
            const VertexInfo& i1 = that->m_vertex_info[v1];
            const VertexInfo& i2 = that->m_vertex_info[v2];
            if (i2.in_some_class)
                return !i1.in_some_class;
            if (i1.in_some_class)
                return false;
            if (i1.num_open_classes < i2.num_open_classes)
                return true;
            if (i2.num_open_classes < i1.num_open_classes)
                return false;
            return i2.degree < i1.degree;
        }

        const SATDSaturCoverer* that;
    };

    struct ColorClass {
        ColorClass(SATDSaturCoverer* that, const SharedDBPropagator& propagator,
                   std::size_t index)
            : propagator(propagator), index(index) {
            p_init_info_from_trail(that);
        }

        ColorClass(SATDSaturCoverer* that, const SharedDBPropagator& propagator,
                   std::size_t index, std::size_t initial_vertex)
            : propagator(propagator), index(index) {
            push_vertex(this->propagator,
                        that->m_subgraph->vertex(initial_vertex));
            p_init_info_from_trail(that);
        }

        void add_vertex(SATDSaturCoverer* that, std::size_t vertex) {
            UniverseSubgraph* g = that->m_subgraph;
            std::size_t tlength = propagator.get_trail().size();
            if (push_vertex(propagator, g->vertex(vertex)) < 0) {
                throw std::logic_error(
                    "Called add_vertex on incompatible class!");
            }
            for (std::size_t vo : g->matrix_row(vertex).ones()) {
                that->m_vertex_info[vo].remove_from(index);
            }
            for (Lit l : IteratorRange{propagator.get_trail().begin() + tlength,
                                       propagator.get_trail().end()})
            {
                Lit lneg = lit::negate(l);
                for (std::size_t v_pot : that->m_vertices_with_literal[l]) {
                    Vertex v = g->vertex(v_pot);
                    if (propagator.is_true(v.first) &&
                        propagator.is_true(v.second))
                    {
                        that->m_vertex_info[v_pot].in_some_class = true;
                    }
                }
                for (std::size_t v_out : that->m_vertices_with_literal[lneg]) {
                    that->m_vertex_info[v_out].remove_from(index);
                }
            }
        }

        bool can_add(SATDSaturCoverer* that, std::size_t vertex) {
            Vertex v = that->m_subgraph->vertex(vertex);
            return can_push(propagator, v);
        }

        void reset(SATDSaturCoverer* that,
                   const DynamicBitset& expl_vertex_is_in) {
            propagator.reset_or_throw();
            p_init_info_from_trail(that);
            for (std::size_t vi_i : expl_vertex_is_in.ones()) {
                std::size_t vi = that->m_explicit_cover_order[vi_i];
                add_vertex(that, vi);
            }
        }

        SharedDBPropagator propagator;
        std::size_t index;

      private:
        void p_init_info_from_trail(SATDSaturCoverer* that) {
            UniverseSubgraph* g = that->m_subgraph;
            for (Lit l : propagator.get_trail()) {
                for (std::size_t v_pot : that->m_vertices_with_literal[l]) {
                    Vertex v = g->vertex(v_pot);
                    if (v.first == l) {
                        if (propagator.is_true(v.second)) {
                            that->m_vertex_info[v_pot].in_some_class = true;
                        }
                    }
                }
                Lit lneg = lit::negate(l);
                for (std::size_t v_out : that->m_vertices_with_literal[lneg]) {
                    that->m_vertex_info[v_out].remove_from(index);
                }
            }
        }
    };

    using SatLit = typename IncrementalSolver::Lit;
    using LitOrFixed = std::variant<SatLit, bool>;

    UniverseSubgraph* m_subgraph;
    std::vector<std::vector<std::size_t>> m_vertices_with_literal;
    std::vector<VertexInfo> m_vertex_info;
    std::vector<ColorClass> m_color_classes;
    BestK<std::size_t, CompareByInfo> m_best_infos;
    std::size_t m_clique_size, m_lower_bound, m_upper_bound;
    std::vector<std::size_t> m_explicit_cover_order;
    std::vector<std::size_t> m_true_zero_vertices;
    IncrementalSolver m_incremental_solver;
    // model variables: variable copies for each class;
    // entry [c][v] is the v-th variable of class c
    std::vector<std::vector<LitOrFixed>> m_class_literals;
    // model variables: vertex-in-class variables
    // entry [c][vi] encodes if vertex m_explicit_cover_order[vi] is in class c
    std::vector<std::vector<LitOrFixed>> m_vertex_in_class;
    // model variables: variables indicating that
    // more than some number of colors are needed;
    // entry [c] says that more than c + m_clique_size
    // colors are necessary
    std::vector<LitOrFixed> m_not_enough_colors;
    // buffer for building clauses
    std::vector<SatLit> m_clause_buffer;
    // flag to abort the search
    std::atomic<bool> m_aborted{false};
    // indicates whether we found an optimal coloring
    // or proved no improvement is possible, rather
    // than having been aborted
    bool m_success{false};

    void p_fill_vertices_with_literal() {
        m_vertices_with_literal.resize(
            2 * m_subgraph->get_propagator().db().num_vars());
        std::size_t index = 0;
        for (Vertex v : m_subgraph->vertex_set()) {
            m_vertices_with_literal[v.first].push_back(index);
            m_vertices_with_literal[v.second].push_back(index);
            ++index;
        }
    }

    void p_init_vertex_info() {
        const auto n = m_subgraph->n();
        VertexInfo info{DynamicBitset(m_lower_bound, true), m_lower_bound,
                        false, 0};
        m_vertex_info.resize(n, info);
        for (std::size_t i = 0, n = m_subgraph->n(); i < n; ++i) {
            m_vertex_info[i].degree = m_subgraph->get_degree(i);
        }
    }

    void p_init_classes(std::vector<std::size_t> clique) {
        m_color_classes.reserve(clique.size());
        std::size_t index = 0;
        UniverseSubgraph* g = m_subgraph;
        SharedDBPropagator& propagator = g->get_propagator();
        propagator.reset_or_throw();
        for (std::size_t vi : clique) {
            m_color_classes.emplace_back(this, propagator, index++, vi);
        }
        m_explicit_cover_order = std::move(clique);
    }

    std::size_t
    p_find_true_next_vertex(const std::vector<std::size_t>& candidates) {
        m_true_zero_vertices.clear();
        for (std::size_t c : candidates) {
            VertexInfo& info = m_vertex_info[c];
            if (info.in_some_class)
                continue;
            std::size_t count = 0;
            for (std::size_t cindex : info.open_classes.ones()) {
                if (m_color_classes[cindex].can_add(this, c)) {
                    ++count;
                } else {
                    info.open_classes[cindex].reset();
                    --info.num_open_classes;
                }
            }
            if (count == 0) {
                m_true_zero_vertices.push_back(c);
            }
        }
        if (!m_true_zero_vertices.empty()) {
            return m_true_zero_vertices.front();
        }
        return *std::min_element(candidates.begin(), candidates.end(),
                                 CompareByInfo{this});
    }

    std::optional<std::size_t> p_find_next_vertex() {
        m_best_infos.clear();
        for (std::size_t i : range(m_subgraph->n())) {
            m_best_infos.push(i);
        }
        auto is_covered = [&](std::size_t v) {
            return m_vertex_info[v].in_some_class;
        };
        const auto& best = m_best_infos.elements();
        if (std::all_of(best.begin(), best.end(), is_covered)) {
            return std::nullopt;
        }
        return p_find_true_next_vertex(best);
    }

    static bool p_is_level_vertex(const SharedDBPropagator& propagator,
                                  std::int32_t level, Vertex vertex) {
        Lit l = *propagator.level_begin(level);
        return vertex.first == l || vertex.second == l;
    }

    auto p_implied_by(const SharedDBPropagator& propagator, Vertex vertex) {
        if (propagator.get_current_level() == 0 ||
            !p_is_level_vertex(propagator, 1, vertex))
        {
            return IteratorRange{propagator.get_trail().begin(),
                                 propagator.level_end(0)};
        }
        if (propagator.get_current_level() < 2 ||
            !p_is_level_vertex(propagator, 2, vertex))
        {
            return IteratorRange{propagator.get_trail().begin(),
                                 propagator.level_end(1)};
        }
        return IteratorRange{propagator.get_trail().begin(),
                             propagator.level_end(2)};
    }

    LitOrFixed p_value_in_class(Lit l, std::size_t class_index) {
        std::vector<LitOrFixed>& class_vars = m_class_literals[class_index];
        const bool negate = lit::negative(l);
        LitOrFixed& var = class_vars[lit::var(l)];
        return std::visit(
            overloaded{[&](bool& fixed) {
                           return LitOrFixed{std::in_place_index<1>,
                                             fixed ^ negate};
                       },
                       [&](SatLit& lit) {
                           return LitOrFixed{std::in_place_index<0>,
                                             negate ? -lit : lit};
                       }},
            var);
    }

    template <typename LitIterator>
    void p_class_literals_embed_clause(LitIterator clause_begin,
                                       LitIterator clause_end,
                                       std::size_t class_index) {
        m_clause_buffer.clear();
        bool found_true = false;
        for (Lit l : IteratorRange{clause_begin, clause_end}) {
            LitOrFixed value = p_value_in_class(l, class_index);
            std::visit(overloaded{[&](bool& fixed) {
                                      if (fixed)
                                          found_true = true;
                                  },
                                  [&](SatLit& lit) {
                                      m_clause_buffer.push_back(lit);
                                  }},
                       value);
            if (found_true)
                return;
        }
        if (m_clause_buffer.empty()) {
            throw std::logic_error("UNSAT class (has empty clause)!");
        }
        m_incremental_solver.add_clause(m_clause_buffer.begin(),
                                        m_clause_buffer.end());
    }

    void p_class_literals_embed_formula(const ClauseDB& clauses,
                                        std::size_t class_index) {
        for (auto [l1, l2] : clauses.binary_clauses()) {
            const Lit l[2] = {l1, l2};
            p_class_literals_embed_clause(+l, l + 2, class_index);
        }
        for (CRef c = 1, ndb = clauses.literal_db_size(); c < ndb;
             c = clauses.next_clause(c))
        {
            auto lits = clauses.lits_of(c);
            p_class_literals_embed_clause(lits.begin(), lits.end(),
                                          class_index);
        }
    }

    template <typename FixedRange>
    void p_init_class_literals_with_fixed_range(std::size_t class_index,
                                                FixedRange fixed_range) {
        const SharedDBPropagator& prop =
            m_color_classes[class_index].propagator;
        const Var num_vars = prop.db().num_vars();
        DynamicBitset fixed{num_vars, false};
        for (Lit l : fixed_range) {
            fixed[lit::var(l)].set();
        }
        std::vector<LitOrFixed> class_literals;
        for (Var v : range(num_vars)) {
            if (fixed[v]) {
                bool value = prop.is_true(lit::positive_lit(v));
                class_literals.emplace_back(std::in_place_index<1>, value);
            } else {
                class_literals.emplace_back(std::in_place_index<0>,
                                            m_incremental_solver.new_var());
            }
        }
        m_class_literals.emplace_back(std::move(class_literals));
        p_class_literals_embed_formula(prop.db(), class_index);
    }

    void p_init_class_literals_clique(std::size_t class_index) {
        const SharedDBPropagator& prop =
            m_color_classes[class_index].propagator;
        std::size_t clique_vertex_index = m_explicit_cover_order[class_index];
        Vertex clique_vertex = m_subgraph->vertex(clique_vertex_index);
        p_init_class_literals_with_fixed_range(
            class_index, p_implied_by(prop, clique_vertex));
    }

    void p_init_class_literals_not_clique(std::size_t class_index) {
        const SharedDBPropagator& prop =
            m_color_classes[class_index].propagator;
        p_init_class_literals_with_fixed_range(
            class_index, IteratorRange{prop.level_begin(0), prop.level_end(0)});
    }

    void p_init_vertex_in_class(std::size_t class_index) {
        if (class_index > m_vertex_in_class.size())
            throw std::logic_error("Error: class_index out of range!");
        if (class_index == m_vertex_in_class.size()) {
            m_vertex_in_class.emplace_back();
        }
        std::vector<LitOrFixed>& v_in_c = m_vertex_in_class[class_index];
        auto one_fixed = [&](bool f1, SatLit l2) {
            if (!f1) {
                v_in_c.emplace_back(std::in_place_index<1>, false);
            } else {
                v_in_c.emplace_back(std::in_place_index<0>, l2);
            }
        };
        for (std::size_t vi_i :
             range(v_in_c.size(), m_explicit_cover_order.size()))
        {
            std::size_t vi = m_explicit_cover_order[vi_i];
            Vertex v = m_subgraph->vertex(vi);
            LitOrFixed val1 = p_value_in_class(v.first, class_index);
            LitOrFixed val2 = p_value_in_class(v.second, class_index);
            std::visit(
                overloaded{
                    [&](bool& f1, bool& f2) {
                        v_in_c.emplace_back(std::in_place_index<1>, f1 && f2);
                    },
                    [&](SatLit& l1, SatLit& l2) {
                        SatLit var = m_incremental_solver.new_var();
                        m_incremental_solver.add_short_clause(-var, l1);
                        m_incremental_solver.add_short_clause(-var, l2);
                        m_incremental_solver.add_short_clause(var, -l1, -l2);
                        v_in_c.emplace_back(std::in_place_index<0>, var);
                    },
                    [&](bool& f1, SatLit& l2) { one_fixed(f1, l2); },
                    [&](SatLit& l1, bool& f2) { one_fixed(f2, l1); }},
                val1, val2);
        }
    }

    void p_init_class_clique(std::size_t class_index) {
        p_init_class_literals_clique(class_index);
        p_init_vertex_in_class(class_index);
    }

    void p_init_class_not_clique(std::size_t class_index) {
        p_init_class_literals_not_clique(class_index);
        p_init_vertex_in_class(class_index);
    }

    void p_init_class_literals() {
        for (std::size_t i : range(m_clique_size)) {
            p_init_class_clique(i);
        }
        for (std::size_t i : range(m_clique_size, m_lower_bound)) {
            p_init_class_not_clique(i);
        }
    }

    void p_setup_vertex_constraints(SatLit or_not_enough,
                                    std::size_t begin_vertex,
                                    std::size_t end_vertex) {
        auto is_fixed = [](const LitOrFixed& f) { return f.index() == 1; };
        for (std::size_t vi_i : range(begin_vertex, end_vertex)) {
            bool found_true = false;
            m_clause_buffer.clear();
            for (std::size_t cindex : range(m_lower_bound)) {
                const LitOrFixed& f = m_vertex_in_class[cindex][vi_i];
                if (is_fixed(f)) {
                    if (std::get<1>(f)) {
                        found_true = true;
                        break;
                    }
                } else {
                    m_clause_buffer.push_back(std::get<0>(f));
                }
            }
            if (found_true)
                continue;
            m_clause_buffer.push_back(or_not_enough);
            m_incremental_solver.add_clause(m_clause_buffer.begin(),
                                            m_clause_buffer.end());
        }
    }

    void p_setup_vertex_constraints(std::size_t old_n) {
        std::size_t not_enough_idx = m_lower_bound - m_clique_size;
        std::size_t new_n = m_vertex_in_class[0].size();
        while (not_enough_idx > m_not_enough_colors.size()) {
            m_not_enough_colors.emplace_back(std::in_place_index<1>, true);
        }
        if (not_enough_idx == m_not_enough_colors.size()) {
            SatLit var = m_incremental_solver.new_var();
            m_not_enough_colors.emplace_back(std::in_place_index<0>, var);
            p_setup_vertex_constraints(var, 0, new_n);
        } else {
            p_setup_vertex_constraints(
                std::get<0>(m_not_enough_colors[not_enough_idx]), old_n, new_n);
        }
    }

    void p_add_class() {
        for (VertexInfo& info : m_vertex_info) {
            info.open_classes.push_back(true);
            ++info.num_open_classes;
        }
        m_color_classes.emplace_back(this, m_subgraph->get_propagator(),
                                     m_lower_bound - 1);
    }

    std::size_t p_update_classes() {
        std::size_t old_num_vertices = 0;
        if (m_class_literals.empty()) {
            p_init_class_literals();
        } else {
            old_num_vertices = m_vertex_in_class[0].size();
            for (std::size_t cindex : range(m_class_literals.size())) {
                p_init_vertex_in_class(cindex);
            }
            for (std::size_t cindex :
                 range(m_class_literals.size(), m_lower_bound))
            {
                p_init_class_not_clique(cindex);
            }
        }
        return old_num_vertices;
    }

    void p_no_class_available() {
        m_explicit_cover_order.insert(m_explicit_cover_order.end(),
                                      m_true_zero_vertices.begin(),
                                      m_true_zero_vertices.end());
        for (;;) {
            if (m_aborted.load())
                return;
            std::size_t old_num_vertices = p_update_classes();
            p_setup_vertex_constraints(old_num_vertices);
            std::vector<SatLit> assumptions;
            assumptions.push_back(-std::get<0>(m_not_enough_colors.back()));
            if (m_aborted.load())
                return;
            auto result = m_incremental_solver.solve(assumptions);
            if (!result) {
                m_aborted.store(true);
            }
            if (m_aborted.load())
                return;
            if (*result) {
                // SAT
                const auto& model = m_incremental_solver.get_model();
                p_handle_sat_assignment(model);
                return;
            } else {
                // UNSAT - repeat with more colors
                m_incremental_solver.fix(
                    std::get<0>(m_not_enough_colors.back()));
                if (++m_lower_bound == m_upper_bound)
                    return;
                p_add_class();
            }
        }
    }

    template <typename ModelType>
    bool p_vertex_in_class(const ModelType& model, std::size_t class_index,
                           std::size_t vi_i) {
        LitOrFixed& lf = m_vertex_in_class[class_index][vi_i];
        return std::visit(
            overloaded{[&](bool& b) -> bool { return b; },
                       [&](SatLit& l) -> bool { return model[l]; }},
            lf);
    }

    void p_reset_info() {
        for (VertexInfo& v : m_vertex_info) {
            v.open_classes.set();
            v.num_open_classes = m_lower_bound;
            v.in_some_class = false;
        }
    }

    template <typename ModelType>
    void p_handle_sat_assignment(const ModelType& model) {
        p_reset_info();
        DynamicBitset vertex_is_in{m_explicit_cover_order.size(), false};
        for (std::size_t cindex : range(m_lower_bound)) {
            vertex_is_in.reset();
            for (std::size_t vi_i : range(m_explicit_cover_order.size())) {
                if (p_vertex_in_class(model, cindex, vi_i)) {
                    vertex_is_in[vi_i].set();
                }
            }
            m_color_classes[cindex].reset(this, vertex_is_in);
        }
    }

    void p_reset() {
        p_reset_info();
        DynamicBitset explicit_is_in{m_explicit_cover_order.size(), false};
        for (std::size_t cindex : range(m_clique_size)) {
            explicit_is_in[cindex].set();
            m_color_classes[cindex].reset(this, explicit_is_in);
            explicit_is_in[cindex].reset();
        }
        for (std::size_t cindex : range(m_clique_size, m_lower_bound)) {
            m_color_classes[cindex].reset(this, explicit_is_in);
        }
        m_explicit_cover_order.resize(m_clique_size);
        m_incremental_solver = IncrementalSolver{};
        m_class_literals.clear();
        m_not_enough_colors.clear();
        m_vertex_in_class.clear();
    }

    bool p_finalize() {
        struct EmptyHandler {};
        EmptyHandler handler;
        Var n_all = m_subgraph->get_propagator().db().num_vars();
        Var n_concrete = m_subgraph->get_n_concrete();
        ClassCompleter<EmptyHandler> completer{n_concrete, n_all, &handler};
        bool result = true;
        for (std::size_t cindex : range(m_color_classes.size())) {
            SharedDBPropagator& prop = m_color_classes[cindex].propagator;
            if (!completer.complete_class(prop)) {
                result = false;
            }
        }
        if (result)
            return true;
        p_reset();
        return false;
    }

  public:
    SATDSaturCoverer(UniverseSubgraph* subgraph,
                     std::vector<std::size_t> clique, std::size_t upper_bound)
        : m_subgraph(subgraph), m_best_infos(20, CompareByInfo{this}),
          m_clique_size(clique.size()), m_lower_bound(clique.size()),
          m_upper_bound(upper_bound) {
        p_fill_vertices_with_literal();
        p_init_vertex_info();
        p_init_classes(std::move(clique));
    }

    void abort() {
        m_aborted.store(true);
        m_incremental_solver.terminate();
    }

    void run_coloring(bool until_complete = true) {
        for (;;) {
            if (m_aborted.load())
                return;
            auto n = p_find_next_vertex();
            if (!n) {
                if (!until_complete || p_finalize()) {
                    m_success = true;
                    return;
                }
                continue;
            }
            const VertexInfo& next = m_vertex_info[*n];
            if (next.num_open_classes == 0) {
                p_no_class_available();
                if (m_lower_bound == m_upper_bound) {
                    m_success = true;
                    return;
                }
            } else {
                for (std::size_t cindex : next.open_classes.ones()) {
                    m_color_classes[cindex].add_vertex(this, *n);
                    m_explicit_cover_order.push_back(*n);
                    break;
                }
            }
        }
    }

    bool was_successful() const { return m_success; }

    bool did_improve() const { return m_lower_bound < m_upper_bound; }

    std::size_t get_lower_bound() const noexcept { return m_lower_bound; }

    std::vector<SharedDBPropagator> get_incomplete_solution() const {
        std::vector<SharedDBPropagator> props;
        props.reserve(m_color_classes.size());
        std::transform(m_color_classes.begin(), m_color_classes.end(),
                       std::back_inserter(props),
                       [](const ColorClass& cc) { return cc.propagator; });
        return props;
    }

    PartialSolution get_partial_solution() {
        if (!was_successful() || !did_improve()) {
            throw std::logic_error("Tried to obtain partial solution from "
                                   "unsucessful or non-improving coverer!");
        }
        PartialSolution solution{m_subgraph->get_propagator().db().num_vars(),
                                 m_subgraph->get_infeasibility_map()};
        for (const auto& cc : m_color_classes) {
            solution.add_class(cc.propagator);
        }
        return solution;
    }
};

} // namespace sammy

#endif
==> ./logging.h <==
#ifndef SAMMY_LOGGING_H_INCLUDED_
#define SAMMY_LOGGING_H_INCLUDED_

#include "literals.h"

namespace sammy {

template <typename StreamType, typename ContainerType>
inline void print_clause_internal(StreamType& stream,
                                  const ContainerType& clause) {
    bool first = true;
    for (Lit l : clause) {
        if (!first)
            stream << ' ';
        stream << l;
        first = false;
    }
}

template <typename StreamType, typename ContainerType>
inline void print_clause_external(StreamType& stream,
                                  const ContainerType& clause) {
    bool first = true;
    for (Lit l : clause) {
        if (!first)
            stream << ' ';
        stream << lit::externalize(l);
        first = false;
    }
}

} // namespace sammy

#endif
==> ./stamp_set.h <==
#ifndef SAMMY_STAMP_SET_H_INCLUDED_
#define SAMMY_STAMP_SET_H_INCLUDED_

#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>

namespace sammy {

template <typename ValueType, typename StampType = std::uint32_t>
class StampSet {
  private:
    static_assert(std::is_integral_v<ValueType>, "Value type must be integral");
    static_assert(std::is_integral_v<StampType> &&
                      std::is_unsigned_v<StampType>,
                  "Stamp type must be unsigned");

    std::vector<StampType> m_stamps;
    StampType m_current_stamp;

  public:
    explicit StampSet(ValueType universe_size)
        : m_stamps(universe_size, StampType(0)), m_current_stamp(1) {}

    StampSet(const StampSet&) = default;
    StampSet& operator=(const StampSet&) = default;
    StampSet(StampSet&&) noexcept = default;
    StampSet& operator=(StampSet&&) noexcept = default;

    std::size_t universe_size() const noexcept { return m_stamps.size(); }

    void clear() noexcept {
        if (++m_current_stamp == 0) {
            std::fill(m_stamps.begin(), m_stamps.end(), StampType(0));
            m_current_stamp = 1;
        }
    }

    template <typename ForwardIterator>
    void assign(ForwardIterator begin, ForwardIterator end) noexcept {
        clear();
        insert(begin, end);
    }

    template <typename ForwardIterator>
    void insert(ForwardIterator begin, ForwardIterator end) noexcept {
        std::for_each(begin, end, [&](ValueType l) { insert(l); });
    }

    void insert(ValueType v) noexcept { m_stamps[v] = m_current_stamp; }

    void erase(ValueType v) noexcept { m_stamps[v] = 0; }

    bool check_insert(ValueType v) noexcept {
        StampType& s = m_stamps[v];
        bool result = (s != m_current_stamp);
        s = m_current_stamp;
        return result;
    }

    bool check_erase(ValueType v) noexcept {
        StampType& s = m_stamps[v];
        bool result = (s == m_current_stamp);
        s = 0;
        return result;
    }

    bool count(ValueType v) const noexcept {
        return m_stamps[v] == m_current_stamp;
    }

    bool contains(ValueType v) const noexcept { return count(v); }
};

} // namespace sammy

#endif
==> ./thread_clauses.h <==
#ifndef SAMMY_THREAD_CLAUSES_H_INCLUDED_
#define SAMMY_THREAD_CLAUSES_H_INCLUDED_

#include "clause_db.h"
#include "literals.h"
#include <mutex>
#include <optional>

namespace sammy {

template <typename T, typename Tag = void> class ThreadLocalManager {
  public:
    using Ticket = std::size_t;

    /**
     * Create a new 'slot' for a thread-local object of type T,
     * identified by the given Ticket value.
     * Threads can access a thread-local mutable copy of the object constructed
     * by the given arguments using the returned Ticket.
     */
    template <typename... Args>
    Ticket new_ticket(Args&&... unshared_object_args) {
        std::unique_lock l{m_mutex};
        Ticket result = m_original_objects.size();
        m_original_objects.emplace_back(
            std::forward<Args>(unshared_object_args)...);
        return result;
    }

    T& local(Ticket ticket) { return p_local(ticket); }

    const T& local(Ticket ticket) const { return p_local(ticket); }

  private:
    T& p_local(Ticket ticket) const {
        thread_local std::vector<std::optional<T>> objects;
        while (ticket > objects.size()) {
            objects.emplace_back(std::nullopt);
        }
        if (ticket == objects.size()) {
            std::unique_lock l{m_mutex};
            objects.emplace_back(m_original_objects.at(ticket));
        } else if (!objects[ticket]) {
            std::unique_lock l{m_mutex};
            objects[ticket] = m_original_objects.at(ticket);
        }
        return *objects[ticket];
    }

    mutable std::mutex m_mutex;
    std::vector<T> m_original_objects;
};

namespace detail {
inline ThreadLocalManager<ClauseDB>& get_tl_clause_manager() {
    static ThreadLocalManager<ClauseDB> manager;
    return manager;
}
} // namespace detail

using ClausesTicket = ThreadLocalManager<ClauseDB>::Ticket;

inline ClausesTicket publish_clauses(const ClauseDB& clause_db) {
    return detail::get_tl_clause_manager().new_ticket(clause_db);
}

inline ClauseDB& local_clauses(ClausesTicket ticket) {
    return detail::get_tl_clause_manager().local(ticket);
}

} // namespace sammy

#endif
==> ./primal_dual_driver.h <==
#ifndef SAMMY_PRIMAL_DUAL_DRIVER_H_INCLUDED_
#define SAMMY_PRIMAL_DUAL_DRIVER_H_INCLUDED_

#include "algorithm_ex.h"
#include "coloring.h"
#include "experiment_flags.h"
#include "fast_clique.h"
#include "gurobi_clique_solver_g2.h"
#include "initial_coloring_heuristic.h"
#include "output.h"
#include "pair_infeasibility_map.h"
#include "simplification.h"
#include "universe_subgraph.h"

namespace sammy {

static constexpr std::size_t MAX_CUTS_FROM_PRIMAL = 5;
static constexpr double MIN_VERTICES_PER_PRICING = 10.0;
static constexpr double MIN_RELATIVE_PER_PRICING = 0.01;
static constexpr std::size_t MAX_CUT_ROUNDS_PER_PRICING = 40;
static constexpr double PRIMAL_TO_DUAL_GOAL_RATIO = 0.1;

class PossiblySimplified {
  public:
    PossiblySimplified(const ClauseDB& original, Var n_concrete)
        : m_original(original), m_original_concrete(n_concrete) {}

    ClauseDB& formula() noexcept {
        return m_simplified ? m_simplified->formula : m_original;
    }

    const ClauseDB& formula() const noexcept {
        return m_simplified ? m_simplified->formula : m_original;
    }

    Var n_concrete() const noexcept {
        return m_simplified ? m_simplified->num_concrete : m_original_concrete;
    }

    bool is_simplified() const noexcept {
        return static_cast<bool>(m_simplified);
    }

    void simplify(OutputObject& output) {
        m_simplification_stats.capture_before(m_original, m_original_concrete);
        m_simplifier.emplace(m_original, m_original_concrete);
        m_simplified.emplace(
            sammy::run_simplifier(*m_simplifier, m_simplification_stats));
        m_simplification_stats.capture_after(formula(), n_concrete());
        export_simplification(output, *m_simplifier, *m_simplified,
                              &statistics());
    }

    const SimplificationStats& statistics() const noexcept {
        return m_simplification_stats;
    }

    std::vector<std::vector<bool>>
    reconstruct_sample(const std::vector<std::vector<bool>>& s) const {
        if (!m_simplified)
            return s;
        return m_simplifier->reconstruct_sample(*m_simplified, s);
    }

    std::vector<Vertex> reconstruct_mes(const std::vector<Vertex>& mes) const {
        if (!m_simplified)
            return mes;
        return m_simplifier->reconstruct_lb(*m_simplified, mes);
    }

  private:
    ClauseDB m_original;
    Var m_original_concrete;
    std::optional<SimplifyDatastructure> m_simplifier;
    std::optional<SimplifiedInstance> m_simplified;
    SimplificationStats m_simplification_stats;
};

class PrimalDualDriver {
  public:
    PrimalDualDriver(EventRecorder& recorder, OutputObject& output,
                     const ExperimentFlagsConfig& config,
                     const ClauseDB& formula, Var n_concrete)
        : m_begin_time(Clock::now()), m_config(config), m_recorder(&recorder),
          m_output(&output), m_thread_pool(), m_instance(formula, n_concrete),
          m_inf_map((p_possibly_simplify(), m_instance.n_concrete())),
          m_fast_col_solver(&m_instance.formula(), m_instance.n_concrete(),
                            &m_thread_pool, &m_inf_map),
          m_fast_clq_solver(SharedDBPropagator(&m_instance.formula()),
                            &m_thread_pool),
          m_state(nullptr), m_last_primal_ended_at(m_begin_time) {
        m_fast_col_solver.set_quiet(true);
    }

    void run() {
        if (!m_state) {
            m_state = p_initial_state();
        }
        while (m_state) {
            m_state->main();
            auto next = m_state->leave();
            if (!next) {
                m_state.reset();
                break;
            } else {
                next->enter();
                m_state = std::move(next);
                if (m_state->is_final())
                    break;
            }
        }
    }

    bool mes_is_optimal() const noexcept {
        if (!m_best_mes || !m_best_sample)
            return false;
        if (m_best_mes->size() >= m_best_sample->size())
            return true;
        if (!m_state)
            return false;
        auto* last_state = m_state.get();
        if (dynamic_cast<LBOptimalState*>(last_state)) {
            return true;
        }
        return false;
    }

    std::size_t get_best_lb() const noexcept {
        std::size_t bound = 0;
        if (get_best_mes()) {
            bound = get_best_mes()->size();
        }
        if (m_best_sat_bound > bound) {
            bound = m_best_sat_bound;
        }
        return bound;
    }

    const std::optional<std::vector<Vertex>>& get_best_mes() const {
        return m_best_mes;
    }

    const std::optional<std::vector<std::vector<bool>>>&
    get_best_sample() const {
        return m_best_sample;
    }

    class Timeout : public std::exception {
        const char* what() const noexcept override { return "timeout"; }
    };

    bool update_best_sample(const std::vector<std::vector<bool>>& sample,
                            const char* source) {
        if (!m_best_sample || sample.size() < m_best_sample->size()) {
            m_recorder->store_event(
                "IMPROVED_SAMPLE",
                {{"size", sample.size()}, {"source", source}}, "size",
                "source");
            m_best_sample = sample;
            return true;
        }
        return false;
    }

    bool update_best_sample(const std::vector<SharedDBPropagator>& propagators,
                            const char* source) {
        if (m_best_sample && m_best_sample->size() <= propagators.size())
            return false;
        std::vector<std::vector<bool>> s;
        s.reserve(propagators.size());
        std::transform(propagators.begin(), propagators.end(),
                       std::back_inserter(s),
                       [](const SharedDBPropagator& prop) {
                           return prop.extract_assignment();
                       });
        return update_best_sample(s, source);
    }

    bool update_best_mes(const std::vector<Vertex>& mes, const char* source) {
        if (!m_best_mes || mes.size() > m_best_mes->size()) {
            m_recorder->store_event("IMPROVED_MUTUALLY_EXCLUSIVE_SET",
                                    {{"size", mes.size()},
                                     {"source", source},
                                     {"simplified_vertices", mes}},
                                    "size", "source");
            m_best_mes = mes;
            return true;
        }
        return false;
    }

    template <typename InputIterator>
    void update_spawners(InputIterator begin, InputIterator end) {
        auto old_size = m_all_spawners_set.size();
        m_all_spawners_set.insert(begin, end);
        if (m_all_spawners_set.size() != old_size) {
            m_all_spawners.assign(m_all_spawners_set.begin(),
                                  m_all_spawners_set.end());
        }
    }

    template <typename Container>
    void update_spawners(const Container& container) {
        using std::begin;
        using std::end;
        update_spawners(begin(container), end(container));
    }

    ClauseDB& formula() noexcept { return m_instance.formula(); }

    const ClauseDB& formula() const noexcept { return m_instance.formula(); }

    const std::vector<Vertex>& valid_interactions() {
        if (!m_valid_interactions) {
            m_recorder->store_event("BEGIN_COLLECTING_VALID_INTERACTIONS");
            m_valid_interactions.emplace(m_inf_map.collect_vertices(1'000'000));
            m_recorder->store_event("END_COLLECTING_VALID_INTERACTIONS",
                                    {{"count", m_valid_interactions->size()}},
                                    "count");
        }
        return *m_valid_interactions;
    }

    double time_remaining() const {
        return m_config.lb_mip_config.total_mip_timeout -
               seconds_between(m_begin_time, Clock::now());
    }

    bool rerun_primal_with_dual_info() {
        Clock::time_point before_prun = Clock::now();
        auto& lbsolver = *m_cnp_clq_solver;
        auto& ubsolver = m_fast_col_solver;
        const auto& best_mes = *m_best_mes;
        auto& rec = *m_recorder;
        double fraction =
            std::uniform_real_distribution<double>{0.0, 1.0}(sammy::rng());
        std::size_t max_num_constraints =
            std::size_t(fraction * (best_mes.size() - 1));
        rec.store_event("BEGIN_PRIMAL_WITH_DUAL_INFORMATION",
                        {{"max_configs_from_dual", max_num_constraints}},
                        "max_configs_from_dual");
        ubsolver.reset_coloring();
        lbsolver.export_highest_weights_to_primal(ubsolver,
                                                  max_num_constraints);
        std::vector<Vertex> vinitial = lbsolver.get_fractional_mes_support();
        vinitial.insert(vinitial.end(), best_mes.begin(), best_mes.end());
        std::sort(vinitial.begin(), vinitial.end());
        vinitial.erase(std::unique(vinitial.begin(), vinitial.end()),
                       vinitial.end());
        ubsolver.color_lazy(vinitial);
        update_spawners(ubsolver.class_spawners());
        update_best_sample(ubsolver.all_classes(),
                           "dual-initialized heuristic");
        std::size_t cuts_extracted = lbsolver.cuts_from_primal(
            ubsolver.all_classes(), MAX_CUTS_FROM_PRIMAL);
        rec.store_event("END_PRIMAL_WITH_DUAL_INFORMATION",
                        {{"solution_found", ubsolver.all_classes().size()},
                         {"best_known_solution", m_best_sample->size()},
                         {"cuts_extracted", cuts_extracted}},
                        "solution_found", "best_known_solution",
                        "cuts_extracted");
        Clock::time_point after_prun = Clock::now();
        m_last_primal_time = seconds_between(before_prun, after_prun);
        m_last_primal_ended_at = after_prun;
        return cuts_extracted > 0;
    }

  private:
    ClauseDB& p_possibly_simplify() {
        if (!m_config.dont_simplify && !m_instance.is_simplified()) {
            m_instance.simplify(*m_output);
        }
        return m_instance.formula();
    }

    class DriverState {
      public:
        explicit DriverState(PrimalDualDriver* driver)
            : driver(driver), successor(nullptr) {}

        void enter() { p_enter(); }

        virtual ~DriverState() = default;
        virtual void main() = 0;

        virtual bool is_final() const { return false; }

        virtual std::unique_ptr<DriverState> leave() {
            p_leave();
            std::unique_ptr<DriverState> result = std::move(successor);
            return result;
        }

      protected:
        PrimalDualDriver* driver;
        std::unique_ptr<DriverState> successor;

        // called with old state in driver->m_state
        virtual void on_enter() {}
        // called with current state in driver->m_state
        virtual void on_leave() {}

      private:
        void p_enter() {
            p_check_time();
            on_enter();
        }

        void p_leave() { on_leave(); }

        void p_check_time() {
            if (!std::isfinite(
                    driver->m_config.lb_mip_config.total_mip_timeout))
            {
                return;
            }
            auto tnow = seconds_between(driver->m_begin_time, Clock::now());
            if (tnow >= driver->m_config.lb_mip_config.total_mip_timeout) {
                throw Timeout();
            }
        }
    };

    /**
     * The typical initial state, where we run
     * the initial heuristic for the first time.
     * This requires some extra actions compared
     * to follow-up runs of the heuristic.
     */
    class InitialState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            auto& driver = *this->driver;
            auto& rec = *driver.m_recorder;
            auto& fcol = driver.m_fast_col_solver;
            auto& fclq = driver.m_fast_clq_solver;
            rec.store_event("BEGIN_FIRST_SOLVE");
            fcol.initialize_feasibilities();
            fcol.color_lazy();
            rec.store_event("DONE_FIRST_SOLVE");
            driver.update_best_sample(fcol.all_classes(),
                                      "initial coloring heuristic");
            rec.store_event("BEGIN_FEASIBILITY_EXTRACTION");
            fcol.extract_feasibilities();
            rec.store_event("DONE_FEASIBILITY_EXTRACTION");
            rec.store_event("BEGIN_LEARN_INFEASIBILITIES");
            learn_infeasibilities(driver.formula(), &driver.m_inf_map);
            rec.store_event("DONE_LEARN_INFEASIBILITIES");
            driver.update_spawners(fcol.class_spawners());
            rec.store_event("BEGIN_FIRST_LB");
            auto r = driver.m_config.initial_heuristic_config
                         .random_clique_restarts_per_iteration;
            auto best_mes =
                fclq.random_multistart_best_clique(r, driver.m_all_spawners);
            driver.update_best_mes(best_mes, "initial clique heuristic");
            rec.store_event("DONE_FIRST_LB");
            ++driver.m_initial_iterations;
        }

        std::unique_ptr<DriverState> leave() override {
            return std::make_unique<InitialHeuristicState>(driver);
        }
    };

    /**
     * State where we repeat runs of the initial heuristic.
     */
    class InitialHeuristicState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            auto& driver = *this->driver;
            auto& rec = *driver.m_recorder;
            auto& fcol = driver.m_fast_col_solver;
            auto& fclq = driver.m_fast_clq_solver;
            ++driver.m_initial_iterations;
            rec.store_event("BEGIN_NEXT_HEURISTIC_SAMPLE",
                            {{"iteration", driver.m_initial_iterations}},
                            "iteration");
            fcol.reset_coloring();
            fcol.color_lazy(*driver.m_best_mes);
            driver.update_best_sample(fcol.all_classes(),
                                      "repeat coloring heuristic");
            driver.update_spawners(fcol.class_spawners());
            rec.store_event("DONE_NEXT_HEURISTIC_SAMPLE",
                            {{"iteration", driver.m_initial_iterations}},
                            "iteration");
            rec.store_event("BEGIN_NEXT_HEURISTIC_MES",
                            {{"iteration", driver.m_initial_iterations}},
                            "iteration");
            auto r = driver.m_config.initial_heuristic_config
                         .random_clique_restarts_per_iteration;
            driver.update_best_mes(fclq.random_multistart_best_clique(
                                       r / 2, driver.m_all_spawners),
                                   "repeat clique heuristic (all spawners)");
            driver.update_best_mes(fclq.random_multistart_best_clique(
                                       r / 2, fcol.class_spawners()),
                                   "repeat clique heuristic (new spawners)");
            rec.store_event("DONE_NEXT_HEURISTIC_MES",
                            {{"iteration", driver.m_initial_iterations}},
                            "iteration");
        }

        std::unique_ptr<DriverState> leave() override {
            const double time_needed =
                seconds_between(driver->m_begin_time, Clock::now());
            auto goal_iter =
                driver->m_config.initial_heuristic_config.goal_iterations;
            if (time_needed >=
                    driver->m_config.initial_heuristic_config.max_time ||
                (time_needed >=
                     driver->m_config.initial_heuristic_config.min_time &&
                 driver->m_initial_iterations >= goal_iter))
            {
                if (driver->m_config.sat_coloring_config.initial_sat_coloring) {
                    return std::make_unique<InitialSATColoringState>(driver);
                }
                return std::make_unique<PrimalDualInitializationState>(driver);
            }
            return std::make_unique<InitialHeuristicState>(driver);
        }
    };

    /**
     * If we want to initialize using SAT coloring,
     * we enter this state.
     */
    class InitialSATColoringState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            tbegin = Clock::now();
            successor = std::make_unique<PrimalDualInitializationState>(driver);
            auto& rec = *driver->m_recorder;
            auto& ubsolver = driver->m_fast_col_solver;
            auto& formula = driver->formula();
            std::vector<Vertex> initial_vertices(*driver->get_best_mes());
            add_and_make_unique(initial_vertices,
                                ubsolver.class_spawners().begin(),
                                ubsolver.class_spawners().end());
            UniverseSubgraph subgraph{&formula, &driver->m_thread_pool,
                                      &driver->m_inf_map, initial_vertices};
            subgraph.extend_matrix_by_propagation();
            p_update_clique_indices(*driver->get_best_mes(), subgraph);
            extra_constraints.update_vertex_count(subgraph);
            GurobiCliqueSolverG2 g2_clique_solver{
                &subgraph,
                formula.num_vars(),
                &rec,
                *driver->get_best_sample(),
                *driver->get_best_mes(),
                driver->m_config.lb_mip_config};
            p_run_satcoloring(subgraph, g2_clique_solver);
        }

      private:
        void p_run_satcoloring(UniverseSubgraph& subgraph,
                               GurobiCliqueSolverG2& clique_solver) {
            std::size_t best_bound = driver->get_best_mes()->size();
            std::size_t best_soln = driver->get_best_sample()->size();
            auto& recorder = *driver->m_recorder;
            auto& ubsolver = driver->m_fast_col_solver;
            while (best_bound < best_soln) {
                double trem = driver->m_config.sat_coloring_config
                                  .initial_sat_coloring_timeout;
                trem -= seconds_between(tbegin, Clock::now());
                if (trem <= 0) {
                    // timeout for initial sat coloring;
                    // return (successor already set)
                    driver->m_best_sat_bound = best_bound;
                    return;
                }
                SATKColoringSolver ksolve{&subgraph, driver->m_recorder,
                                          cindices, best_bound};
                ksolve.set_extra_constraints(&extra_constraints);
                auto ksolve_res = ksolve.solve(trem);
                if (!ksolve_res) {
                    // timeout for initial sat coloring;
                    // return (successor already set)
                    driver->m_best_sat_bound = best_bound;
                    return;
                }
                if (!*ksolve_res) {
                    best_bound += 1;
                    if (p_optimize_on_subgraph(clique_solver)) {
                        driver->update_best_mes(
                            clique_solver.get_best_mes(),
                            "LP rounding during initial SAT coloring");
                        const auto& mes = *driver->get_best_mes();
                        p_update_clique_indices(mes, subgraph);
                        if (mes.size() > best_bound)
                            best_bound = mes.size();
                    }
                } else {
                    std::vector<SharedDBPropagator> classes =
                        extra_constraints.coloring_to_classes(
                            subgraph.get_propagator(), subgraph,
                            ksolve.get_coloring(), best_bound);
                    if (classes.empty()) {
                        recorder.store_event("ADDED_CONSTRAINTS");
                        continue; // new extra constraints!
                    }
                    ubsolver.reset_coloring();
                    for (SharedDBPropagator& c : classes) {
                        ubsolver.add_color_class(std::move(c));
                    }
                    ubsolver.color_lazy(driver->m_all_spawners);
                    if (driver->update_best_sample(
                            ubsolver.all_classes(),
                            "extension of initial SAT coloring"))
                    {
                        best_soln = driver->get_best_sample()->size();
                    }
                    if (best_soln != best_bound) {
                        const auto& new_spawners = ubsolver.class_spawners();
                        clique_solver.add_new_vertices(new_spawners);
                        extra_constraints.update_vertex_count(subgraph);
                    }
                }
            }
        }

        bool p_optimize_on_subgraph(GurobiCliqueSolverG2& clique_solver) {
            std::size_t clique_size_before = driver->get_best_mes()->size();
            for (;;) {
                double trem = driver->m_config.sat_coloring_config
                                  .initial_sat_coloring_timeout;
                trem -= seconds_between(tbegin, Clock::now());
                switch (clique_solver.solve_full_relaxation(trem)) {
                case SolverState::OPTIMUM_ON_SUBGRAPH:
                case SolverState::TIMEOUT_IMPROVEMENT:
                case SolverState::TIMEOUT_NO_IMPROVEMENT:
                    return clique_size_before <
                           clique_solver.get_best_mes().size();

                default:
                    if (clique_solver.greedy_add_to_cutting_planes())
                        continue;
                    if (clique_solver.greedy_generate_cutting_planes())
                        continue;
                    return clique_size_before <
                           clique_solver.get_best_mes().size();
                }
            }
        }

        void p_update_clique_indices(const std::vector<Vertex>& clique,
                                     const UniverseSubgraph& subgraph) {
            cindices.clear();
            std::transform(clique.begin(), clique.end(),
                           std::back_inserter(cindices),
                           [&](Vertex v) { return subgraph.vertex_index(v); });
        }

        Clock::time_point tbegin;
        std::vector<std::size_t> cindices;
        ColoringExtraConstraints extra_constraints;
    };

    /**
     * State in which there is some solution and lower bound
     * from the initial heuristic and we initialize the primal/dual
     * information, build the subgraph, etc.
     */
    class PrimalDualInitializationState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            auto& rec = *driver->m_recorder;
            rec.store_event("BEGIN_INITIAL_SUBGRAPH_CREATION");
            p_create_subgraph();
            auto& sg = *driver->m_subgraph;
            auto& formula = driver->formula();
            rec.store_event("END_INITIAL_SUBGRAPH_CREATION", {{"n", sg.n()}},
                            "n");
            rec.store_event("BEGIN_INITIAL_MODEL_CREATION");
            driver->m_cnp_clq_solver.emplace(
                &sg, formula.num_vars(), driver->m_recorder,
                *driver->m_best_sample, *driver->m_best_mes,
                driver->m_config.lb_mip_config);
            rec.store_event("END_INITIAL_MODEL_CREATION");
        }

        std::unique_ptr<DriverState> leave() override {
            return std::make_unique<PrimalDualBaseState>(driver);
        }

      private:
        std::vector<std::size_t>
        p_extract_coloring(const std::vector<std::vector<bool>>& sample,
                           const std::vector<Vertex>& vertices) {
            const auto n = vertices.size(), nc = sample.size();
            std::vector<std::size_t> result(n, 0);
            for (std::size_t i = 0; i < n; ++i) {
                Vertex v = vertices[i];
                for (std::size_t k = 0; k < nc; ++k) {
                    const std::vector<bool>& config = sample[k];
                    if (config[v.first] && config[v.second]) {
                        result[i] = k;
                        break;
                    }
                }
            }
            return result;
        }

        void p_run_satcoloring() {
            std::vector<Vertex> verts =
                driver->m_fast_col_solver.class_spawners();
            verts.insert(verts.end(), driver->m_best_mes->begin(),
                         driver->m_best_mes->end());
            std::sort(verts.begin(), verts.end());
            verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
            UniverseSubgraph subgraph(&driver->formula(),
                                      &driver->m_thread_pool,
                                      &driver->m_inf_map, verts);
            subgraph.extend_matrix_by_propagation();
            SATColoringSolver solver{
                &subgraph, driver->m_recorder, *driver->m_best_mes,
                p_extract_coloring(*driver->m_best_sample, verts)};
            solver.optimize_coloring();
            auto [classes, spawners] = solver.coloring_to_classes();
            driver->m_recorder->store_event("CLASSES_FROM_COLORING",
                                            {{"num_classes", classes.size()}},
                                            "num_classes");
            auto& fcol = driver->m_fast_col_solver;
            fcol.reset_coloring();
            for (std::size_t ci = 0, cn = classes.size(); ci < cn; ++ci) {
                fcol.add_color_class(std::move(classes[ci]), spawners[ci]);
            }
            classes.clear();
            spawners.clear();
            fcol.color_lazy(driver->m_all_spawners);
            driver->update_best_sample(fcol.all_classes(),
                                       "extended SAT coloring");
            driver->update_spawners(fcol.class_spawners());
            subgraph.add_vertices(fcol.class_spawners());
            SATColoringSolver solver2{
                &subgraph, driver->m_recorder, *driver->m_best_mes,
                p_extract_coloring(*driver->m_best_sample,
                                   subgraph.vertex_set())};
            solver2.optimize_coloring();
            std::tie(classes, spawners) = solver.coloring_to_classes();
            fcol.reset_coloring();
            for (std::size_t ci = 0, cn = classes.size(); ci < cn; ++ci) {
                fcol.add_color_class(std::move(classes[ci]), spawners[ci]);
            }
            fcol.color_lazy(driver->m_all_spawners);
            driver->m_recorder->store_event(
                "LAZY_COLORING_RESULT",
                {{"num_classes", fcol.all_classes().size()}}, "num_classes");
        }

        void p_create_subgraph() {
            auto& opt_sg = driver->m_subgraph;
            auto& formula = driver->formula();
            std::vector<Vertex> initial = *driver->m_best_mes;
            const auto& last_spawners =
                driver->m_fast_col_solver.class_spawners();
            const auto& all_spawners = driver->m_all_spawners;
            initial.insert(initial.end(), last_spawners.begin(),
                           last_spawners.end());
            std::sort(initial.begin(), initial.end());
            initial.erase(std::unique(initial.begin(), initial.end()),
                          initial.end());
            if (initial.size() < 2000 && all_spawners.size() < 20000) {
                initial.insert(initial.end(), all_spawners.begin(),
                               all_spawners.end());
                std::sort(initial.begin(), initial.end());
                initial.erase(std::unique(initial.begin(), initial.end()),
                              initial.end());
            }
            opt_sg.emplace(&formula, &driver->m_thread_pool, &driver->m_inf_map,
                           std::move(initial));
            opt_sg->extend_matrix_by_propagation();
        }
    };

    /**
     * The base state of the primal/dual approach.
     * We enter this state if the relaxation needs
     * to be solved (or re-solved).
     */
    class PrimalDualBaseState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            auto& driver = *this->driver;
            auto& lbsolver = *driver.m_cnp_clq_solver;
            double trem = driver.time_remaining();
            SolverState sol_state = lbsolver.solve_full_relaxation(trem);
            lbsolver.count_compressible_constraints(50);
            ++driver.m_full_relaxations_solved;
            if (sol_state == SolverState::TIMEOUT_IMPROVEMENT) {
                driver.update_best_mes(lbsolver.get_best_mes(), "LP rounding");
            }
            driver.m_relaxation_history.push_back(lbsolver.get_last_value());
            switch (sol_state) {
            default:
                throw Timeout{};

            case SolverState::IMPROVEMENT_FOUND:
                driver.update_best_mes(lbsolver.get_best_mes(), "LP rounding");
                successor =
                    std::make_unique<RerunPrimalHeuristicState>(&driver);
                return;

            case SolverState::NO_IMPROVEMENT_FOUND:
                if (driver.p_check_primal_goal_ratio()) {
                    successor =
                        std::make_unique<RerunPrimalHeuristicState>(&driver);
                } else {
                    successor =
                        std::make_unique<NotOptimalOnSubgraphState>(&driver);
                }
                return;

            case SolverState::OPTIMUM_ON_SUBGRAPH:
                driver.m_relaxation_history.clear();
                driver.m_have_optimum_on_subgraph = true;
                if (driver.update_best_mes(lbsolver.get_best_mes(),
                                           "LP rounding") ||
                    driver.p_check_primal_goal_ratio())
                {
                    successor =
                        std::make_unique<RerunPrimalHeuristicState>(&driver);
                } else {
                    successor =
                        std::make_unique<OptimalOnSubgraphState>(&driver);
                }
                return;
            }
        }
    };

    class RerunPrimalHeuristicState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            auto& driver = *this->driver;
            bool primal_added_constraints =
                driver.rerun_primal_with_dual_info();
            if (primal_added_constraints) {
                successor = std::make_unique<PrimalDualBaseState>(&driver);
            } else if (!driver.m_have_optimum_on_subgraph) {
                successor =
                    std::make_unique<NotOptimalOnSubgraphState>(&driver);
            } else {
                successor = std::make_unique<OptimalOnSubgraphState>(&driver);
            }
        }
    };

    class NotOptimalOnSubgraphState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            auto& driver = *this->driver;
            if (driver.m_relaxation_history.size() >=
                MAX_CUT_ROUNDS_PER_PRICING)
            {
                p_price_suboptimal();
                driver.m_relaxation_history.clear();
                successor = std::make_unique<PrimalDualBaseState>(&driver);
                return;
            }
            auto& lbsolver = *driver.m_cnp_clq_solver;
            if (lbsolver.greedy_add_to_cutting_planes() ||
                lbsolver.greedy_generate_cutting_planes() ||
                driver.rerun_primal_with_dual_info())
            {
                successor = std::make_unique<PrimalDualBaseState>(&driver);
                return;
            }
            successor = std::make_unique<CutFailedState>(&driver);
        }

      private:
        void p_price_suboptimal() {
            if (p_price_suboptimal(driver->m_fast_col_solver.class_spawners()))
                return;
            if (p_price_suboptimal(driver->m_all_spawners))
                return;
            p_price_suboptimal(
                driver->m_fast_col_solver.extract_coloring_order());
        }

        bool p_price_suboptimal(const std::vector<Vertex>& vs) {
            return driver->m_cnp_clq_solver->price_vertices(vs.begin(),
                                                            vs.end()) > 0;
        }
    };

    /**
     * A state we enter when we could not find any
     * violated cutting plane despite not having the
     * optimal clique on our current subgraph.
     * Currently, pricing is the only option in this case;
     * we could also try an exact separation procedure first.
     */
    class CutFailedState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            const auto& valid = driver->valid_interactions();
            if (!driver->m_cnp_clq_solver->price_vertices(valid.begin(),
                                                          valid.end()))
            {
                successor = std::make_unique<CutAndPriceFailedState>(driver);
            } else {
                driver->m_relaxation_history.clear();
                successor = std::make_unique<PrimalDualBaseState>(driver);
            }
        }
    };

    class CutAndPriceFailedState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            auto& rec = *driver->m_recorder;
            rec.store_event(
                "CUT_AND_PRICE_FAILED",
                {{"best_mes", driver->m_best_mes->size()},
                 {"lp_bound", driver->m_cnp_clq_solver->get_last_value()}},
                "best_mes", "lp_bound");
            successor = std::make_unique<CutAndPriceFailedState>(driver);
        }

        bool is_final() const override { return true; }
    };

    class OptimalOnSubgraphState : public DriverState {
      public:
        using DriverState::DriverState;

        void main() override {
            auto& driver = *this->driver;
            auto& lbsolver = *driver.m_cnp_clq_solver;
            auto& ubsolver = driver.m_fast_col_solver;
            const auto& subgraph = *driver.m_subgraph;
            const auto& spawners = ubsolver.class_spawners();
            double percent_of_graph =
                (std::max)(MIN_VERTICES_PER_PRICING,
                           MIN_RELATIVE_PER_PRICING * subgraph.n());
            std::size_t pos_vertices =
                lbsolver.price_vertices(spawners.begin(), spawners.end());
            if (pos_vertices >= percent_of_graph) {
                p_successor_found_new();
                return;
            }
            const auto& all_spawners = driver.m_all_spawners;
            pos_vertices += lbsolver.price_vertices(all_spawners.begin(),
                                                    all_spawners.end());
            if (pos_vertices >= percent_of_graph) {
                p_successor_found_new();
                return;
            }
            auto coloring_order = ubsolver.extract_coloring_order();
            pos_vertices += lbsolver.price_vertices(coloring_order.begin(),
                                                    coloring_order.end());
            if (pos_vertices >= percent_of_graph) {
                p_successor_found_new();
                return;
            }
            const auto& valid = driver.valid_interactions();
            pos_vertices += lbsolver.price_vertices(valid.begin(), valid.end());
            if (pos_vertices) {
                p_successor_found_new();
                return;
            }
            successor = std::make_unique<LBOptimalState>(&driver);
        }

      private:
        void p_successor_found_new() {
            driver->m_relaxation_history.clear();
            driver->m_have_optimum_on_subgraph = false;
            successor = std::make_unique<PrimalDualBaseState>(driver);
        }
    };

    class LBOptimalState : public DriverState {
      public:
        using DriverState::DriverState;

        bool is_final() const override { return true; }
        void main() override {
            auto& rec = *driver->m_recorder;
            rec.store_event(
                "OPTIMUM_CLIQUE_ON_G2",
                {{"lb", driver->m_best_mes->size()},
                 {"ub", driver->m_best_sample->size()},
                 {"lp_bound", driver->m_cnp_clq_solver->get_last_value()}},
                "lb", "ub", "lp_bound");
            successor = std::make_unique<LBOptimalState>(driver);
        }
    };

    std::unique_ptr<DriverState> p_initial_state() {
        std::unique_ptr<DriverState> result;
        if (!m_best_sample || !m_best_mes) {
            result = std::make_unique<InitialState>(this);
        }
        // TODO: state on re-run?
        result->enter();
        return result;
    }

    bool p_check_primal_goal_ratio() const {
        auto now = Clock::now();
        double tlast = m_last_primal_time;
        double tsince = seconds_between(m_last_primal_ended_at, now);
        double ftest = 1.0 / PRIMAL_TO_DUAL_GOAL_RATIO - 1.0;
        return tsince >= tlast * ftest;
    }

    Clock::time_point m_begin_time;
    ExperimentFlagsConfig m_config;
    EventRecorder* m_recorder;
    OutputObject* m_output;
    ThreadGroup<void> m_thread_pool;
    PossiblySimplified m_instance;
    PairInfeasibilityMap m_inf_map;
    ColoringHeuristicSolver m_fast_col_solver;
    ParallelFastCliqueBuilder m_fast_clq_solver;
    std::unique_ptr<DriverState> m_state;
    EdgeSet m_all_spawners_set;
    std::vector<Vertex> m_all_spawners;
    std::optional<std::vector<std::vector<bool>>> m_best_sample;
    std::optional<std::vector<Vertex>> m_best_mes;
    std::optional<UniverseSubgraph> m_subgraph;
    std::optional<GurobiCliqueSolverG2> m_cnp_clq_solver;
    std::optional<std::vector<Vertex>> m_valid_interactions;
    std::deque<double> m_relaxation_history;
    std::size_t m_initial_iterations = 0;
    std::size_t m_full_relaxations_solved = 0;
    std::size_t m_best_sat_bound = 0;
    double m_last_primal_time = 0.0;
    Clock::time_point m_last_primal_ended_at;
    bool m_have_optimum_on_subgraph = false;
};

} // namespace sammy

#endif
==> ./rng.h <==
#ifndef HS_RNG_H_INCLUDED_
#define HS_RNG_H_INCLUDED_

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <random>

namespace sammy {

inline std::mt19937_64& rng() noexcept {
    static std::atomic<std::int32_t> counter(0);
    thread_local std::mt19937_64 res(1337 + counter++);
    return res;
}

template<typename ContainerIn, typename ContainerOut, typename RNG>
inline void sample_from_range(const ContainerIn& container_in,
                              ContainerOut& container_out,
                              std::size_t goal_size, RNG& rng) 
{
    using IndexDist = std::uniform_int_distribution<std::size_t>;
    using Value = typename ContainerIn::value_type;
    if (container_in.empty()) {
        return;
    }
    if (container_in.size() <= goal_size) {
        container_out.insert(container_out.end(), container_in.begin(),
                             container_in.end());
        return;
    }
    std::size_t expected_max_out(1.1 * goal_size + container_out.size());
    container_out.reserve(expected_max_out);
    double probability = double(goal_size) / container_in.size();
    if (probability > 0.1) {
        IndexDist cdist(0, container_in.size() - 1);
        for (const Value& v : container_in) {
            std::size_t x = cdist(rng);
            if (x < goal_size) {
                container_out.push_back(v);
            }
        }
    } else {
        std::geometric_distribution<std::size_t> tdist(probability);
        auto i = container_in.begin(), e = container_in.end();
        i += tdist(rng);
        while (i < e) {
            container_out.push_back(*i);
            i += tdist(rng) + 1;
        }
    }
}

template <typename Container, typename RNG>
inline std::vector<typename Container::value_type>
sample_from_range(const Container& container, std::size_t goal_size, RNG& rng) {
    using Value = typename Container::value_type;
    std::vector<Value> result;
    sample_from_range(container, result, goal_size, rng);
    return result;
}

} // namespace sammy

#endif
==> ./implied_vertices.h <==
#ifndef SAMMY_IMPLIED_VERTICES_H_INCLUDED_
#define SAMMY_IMPLIED_VERTICES_H_INCLUDED_

#include "algorithm_ex.h"
#include "dynamic_bitset.h"
#include "literals.h"
#include "pair_infeasibility_map.h"
#include "range.h"
#include "shared_db_propagator.h"
#include "vertex_operations.h"
#include "output.h"

namespace sammy {

/**
 * On demand elimination of implied vertices.
 * Returns a list of vertices such that covering all of them guarantees a full
 * cover, by removing all vertices that are implied by another vertex.
 */
inline std::vector<Vertex>
eliminate_implied_vertices(const std::vector<Vertex>& vertices,
                           SharedDBPropagator& propagator) {
    if (propagator.get_current_level() > 0)
        throw std::logic_error(
            "eliminate_implied_vertices: decision level > 0");
    DynamicBitset implied(vertices.size(), false);
    std::vector<std::vector<std::pair<Lit, std::size_t>>> partners_of(
        2 * propagator.db().num_vars());
    for (std::size_t i = 0, n = vertices.size(); i < n; ++i) {
        Vertex v = vertices[i];
        if (propagator.is_true(v.first) && propagator.is_true(v.second)) {
            implied[i].set();
            continue;
        }
        partners_of[v.first].emplace_back(v.second, i);
        partners_of[v.second].emplace_back(v.first, i);
    }
    for (std::size_t i = 0, n = vertices.size(); i < n; ++i) {
        if (implied[i])
            continue;
        Vertex v = vertices[i];
        reset_and_push_noresolve(propagator, v);
        for (Lit l : propagator.get_trail()) {
            auto& ps = partners_of[l];
            ps.erase(std::remove_if(ps.begin(), ps.end(),
                                    [&](std::pair<Lit, std::size_t> entry) {
                                        if (implied[entry.second])
                                            return true;
                                        if (entry.second == i)
                                            return false;
                                        if (propagator.is_true(entry.first)) {
                                            implied[entry.second].set();
                                            return true;
                                        }
                                        return false;
                                    }),
                     ps.end());
        }
    }
    std::vector<Vertex> result;
    for (std::size_t i = 0, n = vertices.size(); i < n; ++i) {
        if (!implied[i])
            result.push_back(vertices[i]);
    }
    propagator.reset_or_throw();
    return result;
}

inline std::vector<Vertex>
eliminate_implied_vertices(const std::vector<Vertex>& vertices,
                           const SharedDBPropagator& propagator) {
    SharedDBPropagator propagator_(propagator);
    return eliminate_implied_vertices(vertices, propagator_);
}

/**
 * Lookup class for implied vertices.
 * Serves as storage for the entire
 * (reduced) universe of vertices.
 * Offers four main operations:
 *  - Eliminate all implied vertices,
 *    i.e., vertices that are implicitly covered
 *    if another vertex is covered.
 *  - Check if a single vertex is implied, and if so,
 *    by which non-implied vertex.
 *  - Remove all implied vertices from a given list.
 *    Usually, when this is called, the list
 *    to be reduced will contain the implying vertex as well.
 *    Used primarily to reduce the set of uncovered vertices
 *    in a predictable and consistent manner.
 *  - Replace all implied vertices by their implier.
 *    Usually, when this is called, the list
 *    does not contain both a vertex and its implier.
 *    Used primarily to replace vertices in a mutually exclusive set/clique
 *    by canonical implying vertices.
 * Offers a 'partial'/'limited' elimination version that may
 * let some implied vertices slip through, trading off vs. elimination
 * runtime.
 */
class ImpliedVertexCache {
  public:
    explicit ImpliedVertexCache(const PairInfeasibilityMap* inf_map,
                                std::size_t universe_size) noexcept
        : m_pair_inf_map(inf_map), m_universe_size(universe_size) {}

    /**
     * Compute and cache a reduced version of the universe.
     */
    inline void reduce_universe(ClauseDB& clause_db);

    /**
     * Time-limited version of universe reduction.
     */
    inline void limited_reduce_universe(ClauseDB& clause_db, double time_limit);

    /**
     * Check whether we have computed and stored a reduced universe.
     */
    bool have_reduced_universe() const noexcept {
        return !m_reduced_universe.empty();
    }

    /**
     * Return the size of the reduced universe.
     */
    std::size_t reduced_universe_size() const noexcept {
        return m_reduced_universe.size();
    }

    /**
     * Original universe size.
     */
    std::size_t original_universe_size() const noexcept {
        return m_universe_size;
    }

    /**
     * Get a reference to the reduced universe.
     * Returns an empty vector if we have no reduced universe.
     */
    const std::vector<Vertex>& get_reduced_universe() const noexcept {
        return m_reduced_universe;
    }

    /**
     * Check if the given vertex is implied by another vertex.
     */
    bool is_implied(Vertex v) const noexcept { return m_implied_by.count(v); }

    /**
     * Get the vertex implying the given vertex (or an exception).
     */
    Vertex implying_vertex(Vertex v) const noexcept {
        return m_implied_by.at(v);
    }

    /**
     * Get the vertex implying the given vertex,
     * or the given vertex if the vertex is not implied.
     */
    Vertex implying_or_self(Vertex v) const noexcept {
        auto it = m_implied_by.find(v);
        if (it == m_implied_by.end())
            return v;
        return it->second;
    }

    /**
     * Get a copy of the reduced universe, or the entire universe,
     * if we have no reduced universe.
     */
    std::vector<Vertex> get_universe() const {
        if (have_reduced_universe())
            return m_reduced_universe;
        return m_pair_inf_map->collect_vertices(m_universe_size);
    }

    /**
     * Remove the implied vertices from the given list
     * (throws if universe reduction has not been performed).
     */
    void remove_implied(std::vector<Vertex>& vertices) const {
        if (!have_reduced_universe())
            throw std::logic_error("Universe reduction not performed!");
        auto new_end = std::remove_if(vertices.begin(), vertices.end(),
                                      [&](Vertex v) { return is_implied(v); });
        vertices.erase(new_end, vertices.end());
    }

    /**
     * Remove the implied vertices from the given list (either through the
     * ad-hoc method, or using a cached reduced universe).
     */
    void remove_implied(std::vector<Vertex>& vertices,
                        SharedDBPropagator& propagator) const {
        if (!have_reduced_universe()) {
            vertices = eliminate_implied_vertices(vertices, propagator);
        } else {
            remove_implied(vertices);
        }
    }

    /**
     * Remove the implied vertices from the given list (const propagator
     * version).
     */
    void remove_implied(std::vector<Vertex>& vertices,
                        const SharedDBPropagator& propagator) const {
        if (!have_reduced_universe()) {
            SharedDBPropagator copy(propagator);
            remove_implied(vertices, copy);
        } else {
            remove_implied(vertices);
        }
    }

    /**
     * Replace all implied vertices in the given list
     * by their implier; does nothing if we do not actually
     * have a reduced universe.
     */
    void replace_implied(std::vector<Vertex>& vertices) const {
        if (!have_reduced_universe())
            return;
        std::transform(vertices.begin(), vertices.end(), vertices.begin(),
                       [&](Vertex v) { return implying_or_self(v); });
    }

  private:
    /**
     * Member struct that implements the actual elimination algorithm.
     */
    struct EliminationAlgorithm {
        explicit EliminationAlgorithm(ClauseDB& clauses,
                                      ImpliedVertexCache* that)
            : that(that), propagator(&clauses),
              universe(that->m_pair_inf_map->collect_vertices(
                  that->m_universe_size)),
              implier_of(universe.size(),
                         std::numeric_limits<std::size_t>::max()) {}

        bool is_implied(std::size_t vi) const noexcept {
            return implier_of[vi] != std::numeric_limits<std::size_t>::max();
        }

        inline void
        compute_impliers_handle_trail_literal(Lit trail_literal,
                                              std::size_t pushed_index);
        inline std::size_t mark_and_merge_two_sorted_lists(Lit implier, 
                                                           Lit implied);
        inline void compute_literal_partners_of();
        inline bool handle_level_zero();
        inline bool handle_level_zero_both();
        inline void compute_impliers();
        inline void ensure_pushed(Lit &prev_first, Vertex vertex);
        inline void limited_compute_impliers_const_per_literal(
            std::size_t num_per_literal = 20);
        inline void limited_compute_impliers(double time_limit);
        inline void compute_impliers_single_literal();
        inline void compress_path(std::size_t v);

        /**
         * Compress all paths of length > 1 by
         * replacing each intermediate implier
         * by the corresponding ultimate implier,
         * i.e., non-implied vertex.
         */
        void compress_paths() {
            for (std::size_t v : range(universe.size()))
                compress_path(v);
        }

        inline void export_to_cache();

        /**
         * A reference to the ImpliedVertexCache we are doing computation for.
         */
        ImpliedVertexCache* that;

        /**
         * A propagator that we use to detect implications between vertices.
         */
        SharedDBPropagator propagator;

        /**
         * The full, unreduced universe.
         */
        std::vector<Vertex> universe;

        /**
         * For each literal, the list of vertices that contain that literal
         * and the other literal that makes up the vertex.
         */
        std::vector<std::vector<std::pair<Lit, std::size_t>>> partners_of;

        /**
         * For each vertex (indices into universe), the vertex that implies it,
         * or std::numeric_limits<std::size_t>::max() if there is no implier.
         */
        std::vector<std::size_t> implier_of;

        /**
         * Cache of pointers we need to update during path compression.
         */
        std::vector<std::size_t*> pcompress_cache;
    };

    std::vector<Vertex> m_reduced_universe;
    VertexMapTo<Vertex> m_implied_by;
    const PairInfeasibilityMap* m_pair_inf_map;
    std::size_t m_universe_size;
};

void ImpliedVertexCache::EliminationAlgorithm::compute_literal_partners_of() {
    partners_of.resize(2 * propagator.db().num_vars());
    for (std::size_t i = 0, n = universe.size(); i < n; ++i) {
        Vertex v = universe[i];
        if (propagator.is_true(v.first) && propagator.is_true(v.second)) {
            implier_of[i] = i;
            continue;
        }
        partners_of[v.first].emplace_back(v.second, i);
        partners_of[v.second].emplace_back(v.first, i);
    }
    auto sort_partners_list = [](auto& r) {
        const auto compare_1st = [](const auto& e1, const auto& e2) {
            return e1.first < e2.first;
        };
        std::sort(r.begin(), r.end(), compare_1st);
    };
    std::for_each(partners_of.begin(), partners_of.end(), sort_partners_list);
}

std::size_t ImpliedVertexCache::EliminationAlgorithm::
    mark_and_merge_two_sorted_lists(Lit implier, Lit implied) 
{
    auto& left_list = partners_of[implier];
    auto& right_list = partners_of[implied];
    std::size_t implier_implied_index = std::numeric_limits<std::size_t>::max();
    if (left_list.empty() || right_list.empty())
        return implier_implied_index;

    auto left_in = left_list.begin(), left_out = left_list.begin(),
         left_end = left_list.end();
    auto right_in = right_list.begin(), right_out = right_list.begin(),
         right_end = right_list.end();

    auto left_skip_implied = [&] () {
        while (is_implied(left_in->second)) {
            if (++left_in == left_end)
                return false;
        }
        if(left_in->first == implied) {
            implier_implied_index = left_in->second;
        }
        return true;
    };

    auto right_skip_implied = [&] () {
        while (is_implied(right_in->second)) {
            if (++right_in == right_end)
                return false;
        }
        return true;
    };

    auto advance_left = [&] () {
        *left_out++ = *left_in;
        if (++left_in == left_end) {
            return false;
        }
        return left_skip_implied();
    };

    auto advance_right = [&] () {
        *right_out++ = *right_in;
        if (++right_in == right_end) {
            return false;
        }
        return right_skip_implied();
    };

    auto advance_both = [&] () {
        *left_out++ = *left_in;
        ++left_in, ++right_in;
        if (left_in != left_end) {
            left_skip_implied();
        }
        if (right_in != right_end) {
            right_skip_implied();
        }
        return left_in != left_end && right_in != right_end;
    };

    auto handle_left_remainder = [&]() {
        if (left_in == left_end)
            return;
        while (advance_left())
            ;
    };

    auto handle_right_remainder = [&]() {
        if (right_in == right_end)
            return;
        while (advance_right())
            ;
    };

    if (!left_skip_implied()) {
        left_list.clear();
        if (right_skip_implied()) {
            handle_right_remainder();
        }
        right_list.erase(right_out, right_end);
        return implier_implied_index;
    }

    if (!right_skip_implied()) {
        right_list.clear();
        handle_left_remainder();
        left_list.erase(left_out, left_end);
        return implier_implied_index;
    }

    for (;;) {
        Lit left = left_in->first, right = right_in->first;
        if (left < right) {
            if (!advance_left())
                break;
        } else if (right < left) {
            if (!advance_right())
                break;
        } else {
            assert(left_in->second != right_in->second);
            implier_of[right_in->second] = left_in->second;
            if (!advance_both())
                break;
        }
    }
    handle_right_remainder();
    handle_left_remainder();
    left_list.erase(left_out, left_end);
    right_list.erase(right_out, right_end);
    return implier_implied_index;
}

void ImpliedVertexCache::EliminationAlgorithm::
    compute_impliers_single_literal() {
    for (Lit l = 0, nall = 2 * propagator.db().num_vars(); l < nall; ++l) {
        // find any non-implied interaction involving l;
        // avoid doing propagation if there are no un-implied
        // interactions involving l
        auto& p_of_l = partners_of[l];
        while (!p_of_l.empty() && is_implied(p_of_l.back().second)) {
            p_of_l.pop_back();
        }
        if (p_of_l.empty() || propagator.is_true(l))
            continue;
        if (!propagator.push_level(l)) {
            throw std::logic_error("Infeasible interaction in universe!");
        }
        for (Lit l2 : propagator.get_trail()) {
            if (l2 == l)
                continue;
            // if l implies l2, (l2, x) is implied by (l, x) for any x.
            // so we walk the sorted partner lists of l and l2 and
            // mark (l2,x) for x for which (l,x) is present,
            // stripping both lists of implied elements as we go.
            std::size_t l_l2_index = mark_and_merge_two_sorted_lists(l, l2);
            if(l_l2_index != std::numeric_limits<std::size_t>::max() &&
               !is_implied(l_l2_index)) 
            {
                // also, (l, l2) is implied by any other (l, y).
                if(p_of_l.size() > 1) {
                    auto e_ly = p_of_l[0];
                    if(e_ly.first == l2) e_ly = p_of_l[1];
                    std::size_t implier_index = e_ly.second;
                    while(is_implied(implier_index)) {
                        implier_index = implier_of[implier_index];
                    }
                    if(implier_index != l_l2_index) {
                        implier_of[l_l2_index] = implier_index;
                    }
                }
            }
        }
        propagator.pop_level();
    }
}

void ImpliedVertexCache::EliminationAlgorithm::
    compute_impliers_handle_trail_literal(Lit trail_literal,
                                          std::size_t pushed_index) {
    auto& ps = partners_of[trail_literal];
    ps.erase(std::remove_if(ps.begin(), ps.end(),
                            [&](std::pair<Lit, std::size_t> entry) {
                                // simply drop already-implied vertices
                                if (is_implied(entry.second))
                                    return true;
                                // avoid detecting vertices as implying
                                // themselves
                                if (entry.second == pushed_index)
                                    return false;
                                // if the other literal is true, pushed_index
                                // implies entry.first
                                if (propagator.is_true(entry.first)) {
                                    implier_of[entry.second] = pushed_index;
                                    return true;
                                }
                                // not known to be implied
                                return false;
                            }),
             ps.end());
}

bool ImpliedVertexCache::EliminationAlgorithm::handle_level_zero_both() {
    propagator.reset_or_throw();
    const auto& trail = propagator.get_trail();
    Lit nconc = that->m_pair_inf_map->get_n_concrete();
    std::vector<std::size_t> implied_pairs;
    for(Lit l1 : trail) {
        if(l1 >= nconc) continue;
        auto is_implied_pair = [&] (std::pair<Lit, std::size_t> entry) {
            if(propagator.is_true(entry.first)) {
                if(entry.first > l1) {
                    implied_pairs.push_back(entry.second);
                }
                return true;
            }
            return false;
        };
        auto& p_of_l = partners_of[l1];
        auto new_end = std::remove_if(p_of_l.begin(), p_of_l.end(), 
                                      is_implied_pair);
        p_of_l.erase(new_end, p_of_l.end());
    }
    if(implied_pairs.size() == universe.size()) {
        // all pairs are implied by any configuration;
        // keep only one pair; also allows us to skip the rest
        for(std::size_t i = 1, s = universe.size(); i < s; ++i) {
            implier_of[i] = 0;
        }
        return true;
    }
    if(implied_pairs.empty()) {
        return false;
    }
    std::sort(implied_pairs.begin(), implied_pairs.end());
    std::size_t last = 0;
    --last;
    for(std::size_t implied : implied_pairs) {
        if(implied != ++last) {
            break;
        }
    }
    for(std::size_t implied : implied_pairs) {
        implier_of[implied] = last;
    }
    return false;
}

bool ImpliedVertexCache::EliminationAlgorithm::handle_level_zero() {
    // first: handle pairs both in level 0 (i.e., in ANY configuration)
    if(handle_level_zero_both()) {
        // if that is all pairs, every interaction except one is
        // marked as implied and we can skip the rest
        return true;
    }

    // then handle individual literals true in every configuration
    const auto& trail = propagator.get_trail();
    Lit nconc = 2 * that->m_pair_inf_map->get_n_concrete();
    auto is_concrete = [nconc] (Lit l) { return l < nconc; };
    if(!std::any_of(trail.begin(), trail.end(), is_concrete)) {
        // no concrete literals in the trail
        return false;
    }
    std::vector<std::pair<Lit, std::size_t>> nontrue_partner_of(nconc, {NIL, 0});
    for(Lit l = 0; l < nconc; ++l) {
        auto& p_of_l = partners_of[l];
        for(std::pair<Lit, std::size_t> entry : p_of_l) {
            if(!propagator.is_true(entry.first)) {
                nontrue_partner_of[l] = entry;
                break;
            }
        }
    }
    for(Lit l : trail) {
        if(l >= nconc) continue;
        auto is_implied_pair = [&] (std::pair<Lit, std::size_t> entry) {
            if(is_implied(entry.second)) {
                return true;
            }
            auto nontrue_partner = nontrue_partner_of[entry.first];
            if(nontrue_partner.first == NIL) {
                return false;
            }
            implier_of[entry.second] = nontrue_partner.second;
            return true;
        };
        auto& p_of_l = partners_of[l];
        auto new_end = std::remove_if(p_of_l.begin(), p_of_l.end(),
                                      is_implied_pair);
        p_of_l.erase(new_end, p_of_l.end());
    }
    return false;
}

void ImpliedVertexCache::EliminationAlgorithm::
    ensure_pushed(Lit &previous_first, Vertex v)
{
    // it suffices to check the level 2 part of the trail,
    // since we handle implied-by-single-literal cases earlier.
    // if the previous vertex had the same first literal,
    // we can skip the reset and thus usually the first decision level
    if(previous_first != v.first) {
        previous_first = v.first;
        reset_and_push_noresolve(propagator, v);
    } else {
        auto push_2nd_if_necessary = [&] () {
            if(!propagator.is_true(v.second) && 
                !propagator.push_level(v.second)) 
            {
                throw std::logic_error("Infeasible interaction in universe!");
            }
        };

        if(propagator.get_current_level() == 2) {
            // in this case, we have v.first and the first old
            // literal pushed; push v.second
            propagator.pop_level();
            push_2nd_if_necessary();
        } else if(propagator.get_current_level() == 1) {
            // in this case, we have either v.first pushed or
            // v.second pushed (if v.first is fixed at level 0)
            if(*propagator.current_level_begin() != v.first) {
                // v.second is pushed
                propagator.pop_level();
            }
            push_2nd_if_necessary();
        } else {
            throw std::logic_error("Implausible decision level!");
        }
    }
}

void ImpliedVertexCache::EliminationAlgorithm::compute_impliers() {
    Lit previous_first = NIL;
    for (std::size_t i : range(universe.size())) {
        if (is_implied(i))
            continue;
        Vertex v = universe[i];
        ensure_pushed(previous_first, v);
        if(propagator.get_current_level() < 2)
            continue;
        auto deepest_level = IteratorRange{
            propagator.current_level_begin(),
            propagator.get_trail().end()
        };
        for (Lit l : deepest_level) {
            compute_impliers_handle_trail_literal(l, i);
        }
    }
}

void ImpliedVertexCache::EliminationAlgorithm::
    limited_compute_impliers_const_per_literal(std::size_t num_per_literal)
{
    const Lit nconc = 2 * that->m_pair_inf_map->get_n_concrete();
    std::vector<std::pair<Lit, std::size_t>> random_partners;
    auto& rng = sammy::rng();
    for(Lit l = 0; l < nconc; ++l) {
        auto& p_of_l = partners_of[l];
        auto new_end = std::remove_if(p_of_l.begin(), p_of_l.end(),
                                      [&](std::pair<Lit, std::size_t> entry) {
                                          return is_implied(entry.second);
                                      });
        p_of_l.erase(new_end, p_of_l.end());
        random_partners.clear();
        sample_from_range(p_of_l, random_partners, num_per_literal, rng);
        if(random_partners.empty()) {
            continue;
        }
        propagator.reset_or_throw();
        if(!propagator.is_open(l)) {
            continue;
        }
        if(!propagator.push_level(l)) {
            throw std::logic_error("Infeasible interaction in universe!");
        }
        for(auto&& e : random_partners) {
            if(is_implied(e.second) || !propagator.is_open(e.first)) {
                continue;
            }
            if(!propagator.push_level(e.first)) {
                throw std::logic_error("Infeasible interaction in universe!");
            }
            IteratorRange deepest_level{
                propagator.current_level_begin(),
                propagator.get_trail().end()
            };
            for(Lit l2 : deepest_level) {
                compute_impliers_handle_trail_literal(l2, e.second);
            }
            propagator.pop_level();
        }
    }
    propagator.reset_or_throw();
}

void ImpliedVertexCache::EliminationAlgorithm::
    limited_compute_impliers(double time_limit) 
{
    limited_compute_impliers_const_per_literal();
    auto begin_time = std::chrono::steady_clock::now();
    std::size_t check_count = 0;
    Lit prev_first = NIL;
    for (std::size_t vertex_index = 0, us = universe.size();
         vertex_index < us; ++vertex_index) 
    {
        if (is_implied(vertex_index))
            continue;
        Vertex vertex = universe[vertex_index];
        ensure_pushed(prev_first, vertex);
        if (++check_count == 16384) {
            check_count = 0;
            if (seconds_between(begin_time, std::chrono::steady_clock::now()) >=
                time_limit)
            {
                break;
            }
        }
        if(propagator.get_current_level() < 2)
            continue;
        assert(propagator.get_current_level() == 2);
        IteratorRange deepest_level{
            propagator.current_level_begin(),
            propagator.get_trail().end()
        };
        for (Lit l : deepest_level) {
            compute_impliers_handle_trail_literal(l, vertex_index);
        }
    }
    propagator.reset_or_throw();
}

void ImpliedVertexCache::EliminationAlgorithm::compress_path(std::size_t vi) {
    // check non-impliedness/self-impliedness
    std::size_t* current = &implier_of[vi];
    std::size_t cval = *current;
    if (cval == std::numeric_limits<std::size_t>::max() || cval == vi)
        return;
    // get ready to compress the path
    pcompress_cache.clear();
    // find entries we need to update
    for (;;) {
        std::size_t* next = &implier_of[cval];
        std::size_t nval = *next;
        if (nval == std::numeric_limits<std::size_t>::max())
            break;
        pcompress_cache.push_back(current);
        cval = nval;
        current = next;
    }
    // update all entries to the ultimate implier
    for (std::size_t* ptr : pcompress_cache) {
        *ptr = cval;
    }
}

void ImpliedVertexCache::EliminationAlgorithm::export_to_cache() {
    auto& implied_by = that->m_implied_by;
    auto& unimplied = that->m_reduced_universe;
    for (std::size_t vi : range(universe.size())) {
        Vertex v = universe[vi];
        if (implier_of[vi] != std::numeric_limits<std::size_t>::max()) {
            Vertex implier = universe[implier_of[vi]];
            implied_by.insert_or_assign(v, implier);
        } else {
            unimplied.push_back(v);
        }
    }
}

void ImpliedVertexCache::reduce_universe(ClauseDB& clause_db) {
    EliminationAlgorithm algorithm(clause_db, this);
    algorithm.compute_literal_partners_of();
    if(!algorithm.handle_level_zero()) {
        algorithm.compute_impliers_single_literal();
        algorithm.compute_impliers();
        algorithm.compress_paths();
    }
    algorithm.export_to_cache();
}

void ImpliedVertexCache::limited_reduce_universe(ClauseDB& clause_db,
                                                 double time_limit) {
    if (time_limit >= 0.0 && !std::isfinite(time_limit)) {
        reduce_universe(clause_db);
        return;
    }
    auto begin_time = std::chrono::steady_clock::now();
    EliminationAlgorithm algorithm(clause_db, this);
    algorithm.compute_literal_partners_of();
    if(!algorithm.handle_level_zero()) {
        algorithm.compute_impliers_single_literal();
        double trem = time_limit - 
            seconds_between(begin_time, std::chrono::steady_clock::now());
        if (trem > 0.0) {
            algorithm.limited_compute_impliers(trem);
        }
        algorithm.compress_paths();
    }
    algorithm.export_to_cache();
}

} // namespace sammy

#endif
==> ./literals.h <==
#ifndef SAMMY_LITERALS_H_INCLUDED_
#define SAMMY_LITERALS_H_INCLUDED_

#include "time.h"
#include <algorithm>
#include <boost/container/small_vector.hpp>
#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/unordered/unordered_flat_set.hpp>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace sammy {

/**
 * A literal in the external sense:
 *  - Positive literals are positive numbers,
 *  - Negative literals are negative numbers,
 *  - 0 is not a valid literal.
 */
using ExternalLit = std::int32_t;

/**
 * An external clause, represented as vector of literals.
 */
using ExternalClause = std::vector<ExternalLit>;

/**
 * A literal in the internal sense:
 *  - Even numbers are positive literals,
 *  - Odd numbers are negative literals.
 *  - For variables x, we thus have 2 * x and 2 * x + 1 as literals.
 */
using Lit = std::uint32_t;

/**
 * An internal clause represented as vector of internal literals.
 */
using CVec = std::vector<Lit>;

/**
 * An internal clause represented as small vector of internal literals.
 */
using SCVec = boost::container::small_vector<Lit, 4>;

/**
 * A variable in the internal sense (technically equal type to Lit).
 * Only used to make types clearer in the interfaces.
 */
using Var = Lit;

/**
 * A clause reference to a clause in the database.
 */
using CRef = Lit;

/**
 * A vertex in the universe of pairwise interactions.
 */
using Vertex = std::pair<Lit, Lit>;

/**
 * A vertex in the universe of pairwise interactions
 * in external representation.
 */
using ExternalVertex = std::pair<ExternalLit, ExternalLit>;

/**
 * A value that indicates an invalid variable/literal/clause.
 */
static constexpr Lit NIL = std::numeric_limits<Lit>::max();

namespace lit {

/**
 * @brief External to internal literal conversion.
 *
 * @param l
 * @return Lit
 */
static inline constexpr Lit internalize(ExternalLit l) noexcept {
    return Lit(2) * static_cast<Lit>((std::abs)(l)-1) + Lit(l < 0);
}

static inline constexpr Vertex internalize(ExternalVertex v) noexcept {
    return {internalize(v.first), internalize(v.second)};
}

static inline std::vector<Vertex>
internalize(const std::vector<ExternalVertex>& vertices) {
    std::vector<Vertex> result;
    result.reserve(vertices.size());
    std::transform(vertices.begin(), vertices.end(), std::back_inserter(result),
                   [](ExternalVertex v) { return internalize(v); });
    return result;
}

/**
 * @brief Internal to external literal conversion.
 *
 * @param l
 * @return constexpr ExternalLit
 */
static inline constexpr ExternalLit externalize(Lit l) noexcept {
    ExternalLit result = static_cast<ExternalLit>(l >> 1) + 1;
    if (l & 1)
        result = -result;
    return result;
}

static inline std::vector<std::pair<ExternalLit, ExternalLit>>
externalize(const std::vector<Vertex>& vertices) {
    std::vector<std::pair<ExternalLit, ExternalLit>> result;
    result.reserve(vertices.size());
    std::transform(vertices.begin(), vertices.end(), std::back_inserter(result),
                   [](const auto& v) {
                       return std::pair<ExternalLit, ExternalLit>{
                           externalize(v.first), externalize(v.second)};
                   });
    return result;
}

/**
 * @brief Internal to external literal conversion.
 */
template <
    typename ClauseContainer,
    std::enable_if_t<std::is_integral_v<typename ClauseContainer::value_type>,
                     int> = 0>
static inline ExternalClause externalize(const ClauseContainer& clause) {
    ExternalClause result;
    std::transform(clause.begin(), clause.end(), std::back_inserter(result),
                   [](Lit l) -> ExternalLit { return externalize(l); });
    return result;
}

/**
 * @brief Internal literal negation.
 *
 * @param l
 * @return constexpr Lit
 */
static inline constexpr Lit negate(Lit l) noexcept { return l ^ Lit(1); }

/**
 * @brief Extract the variable from a literal.
 *
 * @param l
 * @return constexpr Lit
 */
static inline constexpr Var var(Lit l) noexcept { return l >> 1; }

/**
 * @brief Turn a variable into its positive literal.
 */
static inline constexpr Lit positive_lit(Var v) noexcept { return v << 1; }

/**
 * @brief Turn a variable into its negative literal.
 */
static inline constexpr Lit negative_lit(Var v) noexcept {
    return (v << 1) + 1;
}

/**
 * @brief Check for negative literal.
 *
 * @param l
 * @return true if the literal is negative.
 * @return false if the literal is positive.
 */
static inline constexpr bool negative(Lit l) noexcept { return l & Lit(1); }

template <typename BitsetType>
static inline bool is_true_in(Lit l, const BitsetType& assignment) noexcept {
    bool value(assignment[var(l)]);
    return negative(l) ? !value : value;
}

template <typename BitsetType>
static inline bool is_false_in(Lit l, const BitsetType& assignment) noexcept {
    return !is_true_in(l, assignment);
}

} // namespace lit

namespace simplify {

static constexpr Lit fixed_positive() noexcept { return NIL - 1; }

static constexpr Lit fixed_negative() noexcept { return NIL; }

static constexpr Lit eliminated() noexcept { return NIL - 2; }

} // namespace simplify

struct PairHash {
    std::size_t operator()(std::pair<int, int> v) const noexcept {
        std::size_t x1 = static_cast<unsigned>(v.first);
        std::size_t x2 = static_cast<unsigned>(v.second);
        x1 = ((x1 << ((CHAR_BIT / 2) * sizeof(std::size_t))) |
              (x1 >> ((CHAR_BIT / 2) * sizeof(std::size_t)))) *
             7;
        x1 += 17 * x2;
        return x1;
    }

    std::size_t
    operator()(std::pair<std::uint32_t, std::uint32_t> v) const noexcept {
        std::size_t x1 = v.first;
        std::size_t x2 = v.second;
        x1 = ((x1 << ((CHAR_BIT / 2) * sizeof(std::size_t))) |
              (x1 >> ((CHAR_BIT / 2) * sizeof(std::size_t)))) *
             7;
        x1 += 17 * x2;
        return x1;
    }

    std::size_t
    operator()(std::pair<std::size_t, std::size_t> v) const noexcept {
        std::size_t x = 7 * v.first;
        x = (x << 25) | (x >> (sizeof(x) * CHAR_BIT - 25));
        std::size_t y = 17 * v.second;
        return x + y;
    }

    using value = std::size_t;
};

template <typename Element> using HashSet = boost::unordered_flat_set<Element>;

template <typename Key, typename Value>
using HashMap = boost::unordered_flat_map<Key, Value>;

using EdgeSet = boost::unordered_flat_set<std::pair<Lit, Lit>>;

template <typename Value>
using VertexMapTo = boost::unordered_flat_map<std::pair<Lit, Lit>, Value>;

template <typename Key> using PairHashSet = boost::unordered_flat_set<Key>;

template <typename Key, typename Value>
using PairMapTo = boost::unordered_flat_map<Key, Value>;

} // namespace sammy

#endif
==> ./detect_equalities.h <==
#ifndef SAMMY_DETECT_EQUALITIES_H_INCLUDED_
#define SAMMY_DETECT_EQUALITIES_H_INCLUDED_

#include "equality_graph.h"
#include "range.h"
#include "shared_db_propagator.h"
#include "simplification_stats.h"
#include "simplify_datastructure.h"
#include "stamp_set.h"
#include <cassert>

namespace sammy {

namespace detail {

static inline auto implied_literals(const SharedDBPropagator& prop) noexcept {
    return IteratorRange{prop.current_level_begin() + 1,
                         prop.get_trail().end()};
}

static inline bool fix_unaries(const SharedDBPropagator& propagator,
                               EqualityGraph& equality_graph,
                               SimplificationStats* stats) {
    bool last_changed = false;
    for (Lit l : propagator.get_trail()) {
        if (equality_graph.make_true(l)) {
            stats->variables_fixed += 1;
            last_changed = true;
        }
    }
    return last_changed;
}

} // namespace detail

class ImplicationGraphBuilder {
  public:
    explicit ImplicationGraphBuilder(SimplifyDatastructure* simplify_ds)
        : simplify_ds(simplify_ds), clauses(simplify_ds->extract_clause_db()),
          old_counts(clauses.get_clause_counts()), propagator(&clauses),
          equality_graph(simplify_ds->original_num_vars()),
          implications_of(2 * simplify_ds->original_num_vars()),
          additional_implications_of(2 * simplify_ds->original_num_vars()),
          stats(&stats_buffer),
          reachable(2 * simplify_ds->original_num_vars()) {}

    void set_stats(SimplificationStats* stats) noexcept { this->stats = stats; }

    bool run_fixing_and_equality_iteration() {
        bool changed = false;
        bool last_changed;
        do {
            last_changed = p_init_iteration();
            for (Var v = 0, nv = simplify_ds->original_num_vars(); v < nv; ++v)
            {
                if (simplify_ds->is_eliminated(v) ||
                    !propagator.is_open(lit::positive_lit(v)))
                {
                    continue;
                }
                last_changed |= p_push_var(v);
            }
            p_add_additional_arcs();
            last_changed |= p_tarjan_scc_find_equalities();
            changed |= last_changed;
        } while (last_changed);
        if (changed) {
            simplify_ds->add_clauses(clauses, old_counts);
            simplify_ds->apply_fixes_and_equalities(
                equality_graph.compute_old_to_new_map());
        }
        return changed;
    }

    bool run_contrapositive_and_strengthening_iteration() {
        return p_add_missing_contrapositivities() ||
               p_strengthen_with_implications();
    }

  private:
    bool p_strengthen_with_implications() {
        auto clauses_with_lit = p_make_clauses_with_lit();
        auto& clauses = simplify_ds->clauses();
        auto is_implied = [&](Lit l) { return reachable.contains(l); };
        bool result = false;
        for (Lit l = 0, nl = 2 * simplify_ds->original_num_vars(); l < nl; ++l)
        {
            if (simplify_ds->is_eliminated(lit::var(l)))
                continue;
            if (clauses_with_lit[l].empty())
                continue;
            p_stamp_reachable_black(l, true);
            reachable.erase(l);
            for (CRef cr : clauses_with_lit[l]) {
                auto& clause = clauses[cr];
                if (std::any_of(clause.begin(), clause.end(), is_implied)) {
                    stats->clauses_strengthened_by_binary_resolution += 1;
                    auto iter = std::find(clause.begin(), clause.end(), l);
                    std::swap(*iter, clause.back());
                    clause.pop_back();
                    if (clause.size() <= 2) {
                        result = true;
                    }
                }
            }
        }
        return result;
    }

    void p_stamp_reachable_black(Lit v, bool clear_first) {
        tj_stack.clear();
        if (clear_first)
            reachable.clear();
        reachable.insert(v);
        tj_stack.push_back(v);
        while (!tj_stack.empty()) {
            Lit w = tj_stack.back();
            tj_stack.pop_back();
            for (Lit x : implications_of[w]) {
                if (reachable.check_insert(x)) {
                    tj_stack.push_back(x);
                }
            }
        }
    }

    bool p_add_missing_contrapositivities() {
        Lit rlast = NIL;
        auto is_last_pushed = [&](Lit w) { return w == rlast; };
        std::size_t added_clauses = 0;
        auto& clauses = simplify_ds->clauses();
        for (Lit v : tj_order) {
            auto& aio = additional_implications_of[v];
            if (aio.empty())
                continue;
            auto& io = implications_of[v];
            p_stamp_reachable_black(
                v, !std::any_of(io.begin(), io.end(), is_last_pushed));
            rlast = v;
            for (Lit w : aio) {
                if (!reachable.contains(w)) {
                    if (arc_set.emplace(v, w).second) {
                        Lit lits[2] = {lit::negate(v), w};
                        io.push_back(w);
                        clauses.emplace_back(+lits, lits + 2);
                        ++added_clauses;
                    }
                }
            }
            aio.clear();
        }
        stats->binary_contrapositivities_learned += added_clauses;
        return added_clauses != 0;
    }

    void p_tarjan_scc_strongconnect(Lit v) {
        auto& data = tj_dfs_data[v];
        data.index = tj_index;
        data.lowlink = tj_index++;
        tj_stack.push_back(v);
        tj_on_stack[v] = true;
        auto handle_neighbor = [&](Lit w) {
            if (tj_dfs_data[w].index == NIL) {
                p_tarjan_scc_strongconnect(w);
                data.lowlink = (std::min)(data.lowlink, tj_dfs_data[w].lowlink);
            } else if (tj_on_stack[w]) {
                data.lowlink = (std::min)(data.lowlink, tj_dfs_data[w].index);
            }
        };
        for (Lit w : implications_of[v]) {
            handle_neighbor(w);
        }
        for (Lit w : additional_implications_of[v]) {
            handle_neighbor(w);
        }
        if (data.lowlink == data.index) {
            // v is root of a strongly connected component
            Lit w;
            while ((w = tj_stack.back()) != v) {
                if (equality_graph.make_equal(v, w)) {
                    tj_changes += 1;
                    stats->variables_eliminated_by_equality += 1;
                }
                tj_stack.pop_back();
                tj_on_stack[w] = false;
            }
            tj_stack.pop_back();
            tj_on_stack[v] = false;
            tj_order.push_back(v);
        }
    }

    bool p_tarjan_scc_find_equalities() {
        std::size_t nl = 2 * simplify_ds->original_num_vars();
        tj_dfs_data.assign(nl, TarjanSCCDFSInfo{});
        tj_on_stack.assign(nl, false);
        tj_stack.clear();
        tj_order.clear();
        tj_index = tj_changes = 0;
        for (Lit v = 0; v < nl; ++v) {
            if (tj_dfs_data[v].index != NIL)
                continue;
            p_tarjan_scc_strongconnect(v);
        }
        return tj_changes != 0;
    }

    bool p_init_iteration() {
        arc_set.clear();
        auto& io = implications_of;
        auto& aio = additional_implications_of;
        auto clear_vector = [](std::vector<Lit>& v) { v.clear(); };
        std::for_each(io.begin(), io.end(), clear_vector);
        std::for_each(aio.begin(), aio.end(), clear_vector);
        for (auto [l1, l2] : clauses.binary_clauses()) {
            if (propagator.is_open(l1) && propagator.is_open(l2)) {
                p_add_double_arc(lit::negate(l1), l2);
            }
        }
        return detail::fix_unaries(propagator, equality_graph, stats);
    }

    bool p_push_var(Var v) {
        Lit p = lit::positive_lit(v), n = lit::negative_lit(v);
        if (!p_push_lit(p))
            return true;
        p_implications_to_reachable();
        p_analyze_trail(p);
        propagator.pop_level();
        if (!p_push_lit(n))
            return true;
        p_analyze_trail(n);
        bool result = p_infer_fixed_from_doubly_implied();
        propagator.pop_level();
        if (result)
            propagator.incorporate_or_throw();
        return result;
    }

    void p_implications_to_reachable() {
        reachable.clear();
        for (Lit l : detail::implied_literals(propagator)) {
            reachable.insert(l);
        }
    }

    bool p_infer_fixed_from_doubly_implied() {
        bool result = false;
        for (Lit l : detail::implied_literals(propagator)) {
            if (reachable.count(l)) {
                result = true;
                clauses.add_clause(&l, &l + 1);
                stats->variables_fixed += 1;
            }
        }
        return result;
    }

    bool p_push_lit(Lit l) {
        if (!propagator.push_level(l)) {
            propagator.resolve_or_throw();
            stats->conflict_clauses_learned += 1;
            stats->variables_fixed += 1;
            stats->failed_literals += 1;
            equality_graph.make_false(l);
            return false;
        }
        return true;
    }

    void p_analyze_trail(Lit pushed) {
        auto iter_lits = propagator.current_level_begin() + 1,
             end = propagator.get_trail().end();
        auto iter_reasons = propagator.current_level_reasons_begin() + 1;
        while (iter_lits < end) {
            Lit limplied = *iter_lits++;
            Reason rimplied = *iter_reasons++;
            assert(rimplied.reason_length >= 2);
            if (rimplied.reason_length == 2)
                continue;
            p_buffer_clause(rimplied.clause);
            if (clause_analysis_buffer.size() == 2) {
                p_add_double_arc(lit::negate(clause_analysis_buffer[0]),
                                 clause_analysis_buffer[1]);
            } else {
                p_add_arc(pushed, limplied);
            }
        }
    }

    void p_add_additional_arcs() {
        for (Lit l = 0, nl = 2 * simplify_ds->original_num_vars(); l < nl; ++l)
        {
            Lit lneg = lit::negate(l);
            for (Lit o : implications_of[l]) {
                Lit oneg = lit::negate(o);
                if (!arc_set.count(std::pair<Lit, Lit>{oneg, lneg})) {
                    additional_implications_of[oneg].push_back(lneg);
                }
            }
        }
    }

    void p_buffer_clause(CRef clause) {
        clause_analysis_buffer.clear();
        auto lits = clauses.lits_of(clause);
        std::copy_if(lits.begin(), lits.end(),
                     std::back_inserter(clause_analysis_buffer), [&](Lit l) {
                         return propagator.get_decision_level(l) != 0;
                     });
    }

    void p_add_arc(Lit from, Lit to) {
        if (arc_set.emplace(from, to).second) {
            implications_of[from].push_back(to);
        }
    }

    void p_add_double_arc(Lit from, Lit to) {
        if (arc_set.emplace(from, to).second) {
            implications_of[from].push_back(to);
        }
        Lit bfrom = lit::negate(to);
        Lit bto = lit::negate(from);
        if (arc_set.emplace(bfrom, bto).second) {
            implications_of[bfrom].push_back(bto);
        }
    }

    std::vector<std::vector<CRef>> p_make_clauses_with_lit() {
        Lit nl = 2 * simplify_ds->original_num_vars();
        std::vector<std::vector<CRef>> clauses_with_lit(nl);
        auto& clauses = simplify_ds->clauses();
        for (std::size_t i = 0, n = clauses.size(); i < n; ++i) {
            const auto& clause = clauses[i];
            if (clause.size() <= 2)
                continue;
            for (Lit l : clause) {
                clauses_with_lit[l].push_back(CRef(i));
            }
        }
        return clauses_with_lit;
    }

    struct TarjanSCCDFSInfo {
        Lit index;
        Lit lowlink;

        TarjanSCCDFSInfo() noexcept {
            index = NIL;
            lowlink = 0;
        }
    };

    SimplifyDatastructure* simplify_ds;
    ClauseDB clauses;
    ClauseCounts old_counts;
    SharedDBPropagator propagator;
    EqualityGraph equality_graph;
    std::vector<std::vector<Lit>> implications_of, additional_implications_of;
    SimplificationStats stats_buffer;
    SimplificationStats* stats;
    EdgeSet arc_set;
    std::vector<Lit> clause_analysis_buffer;
    std::vector<TarjanSCCDFSInfo> tj_dfs_data;
    std::vector<bool> tj_on_stack;
    Lit tj_index, tj_changes;
    std::vector<Lit> tj_stack;
    std::vector<Lit> tj_order;
    StampSet<Lit> reachable;
};

/**
 * A very fast method to detect literals fixed at level 0,
 * and to eliminate those literals from the formula.
 * This method is mainly used for testing & debugging purposes.
 */
inline bool fix_level0_literals(SimplifyDatastructure& simplify_ds,
                                SimplificationStats* stats = nullptr) {
    ClauseDB exported{simplify_ds.extract_clause_db()};
    SharedDBPropagator propagator{&exported};
    EqualityGraph eqgraph{exported.num_vars()};
    SimplificationStats stat_buffer;
    if (!stats)
        stats = &stat_buffer;
    std::size_t count = 0;
    for (Lit l : propagator.get_trail()) {
        eqgraph.make_true(l);
        count += 1;
    }
    stats->variables_fixed += count;
    if (count != 0) {
        simplify_ds.apply_fixes_and_equalities(
            eqgraph.compute_old_to_new_map());
        return true;
    }
    return false;
}

/**
 * Detect fixed, failed and equal literals (by analyzing consequences from
 * pushing and propagating each literal at least once); this identifies most
 * equalities between literals and literals implicitly fixed to true.
 */
inline bool
detect_failed_and_equal_literals(SimplifyDatastructure& simplify_ds,
                                 SimplificationStats* stats = nullptr) {
    ImplicationGraphBuilder detector{&simplify_ds};
    if (stats)
        detector.set_stats(stats);
    if (detector.run_fixing_and_equality_iteration()) {
        return true;
    }
    return detector.run_contrapositive_and_strengthening_iteration();
}

} // namespace sammy

#endif
==> ./clause_reduction.h <==
#ifndef CLAUSE_REDUCTION_H_INCLUDED_
#define CLAUSE_REDUCTION_H_INCLUDED_

#include "literals.h"
#include <algorithm>
#include <optional>
#include <unordered_set>
#include <utility>
#include <vector>

namespace sammy {

/**
 * @brief Reduce an external clause by removing duplicates.
 * Detects tautological clauses (which contain x and -x); returns
 * an empty vector for them. Also sorts the literals.
 *
 * @param x
 * @return std::vector<int>
 */
inline ExternalClause reduce_external_clause(const ExternalClause& x) {
    std::vector<int> result;
    result.reserve(x.size());
    if (x.size() > 256) {
        // use an external set for large clauses.
        std::unordered_set<int> entries;
        for (int l : x) {
            if (entries.count(-l)) {
                // x is a tautology; drop it.
                result.clear();
                return result;
            }
            if (entries.insert(l).second) {
                // filter duplicates.
                result.push_back(l);
            }
        }
    } else {
        // use the result vector for small clauses.
        for (int l : x) {
            bool found = false;
            for (int rl : result) {
                if (rl == -l) {
                    // x is a tautology.
                    result.clear();
                    return result;
                }
                if (rl == l) {
                    found = true;
                }
            }
            if (!found) {
                // filter duplicates.
                result.push_back(l);
            }
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

namespace impl {

static inline void
reduced_remove_identical_pairs(std::vector<ExternalClause>& reduced) {
    auto comp_size = [](const auto& v1, const auto& v2) {
        return v1.size() < v2.size();
    };
    std::sort(reduced.begin(), reduced.end(), comp_size);
    ExternalClause l2(2, 0);
    auto begin_binary =
        std::lower_bound(reduced.begin(), reduced.end(), l2, comp_size);
    auto end_binary =
        std::upper_bound(reduced.begin(), reduced.end(), l2, comp_size);
    std::unordered_set<std::pair<int, int>, PairHash> binaries;
    auto out_binary = begin_binary;
    for (; begin_binary != end_binary; ++begin_binary) {
        std::pair<int, int> p{(*begin_binary)[0], (*begin_binary)[1]};
        if (binaries.insert(p).second) {
            if (out_binary != begin_binary) {
                *out_binary = std::move(*begin_binary);
            }
            ++out_binary;
        }
    }
    if (out_binary != begin_binary) {
        reduced.erase(std::move(end_binary, reduced.end(), out_binary),
                      reduced.end());
    }
}

} // namespace impl

/**
 * @brief Reduce external clauses, removing duplicates, tautologies and subsumed
 * clauses.
 *
 * @param clauses
 * @return std::vector<std::vector<int>>
 */
inline std::vector<ExternalClause>
reduce_external_clauses(const std::vector<ExternalClause>& clauses) {
    std::vector<std::vector<int>> reduced;
    reduced.reserve(clauses.size());
    for (const auto& c : clauses) {
        auto res = reduce_external_clause(c);
        if (!res.empty()) {
            reduced.emplace_back(std::move(res));
        }
    }
    impl::reduced_remove_identical_pairs(reduced);
    return reduced;
}

} // namespace sammy

#endif
==> ./range.h <==
#ifndef SAMMY_RANGE_H_INCLUDED_
#define SAMMY_RANGE_H_INCLUDED_

#include <boost/iterator/iterator_facade.hpp>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <random>
#include <type_traits>
#include <utility>

namespace sammy {

/**
 * @brief A range spanned by two iterators.
 *
 * @tparam Iterator
 */
template <typename Iterator> class IteratorRange {
  public:
    IteratorRange(Iterator b, Iterator e) : m_beg(b), m_end(e) {}

    explicit IteratorRange(const std::pair<Iterator, Iterator>& p)
        : m_beg(p.first), m_end(p.second) {}

    Iterator begin() const noexcept { return m_beg; }
    Iterator end() const noexcept { return m_end; }
    std::size_t size() const noexcept { return std::distance(m_beg, m_end); }

  private:
    Iterator m_beg, m_end;
};

/**
 * @brief A range of random integers of a certain length.
 *
 * @tparam Value
 * @tparam RngType
 */
template <typename Value, typename RngType> class RandomIntRange {
    RngType* m_rng;
    std::size_t m_length;
    std::uniform_int_distribution<Value> m_dist;

  public:
    class Iterator {
        friend class RandomIntRange;
        RngType* m_rng;
        std::uniform_int_distribution<Value> m_dist;
        std::size_t m_offs;
        Value m_curr;

        Iterator(RngType* rng, std::uniform_int_distribution<Value> dist,
                 std::size_t offs, Value curr) noexcept
            : m_rng(rng), m_dist(dist), m_offs(offs), m_curr(curr) {}

      public:
        using iterator_category = std::input_iterator_tag;
        using value_type = Value;
        using difference_type = std::ptrdiff_t;
        using pointer = const Value*;
        using reference = Value;

        Iterator() = default;

        bool operator==(const Iterator& o) const noexcept {
            return m_offs == o.m_offs;
        }

        bool operator!=(const Iterator& o) const noexcept {
            return m_offs != o.m_offs;
        }

        Iterator& operator++() noexcept {
            ++m_offs;
            m_curr = m_dist(*m_rng);
            return *this;
        }

        Iterator operator++(int) noexcept {
            Iterator result(*this);
            ++*this;
            return result;
        }

        Value operator*() const noexcept { return m_curr; }

        const Value* operator->() const noexcept { return &m_curr; }
    };

    Iterator begin() { return Iterator(m_rng, m_dist, 0, m_dist(*m_rng)); }

    Iterator end() { return Iterator(m_rng, m_dist, m_length, 0); }

    RandomIntRange(RngType& rng, std::size_t length,
                   std::uniform_int_distribution<Value> dist = {})
        : m_rng(&rng), m_length(length), m_dist(dist) {}
};

template <typename ValueType, typename RngType>
inline auto random_int_range(std::size_t length, RngType& rng) noexcept {
    return RandomIntRange<ValueType, RngType>(rng, length);
}

template <typename ValueType, typename RngType>
inline auto random_int_range(
    std::size_t length, RngType& rng,
    const std::uniform_int_distribution<ValueType>& dist) noexcept {
    return RandomIntRange<ValueType, RngType>(rng, length, dist);
}

// deduction guide
template <typename Iterator1, typename Iterator2>
explicit IteratorRange(Iterator1, Iterator2)
    -> IteratorRange<std::remove_cv_t<std::remove_reference_t<Iterator1>>>;

template <typename RandomAccessIterator, typename Callable>
static inline void
split_range(RandomAccessIterator begin, RandomAccessIterator end,
            std::size_t max_partitions, Callable&& callback) {
    auto count = std::size_t(end - begin);
    if (max_partitions > count)
        max_partitions = count;
    std::size_t elements_per_partition = count / max_partitions;
    std::size_t rem_mod = count % max_partitions;
    RandomAccessIterator current = begin;
    for (std::size_t i = 0; i < max_partitions; ++i) {
        std::size_t current_count = elements_per_partition;
        if (rem_mod > 0) {
            --rem_mod;
            ++current_count;
        }
        RandomAccessIterator current_end = current + current_count;
        callback(RandomAccessIterator(current),
                 RandomAccessIterator(current_end));
        current = current_end;
    }
}

template <typename RNGType, typename Iterator>
class RandomSkipIterator
    : boost::iterator_facade<RandomSkipIterator<RNGType, Iterator>,
                             typename std::iterator_traits<Iterator>::value,
                             std::forward_iterator_tag,
                             typename std::iterator_traits<Iterator>::reference,
                             std::ptrdiff_t> {
  public:
    RandomSkipIterator(RNGType* rng, Iterator iter, Iterator end,
                       double prob) noexcept
        : rng(rng), iter(iter), end(end), prob(prob) {}

    RandomSkipIterator() noexcept : rng(nullptr), iter(), end(), prob(0.0) {}

  private:
    friend class boost::iterator_core_access;

    typename std::iterator_traits<Iterator>::reference
    dereference() const noexcept {
        return *iter;
    }

    bool equal(const RandomSkipIterator& other) const noexcept {
        return iter == other.iter;
    }

    void increment() noexcept {
        std::geometric_distribution<std::size_t> jump_dist(prob);
        std::size_t dist = jump_dist(*rng) + 1;
        if (std::size_t(std::distance(iter, end)) <= dist) {
            iter = end;
        } else {
            std::advance(iter, dist);
        }
    }

    RNGType* rng;
    Iterator iter;
    Iterator end;
    double prob;
};

} // namespace sammy

#endif
==> ./process.h <==
#ifndef SAMMY_PROCESS_H_INCLUDED_
#define SAMMY_PROCESS_H_INCLUDED_

#include <boost/filesystem.hpp>
#include <memory>
#include <optional>
#include <cstddef>
#include <cstdlib>
#include <vector>
#include <string>
#include <numeric>

namespace sammy {

/**
 * A very simplistic process management class,
 * since Boost.Process is a complete mess when it comes to managing and
 * waiting on processes.
 * It can't even wait on a process in a thread-safe manner,
 * even without any time limits, and with completely unrelated
 * boost::process::child objects.
 */
class Process {
  public:
    Process() noexcept;

    Process(const boost::filesystem::path& path,
            const std::vector<std::string>& args);

    Process(const Process&) = delete;
    Process& operator=(const Process&) = delete;

    Process(Process&&) noexcept;
    Process& operator=(Process&&) noexcept;

    /**
     * Destructor aborts if the process
     * is still not waited for at destruction time.
     */
    ~Process();

    /**
     * Write data to the process.
     * Throws an error on failure.
     */
    void write_to_process(const char* data, std::size_t size);

    /**
     * Send EOF to the process.
     */
    void send_eof();

    /**
     * Read data from the process.
     * Returns the number of bytes read;
     * only less than size on reaching EOF.
     */
    std::size_t read_from_process(char* data, std::size_t size);

    /**
     * Terminate the process, using
     * SIGKILL on unix-like systems and TerminateProcess on windows.
     */
    void terminate();

    /**
     * Wait for the process to exit and grab its exit code.
     */
    int wait();

    /**
     * Get the exit code of the process, if it has exited and
     * we have waited on it; otherwise, throws an error.
     */
    int exit_code() const;

    /**
     * Check if the process has ever been started.
     */
    bool valid() const noexcept { return m_impl != nullptr; }

    /**
     * Reset the process to represent not-a-process.
     */
    void reset();

    /**
     * Check if we have created and waited on the process.
     */
    bool have_waited() const noexcept;

  private:
    class Impl;

    std::unique_ptr<Impl> m_impl;
};

} // namespace sammy

#endif
==> ./time.h <==
#ifndef SAMMY_TIME_H_INCLUDED_
#define SAMMY_TIME_H_INCLUDED_

#include <chrono>

namespace sammy {

using Clock = std::chrono::steady_clock;

template <typename TP1, typename TP2>
inline double seconds_between(const TP1& before, const TP2& after) noexcept {
    return std::chrono::duration_cast<std::chrono::duration<double>>(after -
                                                                     before)
        .count();
}

} // namespace sammy

#endif
==> ./global_stat_info.h <==
#ifndef SAMMY_GLOBAL_STAT_INFO_H_INCLUDED_
#define SAMMY_GLOBAL_STAT_INFO_H_INCLUDED_

#include "literals.h"

namespace sammy {

struct GlobalStatInfo {
    mutable std::mutex lock;
    HashMap<std::string, double> double_stats;
    HashMap<std::string, std::int64_t> int_stats;

    void double_stat_add(const std::string& name, double value) {
        std::unique_lock l{lock};
        auto it = double_stats.find(name);
        if (it == double_stats.end()) {
            double_stats[name] = value;
        } else {
            it->second += value;
        }
    }

    void int_stat_add(const std::string& name, std::int64_t value) {
        std::unique_lock l{lock};
        auto it = int_stats.find(name);
        if (it == int_stats.end()) {
            int_stats[name] = value;
        } else {
            it->second += value;
        }
    }
};

inline std::ostream& operator<<(std::ostream& o, const GlobalStatInfo& l) {
    std::unique_lock lck{l.lock};
    o << "GlobalStatInfo{\n";
    for (const auto& [k, v] : l.double_stats) {
        o << "\t" << k << ": " << v << ",\n";
    }
    for (const auto& [k, v] : l.int_stats) {
        o << "\t" << k << ": " << v << ",\n";
    }
    o << "}";
    return o;
}

inline GlobalStatInfo& get_global_stats() {
    static GlobalStatInfo stats;
    return stats;
}

} // namespace sammy

#endif
==> ./randsat.h <==
#ifndef SAMMY_RANDSAT_H_INCLUDED_
#define SAMMY_RANDSAT_H_INCLUDED_

#include "error.h"
#include "rng.h"
#include "shared_db_propagator.h"
#include <chrono>

namespace sammy {

struct RandSATStats {
    double solve_seconds = 0.0;
    std::size_t propagations = 0;
    std::size_t conflicts = 0;
    std::size_t total_learned_clause_length = 0;
    std::size_t conflict_level_sum = 0;
    std::size_t total_decisions = 0;

    template <typename OutObject> void export_to(OutObject& object) const {
        object["solve_time"] = solve_seconds;
        object["propagations"] = propagations;
        object["conflicts"] = conflicts;
        object["total_learned_clause_length"] = total_learned_clause_length;
        object["conflict_level_sum"] = conflict_level_sum;
        object["total_decisions"] = total_decisions;
    }
};

template <typename RNGType> class RandSAT {
  public:
    struct TrivialClauseHandler {
        bool new_unary(Lit) { return true; }
        bool new_binary(Lit, Lit) { return true; }
        bool new_clause(ClauseDB::Lits) { return true; }
    };

    struct RecordSizeClauseHandler {
        bool new_unary(Lit) {
            that->stats->total_learned_clause_length += 1;
            return true;
        }

        bool new_binary(Lit, Lit) {
            that->stats->total_learned_clause_length += 2;
            return true;
        }

        bool new_clause(ClauseDB::Lits lits) {
            that->stats->total_learned_clause_length +=
                std::distance(lits.begin(), lits.end());
            return true;
        }

        RandSAT* that;
    };

    explicit RandSAT(RNGType& rng, const ClauseDB& clauses)
        : rng(&rng), clauses(clauses), view(&this->clauses),
          propagator(&this->clauses), stats(&stats_buffer),
          literal_dist(0, 2 * clauses.num_vars() - 1), handler{this} {
        TrivialClauseHandler trivial;
        view.handle_new_clauses(trivial);
    }

    void set_stats(RandSATStats* stats) noexcept { this->stats = stats; }

    std::vector<bool> solve() {
        auto before = std::chrono::steady_clock::now();
        while (p_push_next_choice()) {
            stats->total_decisions += 1;
        }
        auto after = std::chrono::steady_clock::now();
        stats->solve_seconds =
            std::chrono::duration_cast<std::chrono::duration<double>>(after -
                                                                      before)
                .count();
        std::vector<bool> result(clauses.num_vars(), false);
        for (Lit l : propagator.get_trail()) {
            result[lit::var(l)] = !lit::negative(l);
        }
        return result;
    }

  private:
    bool p_push_next_choice() {
        auto c = p_pick_next_choice();
        if (!c)
            return false;
        bool presult = propagator.push_level(*c);
        stats->propagations +=
            (propagator.get_trail().end() - propagator.current_level_begin()) -
            1;
        if (!presult) {
            stats->conflicts += 1;
            stats->conflict_level_sum += propagator.get_current_level();
            propagator.resolve_or_throw();
            view.handle_new_clauses(handler);
        }
        return true;
    }

    std::optional<Lit> p_pick_next_choice() {
        if (propagator.get_trail().size() == clauses.num_vars())
            return std::nullopt;
        for (int i = 0; i < 32; ++i) {
            Lit l = literal_dist(*rng);
            if (propagator.is_open(l)) {
                return l;
            }
        }
        open_buffer.clear();
        for (Var v = 0, nv = clauses.num_vars(); v < nv; ++v) {
            if (propagator.is_open(lit::positive_lit(v)))
                open_buffer.push_back(v);
        }
        std::uniform_int_distribution<std::size_t> range{0, open_buffer.size() -
                                                                1};
        Var open_var = open_buffer[range(*rng)];
        return std::uniform_int_distribution<int>{0, 1}(*rng)
                   ? lit::positive_lit(open_var)
                   : lit::negative_lit(open_var);
    }

    RNGType* rng;
    ClauseDB clauses;
    ClauseDBView view;
    SharedDBPropagator propagator;
    RandSATStats stats_buffer;
    RandSATStats* stats;
    std::uniform_int_distribution<Lit> literal_dist;
    std::vector<Var> open_buffer;
    RecordSizeClauseHandler handler;
};

template <class RNGType>
RandSAT(RNGType& rng, const ClauseDB& clauses) -> RandSAT<RNGType>;

inline std::vector<bool> randsolve(const ClauseDB& clauses) {
    RandSAT solver{sammy::rng(), clauses};
    return solver.solve();
}

} // namespace sammy

#endif
==> ./coloring.h <==
#ifndef SAMMY_COLORING_H_INCLUDED_
#define SAMMY_COLORING_H_INCLUDED_

#include "experiment_flags.h"
#include "kissat_solver.h"
#include "output.h"
#include "universe_subgraph.h"
#include <boost/iterator/counting_iterator.hpp>

namespace sammy {

inline std::pair<std::vector<SharedDBPropagator>, std::vector<Vertex>>
coloring_to_classes(const SharedDBPropagator& prop_template,
                    const std::vector<Vertex>& vertices,
                    const std::vector<std::size_t>& coloring,
                    std::size_t num_colors) {
    std::cout << "To classes: " << num_colors << std::endl;
    std::vector<std::vector<std::size_t>> color_classes(num_colors);
    std::vector<SharedDBPropagator> classes(num_colors, prop_template);
    std::vector<Vertex> spawners;
    for (std::size_t i = 0, n = vertices.size(); i < n; ++i) {
        auto c = coloring[i];
        color_classes[c].push_back(i);
    }
    for (std::size_t i = 0, n = color_classes.size(); i < n; ++i) {
        std::cout << "Class " << i << ": " << color_classes[i].size()
                  << std::endl;
    }
    for (std::size_t i = 0; i < num_colors; ++i) {
        auto& cc = color_classes[i];
        auto& ci = classes[i];
        spawners.push_back(vertices[cc.front()]);
        cc.erase(std::remove_if(cc.begin(), cc.end(),
                                [&](std::size_t vi) {
                                    return push_vertex(ci, vertices[vi]) >= 0;
                                }),
                 cc.end());
    }
    for (std::size_t i = 0; i < num_colors; ++i) {
        auto& cc = color_classes[i];
        for (std::size_t vi : cc) {
            bool found = false;
            for (std::size_t k = 0, cc = classes.size(); k < cc; ++k) {
                if (push_vertex(classes[k], vertices[vi]) >= 0) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                classes.push_back(prop_template);
                push_vertex(classes.back(), vertices[vi]);
                spawners.push_back(vertices[vi]);
            }
        }
    }
    return {std::move(classes), std::move(spawners)};
}

struct ColoringExtraConstraints {
    using ExtraConstraint = boost::container::small_vector<std::size_t, 4>;
    using ExtraConstraints = std::vector<ExtraConstraint>;

    ExtraConstraints constraints;
    std::vector<std::vector<std::size_t>> extra_constraints_by_vertex;

    void update_vertex_count(const UniverseSubgraph& subgraph) {
        extra_constraints_by_vertex.resize(subgraph.n());
    }

    void add_constraint(const ExtraConstraint& extra_constraint) {
        constraints.push_back(extra_constraint);
        p_update_extra_constraints_by_vertex();
    }

    void add_constraint(ExtraConstraint&& extra_constraint) {
        constraints.emplace_back(std::move(extra_constraint));
        p_update_extra_constraints_by_vertex();
    }

    template <typename Iterator>
    void add_constraint(Iterator begin, Iterator end) {
        constraints.emplace_back(begin, end);
        p_update_extra_constraints_by_vertex();
    }

    std::vector<SharedDBPropagator>
    coloring_to_classes(SharedDBPropagator& template_propagator,
                        const UniverseSubgraph& subgraph,
                        const std::vector<std::size_t>& coloring,
                        std::size_t num_colors) {
        std::vector<std::vector<std::size_t>> color_classes(num_colors);
        for (std::size_t vi = 0; vi < coloring.size(); ++vi) {
            color_classes[coloring[vi]].push_back(vi);
        }
        std::vector<SharedDBPropagator> result;
        for (std::size_t cc = 0; cc < num_colors; ++cc) {
            auto& cclass = color_classes[cc];
            if (cclass.empty())
                continue;
            template_propagator.reset_or_throw();
            for (auto vi_it = cclass.begin(), vi_end = cclass.end();
                 vi_it != vi_end; ++vi_it)
            {
                Vertex v = subgraph.vertex(*vi_it);
                if (push_vertex(template_propagator, v) < 0) {
                    p_constraint_from_prop(template_propagator, cclass.begin(),
                                           vi_it + 1, subgraph);
                    p_coloring_to_classes_failed(template_propagator, subgraph,
                                                 color_classes, cc + 1);
                    result.clear();
                    return result;
                }
            }
            result.push_back(template_propagator);
        }
        return result;
    }

  private:
    void p_coloring_to_classes_failed(
        SharedDBPropagator& template_propagator,
        const UniverseSubgraph& subgraph,
        const std::vector<std::vector<std::size_t>>& color_classes,
        std::size_t start_with_class) {
        for (std::size_t i = start_with_class; i < color_classes.size(); ++i) {
            auto& cclass = color_classes[i];
            if (cclass.size() <= 2)
                continue;
            template_propagator.reset_or_throw();
            for (auto vi_it = cclass.begin(), vi_end = cclass.end();
                 vi_it != vi_end; ++vi_it)
            {
                Vertex v = subgraph.vertex(*vi_it);
                if (push_vertex(template_propagator, v) < 0) {
                    p_constraint_from_prop(template_propagator, cclass.begin(),
                                           vi_it + 1, subgraph);
                    break;
                }
            }
        }
    }

    template <typename RandomAccessIterator>
    void p_constraint_from_prop(SharedDBPropagator& propagator,
                                RandomAccessIterator pushed_begin,
                                RandomAccessIterator pushed_end,
                                const UniverseSubgraph& subgraph) {
        // pushing this vertex failed;
        // the propagator has all previous vertices pushed
        std::size_t failed_index = *(pushed_end - 1);
        Vertex failed = subgraph.vertex(failed_index);
        if (propagator.is_false(failed.first)) {
            // case 1a: the first literal of the vertex is false
            p_learn_constraint(propagator, pushed_begin, pushed_end - 1,
                               lit::negate(failed.first), failed_index,
                               subgraph);
            return;
        }
        if (propagator.is_false(failed.second)) {
            // case 1b: the second literal of the vertex is false
            p_learn_constraint(propagator, pushed_begin, pushed_end - 1,
                               lit::negate(failed.second), failed_index,
                               subgraph);
            return;
        }
        if (propagator.is_open(failed.first)) {
            if (!propagator.push_level(failed.first)) {
                // case 2a: pushing the first literal from the vertex causes a
                // conflict
                p_learn_constraint_from_conflict(propagator, pushed_begin,
                                                 pushed_end, subgraph);
                return;
            }
            if (propagator.is_false(failed.second)) {
                // case 3: pushing the first literal from the vertex causes no
                // conflict, but implies the negation of the second literal
                // TODO: this case needs more thought...
                add_constraint(pushed_begin, pushed_end);
                return;
            }
        }
        // case 2b: pushing the second literal from the vertex causes a conflict
        if (propagator.is_open(failed.second) &&
            !propagator.push_level(failed.second))
        {
            p_learn_constraint_from_conflict(propagator, pushed_begin,
                                             pushed_end, subgraph);
            return;
        }
        throw std::logic_error("Invalid case in p_constraint_from_prop!");
    }

    template <typename RandomAccessIterator>
    std::vector<DynamicBitset>
    p_compute_vertex_consequences(const SharedDBPropagator& propagator_template,
                                  RandomAccessIterator vbegin,
                                  RandomAccessIterator vend,
                                  const UniverseSubgraph& subgraph) {
        SharedDBPropagator propagator{propagator_template};
        std::size_t n = propagator.db().num_vars();
        std::vector<DynamicBitset> vertex_consequences;
        for (std::size_t vindex : IteratorRange{vbegin, vend}) {
            propagator.reset_or_throw();
            Vertex v = subgraph.vertex(vindex);
            if (propagator.is_open(v.first))
                propagator.push_level(v.first);
            if (propagator.is_open(v.second))
                propagator.push_level(v.second);
            vertex_consequences.emplace_back(2 * n, false);
            auto& b = vertex_consequences.back();
            for (Lit l : propagator.get_trail()) {
                b[l] = true;
            }
        }
        return vertex_consequences;
    }

    template <typename RandomAccessIterator>
    void p_learn_constraint_from_conflict(SharedDBPropagator& propagator,
                                          RandomAccessIterator pushed_begin,
                                          RandomAccessIterator pushed_end,
                                          const UniverseSubgraph& subgraph) {
        auto vertex_consequences = p_compute_vertex_consequences(
            propagator, pushed_begin, pushed_end, subgraph);
        auto [conflict_literal, conflict_reason] = propagator.get_conflict();
        DynamicBitset active(2 * propagator.db().num_vars(), false);
        Lit coverage_block = lit::negate(conflict_literal);
        active[coverage_block] = true;
        for (Lit l : conflict_reason.lits(propagator.db())) {
            if (l != conflict_literal) {
                auto lneg = lit::negate(l);
                active[lneg] = true;
            }
        }
        auto best_cover = p_walk_trail(
            propagator, propagator.get_trail().size() - 1, coverage_block,
            active, vertex_consequences, pushed_begin, pushed_end);
        add_constraint(best_cover);
    }

    template <typename RandomAccessIterator>
    ExtraConstraint
    p_walk_trail(SharedDBPropagator& propagator, std::size_t start_index,
                 Lit initial_blocker, DynamicBitset& active,
                 const std::vector<DynamicBitset>& vertex_consequences,
                 RandomAccessIterator pushed_begin,
                 RandomAccessIterator pushed_end) {
        ExtraConstraint best_cover;
        const auto& trail = propagator.get_trail();
        const auto& reasons = propagator.get_reasons();
        Lit coverage_blocker = initial_blocker;
        DynamicBitset tmp;
        DynamicBitset coverable = vertex_consequences.front();
        std::for_each(vertex_consequences.begin() + 1,
                      vertex_consequences.end(),
                      [&](const DynamicBitset& b) { coverable |= b; });
        for (std::size_t tindex = start_index; tindex <= start_index; --tindex)
        {
            Lit tlit = trail[tindex];
            if (!active[tlit])
                continue;
            Reason r = reasons[tindex];
            switch (r.reason_length) {
            case 0:
                break;
            case 1:
                active[tlit] = false;
                break;
            default:
                for (Lit l : r.lits(propagator.db())) {
                    if (l != tlit) {
                        active[lit::negate(l)] = true;
                    }
                }
                active[tlit] = false;
                break;
            }
            if (!active[coverage_blocker]) {
                tmp = active;
                tmp -= coverable;
                auto blocker = tmp.ones_begin();
                if (blocker == tmp.ones_end()) {
                    p_compute_cover(active, tmp, vertex_consequences,
                                    pushed_begin, pushed_end, best_cover);
                    if (best_cover.size() <= 2)
                        break;
                } else {
                    coverage_blocker = Lit(*blocker);
                }
            }
        }
        if (best_cover.empty())
            throw std::logic_error(
                "Could not cover conflict literals with vertices!");
        return best_cover;
    }

    template <typename RandomAccessIterator>
    void p_learn_constraint(SharedDBPropagator& propagator,
                            RandomAccessIterator pushed_begin,
                            RandomAccessIterator pushed_end,
                            Lit problem_literal, std::size_t extra_vertex,
                            const UniverseSubgraph& subgraph) {
        auto vertex_consequences = p_compute_vertex_consequences(
            propagator, pushed_begin, pushed_end, subgraph);
        DynamicBitset active(2 * propagator.db().num_vars(), false);
        active[problem_literal] = true;
        std::size_t problem_index = propagator.get_trail_index(problem_literal);
        auto best_cover =
            p_walk_trail(propagator, problem_index, problem_literal, active,
                         vertex_consequences, pushed_begin, pushed_end);
        best_cover.push_back(extra_vertex);
        add_constraint(best_cover);
    }

    /**
     * Cover the current conflict literals
     * with vertices in a greedy fashion
     * (their decisions & consequences).
     */
    template <typename RandomAccessIterator>
    void p_compute_cover(DynamicBitset& active, DynamicBitset& tmp,
                         const std::vector<DynamicBitset>& vconsq,
                         RandomAccessIterator begin, RandomAccessIterator end,
                         ExtraConstraint& best_cover) {
        DynamicBitset still_active = active;
        ExtraConstraint cover;
        while (still_active.any()) {
            std::size_t best_num_covered = 0;
            std::size_t best_vertex = 0;
            std::size_t pindex = 0, best_pindex = 0;
            for (RandomAccessIterator p = begin; p != end; ++p, ++pindex) {
                tmp = still_active;
                tmp &= vconsq[pindex];
                auto count = tmp.count();
                if (count > best_num_covered) {
                    best_num_covered = count;
                    best_vertex = *p;
                    best_pindex = pindex;
                }
            }
            cover.push_back(best_vertex);
            still_active -= vconsq[best_pindex];
        }
        if (best_cover.empty() || cover.size() < best_cover.size()) {
            best_cover = std::move(cover);
        }
    }

    void p_update_extra_constraints_by_vertex() {
        ExtraConstraint& last = constraints.back();
        std::cout << "CONSTRAINT ADDED:";
        for (std::size_t vi : last) {
            std::cout << " " << vi;
        }
        std::cout << std::endl;
        std::size_t idx = constraints.size() - 1;
        for (std::size_t v : last) {
            extra_constraints_by_vertex[v].push_back(idx);
        }
    }
};

class SATKColoringSolver {
  public:
    SATKColoringSolver(UniverseSubgraph* subgraph, EventRecorder* recorder,
                       std::vector<std::size_t> clique, std::size_t k)
        : m_subgraph(subgraph), m_recorder(recorder),
          m_clique(std::move(clique)), m_k(k) {}

    void set_extra_constraints(ColoringExtraConstraints* extra) {
        m_extra_cons = extra;
    }

    std::optional<bool>
    solve(double time_limit = std::numeric_limits<double>::infinity()) {
        m_recorder->store_event("BEGIN_SOLVE_K_COLORING",
                                {{"k", m_k}, {"n", m_subgraph->n()}}, "k", "n");
        p_reduce_graph();
        std::vector<std::size_t> core_coloring;
        if (!m_new_to_old.empty()) {
            std::vector<bool> solution;
            auto optres = p_sat_solve(solution, time_limit);
            if (!optres) {
                m_recorder->store_event("DONE_SOLVE_K_COLORING",
                                        {{"result", "aborted"}}, "result");
                return std::nullopt;
            }
            m_recorder->store_event(
                "DONE_SOLVE_K_COLORING",
                {{"result", *optres ? "coloring" : "bound"}}, "result");
            if (*optres) {
                core_coloring = p_core_coloring_from_solution(solution);
            } else {
                return false;
            }
        }
        p_compact_core_coloring(core_coloring);
        p_extend_core_coloring(core_coloring);
        return true;
    }

    const std::vector<std::size_t>& get_coloring() const noexcept {
        return m_coloring;
    }

  private:
    std::vector<std::size_t>
    p_core_coloring_from_solution(const std::vector<bool>& sol) {
        std::size_t nnew = m_new_to_old.size();
        std::vector<std::size_t> result(
            nnew, std::numeric_limits<std::size_t>::max());
        std::size_t offs = 0;
        for (std::size_t vi = 0; vi < nnew; ++vi) {
            for (std::size_t c = 0; c < m_k; ++c, ++offs) {
                if (sol[offs]) {
                    result[vi] = c;
                    offs += m_k - c;
                    break;
                }
            }
        }
        return result;
    }

    auto p_var_range(std::size_t new_vertex_index) {
        using Lit = KissatSolver::Lit;
        return IteratorRange{
            boost::make_counting_iterator(Lit(new_vertex_index * m_k + 1)),
            boost::make_counting_iterator(
                Lit(new_vertex_index * m_k + m_k + 1))};
    }

    void p_add_clique(KissatSolver& solver, std::size_t nnew) {
        std::size_t clique_next_color = 0;
        for (std::size_t cvi : m_clique) {
            std::size_t new_i = m_old_to_new[cvi];
            if (new_i < nnew) {
                solver.add_short_clause(
                    p_var_range(new_i).begin()[clique_next_color++]);
            }
        }
    }

    void p_add_edges(KissatSolver& solver, std::size_t nold, std::size_t nnew) {
        using Lit = KissatSolver::Lit;
        for (std::size_t i = 0; i < nold; ++i) {
            std::size_t new_i = m_old_to_new[i];
            if (new_i >= nnew)
                continue;
            auto vri = p_var_range(new_i);
            for (std::size_t j : m_subgraph->matrix_row(i).ones_from(i + 1)) {
                std::size_t new_j = m_old_to_new[j];
                if (new_j >= nnew)
                    continue;
                auto vrj = p_var_range(new_j);
                auto vi = vri.begin();
                for (Lit lj : vrj) {
                    solver.add_short_clause(-lj, -*vi);
                    ++vi;
                }
            }
        }
    }

    void p_add_extra(KissatSolver& solver) {
        ColoringExtraConstraints::ExtraConstraint extra;
        for (std::size_t i = 0, ncons = m_extra_cons->constraints.size();
             i < ncons; ++i)
        {
            if (!m_extra_cons_active[i])
                continue;
            extra.clear();
            const auto& cons = m_extra_cons->constraints[i];
            std::transform(
                cons.begin(), cons.end(), std::back_inserter(extra),
                [&](std::size_t vold) { return m_old_to_new[vold]; });
            for (std::size_t c = 0; c < m_k; ++c) {
                std::for_each(
                    extra.begin(), extra.end(), [&](std::size_t vnew) {
                        solver.add_literal(-p_var_range(vnew).begin()[c]);
                    });
                solver.finish_clause();
            }
        }
    }

    std::optional<bool> p_sat_solve(std::vector<bool>& solution,
                                    double time_limit) {
        std::size_t nnew = m_new_to_old.size();
        std::size_t nold = m_subgraph->n();
        KissatSolver solver;
        solver.reserve(m_k * nnew);
        for (std::size_t i = 0; i < nnew; ++i) {
            auto r = p_var_range(i);
            solver.add_clause(r.begin(), r.end());
            for (auto c1i = r.begin(), e = r.end(); c1i != e; ++c1i) {
                for (auto c2i = c1i + 1; c2i != e; ++c2i) {
                    solver.add_short_clause(-*c1i, -*c2i);
                }
            }
        }
        p_add_clique(solver, nnew);
        p_add_edges(solver, nold, nnew);
        if (m_extra_cons)
            p_add_extra(solver);
        auto r = solver.solve(time_limit);
        if (r && *r)
            solution = std::move(solver.get_model().raw());
        return r;
    }

    void p_compact_core_coloring(std::vector<std::size_t>& core_coloring) {
        DynamicBitset unused(m_k, true);
        for (std::size_t c : core_coloring) {
            unused[c].reset();
        }
        std::vector<std::size_t> uco;
        std::copy(unused.ones_begin(), unused.ones_end(),
                  std::back_inserter(uco));
        if (uco.empty())
            return;
        std::vector<std::size_t> remap;
        auto uco_iter = uco.begin();
        std::size_t next_color = 0;
        for (std::size_t c = 0; c < m_k; ++c) {
            if (c == *uco_iter) {
                remap.push_back(0);
                if (++uco_iter == uco.end()) {
                    for (++c; c < m_k; ++c) {
                        remap.push_back(next_color++);
                    }
                    break;
                }
            } else {
                remap.push_back(next_color++);
            }
        }
        for (std::size_t& c : core_coloring) {
            c = remap[c];
        }
    }

    void p_extend_core_coloring(const std::vector<std::size_t>& core_coloring) {
        m_coloring.assign(m_subgraph->n(), 0);
        for (std::size_t i = 0, nnew = core_coloring.size(); i < nnew; ++i) {
            std::size_t old_index = m_new_to_old[i];
            m_coloring[old_index] = core_coloring[i];
        }
        DynamicBitset avail(m_k, true);
        for (std::size_t deleted :
             IteratorRange{m_deletion_order.rbegin(), m_deletion_order.rend()})
        {
            avail.assign(m_k, true);
            for (std::size_t neigh : m_subgraph->matrix_row(deleted).ones()) {
                if (m_degrees[neigh] >= m_k) {
                    avail[m_coloring[neigh]].reset();
                }
            }
            if (m_extra_cons) {
                for (std::size_t neigh : m_extra_pseudo_edges[deleted]) {
                    if (m_degrees[neigh] >= m_k) {
                        avail[m_coloring[neigh]].reset();
                    }
                }
            }
            assert(avail.count() != 0);
            m_coloring[deleted] = *avail.ones_begin();
            m_degrees[deleted] =
                m_k; // make sure the new vertex now counts as present
        }
    }

    bool p_degree_too_low(std::size_t vindex) const noexcept {
        return m_degrees[vindex] < m_k;
    }

    void p_reduce_graph() {
        p_compute_degrees();
        p_init_deletion_order();
        if (m_deletion_order.empty()) {
            m_new_to_old.assign(boost::make_counting_iterator(std::size_t(0)),
                                boost::make_counting_iterator(m_subgraph->n()));
            m_old_to_new = m_new_to_old;
            m_extra_cons_active.assign(m_extra_cons->constraints.size(), true);
            m_extra_pseudo_edges.assign(
                m_subgraph->n(),
                boost::container::small_vector<std::size_t, 4>{});
        } else {
            p_bfs_deletion_order();
            p_extract_new_to_old();
        }
    }

    void p_init_deletion_order() {
        m_deletion_order.clear();
        std::copy_if(boost::make_counting_iterator(std::size_t(0)),
                     boost::make_counting_iterator(m_subgraph->n()),
                     std::back_inserter(m_deletion_order),
                     [&](std::size_t vi) { return p_degree_too_low(vi); });
    }

    void p_bfs_deletion_order() {
        if (m_extra_cons) {
            m_extra_cons_active.assign(m_extra_cons->constraints.size(), true);
            m_extra_pseudo_edges.assign(
                m_subgraph->n(),
                boost::container::small_vector<std::size_t, 4>{});
        }
        std::size_t queue_pos = 0;
        while (queue_pos < m_deletion_order.size()) {
            std::size_t deleted = m_deletion_order[queue_pos++];
            for (std::size_t neigh : m_subgraph->matrix_row(deleted).ones()) {
                if (--m_degrees[neigh] == m_k - 1) {
                    m_deletion_order.push_back(neigh);
                }
            }
            if (m_extra_cons) {
                for (std::size_t cons_idx :
                     m_extra_cons->extra_constraints_by_vertex[deleted])
                {
                    if (m_extra_cons_active[cons_idx]) {
                        p_deactivate_cons(deleted, cons_idx);
                    }
                }
            }
        }
    }

    void p_deactivate_cons(std::size_t deleted_vertex,
                           std::size_t deactivated_cons) {
        m_extra_cons_active[deactivated_cons] = false;
        const auto& cons = m_extra_cons->constraints[deactivated_cons];
        std::size_t other_vertex = cons[0];
        if (other_vertex == deleted_vertex) {
            other_vertex = cons[1];
        }
        m_extra_pseudo_edges[deleted_vertex].push_back(other_vertex);
        m_extra_pseudo_edges[other_vertex].push_back(deleted_vertex);
        for (const auto& vertex : cons) {
            if (--m_degrees[vertex] == m_k - 1) {
                m_deletion_order.push_back(vertex);
            }
        }
    }

    void p_extract_new_to_old() {
        m_new_to_old.clear();
        m_old_to_new.clear();
        std::size_t num_used = 0;
        for (std::size_t i = 0, n = m_subgraph->n(); i < n; ++i) {
            if (m_degrees[i] >= m_k) {
                m_new_to_old.push_back(i);
                m_old_to_new.push_back(num_used++);
            } else {
                m_old_to_new.push_back(std::numeric_limits<std::size_t>::max());
            }
        }
    }

    void p_compute_degrees() {
        m_degrees = m_subgraph->get_degrees();
        if (m_extra_cons) {
            for (std::size_t i = 0, n = m_subgraph->n(); i < n; ++i) {
                m_degrees[i] +=
                    m_extra_cons->extra_constraints_by_vertex[i].size();
            }
        }
    }

    UniverseSubgraph* m_subgraph;
    EventRecorder* m_recorder;
    std::vector<std::size_t> m_clique;
    std::vector<std::size_t> m_deletion_order;
    std::vector<std::size_t> m_new_to_old;
    std::vector<std::size_t> m_old_to_new;
    std::vector<std::size_t> m_degrees;
    std::vector<std::size_t> m_coloring;
    std::size_t m_k;
    ColoringExtraConstraints* m_extra_cons = nullptr;
    std::vector<bool> m_extra_cons_active;
    std::vector<boost::container::small_vector<std::size_t, 4>>
        m_extra_pseudo_edges;
};

class SATColoringSolver {
  public:
    SATColoringSolver(UniverseSubgraph* subgraph, EventRecorder* recorder,
                      std::vector<Vertex> clique,
                      std::vector<std::size_t> best_known_coloring)
        : m_subgraph(subgraph), m_recorder(recorder),
          best_known_lower_bound(clique.size()), best_clique(std::move(clique)),
          best_coloring(std::move(best_known_coloring)),
          best_num_colors(
              *std::max_element(best_coloring.begin(), best_coloring.end()) +
              1) {
        clique_indices.reserve(best_clique.size());
        std::transform(best_clique.begin(), best_clique.end(),
                       std::back_inserter(clique_indices),
                       [&](Vertex v) { return m_subgraph->vertex_index(v); });
        best_num_colors = p_compact_coloring(best_coloring, best_num_colors);
    }

    void added_vertices() {
        const auto old_n = best_coloring.size();
        const auto new_n = m_subgraph->n();
        if (old_n >= new_n)
            return;
        for (std::size_t i = old_n; i < new_n; ++i) {
        }
    }

    std::optional<bool> run_with_num_colors(
        std::size_t num_colors,
        double time_limit = std::numeric_limits<double>::infinity()) {
        if (num_colors >= best_num_colors)
            return true;
        if (num_colors < best_known_lower_bound)
            return false;
        const std::size_t n = m_subgraph->n();
        auto x_vc = [&](std::size_t vindex, std::size_t color) -> Lit {
            return Lit(color * n + vindex + 1);
        };
        KissatSolver solver;
        solver.reserve(num_colors * n);
        for (std::size_t i = 0; i < clique_indices.size(); ++i) {
            solver.add_short_clause(x_vc(clique_indices[i], i));
        }
        std::vector<Lit> clause_buffer(num_colors, 0);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j : m_subgraph->matrix_row(i).ones_from(i + 1)) {
                for (std::size_t k = 0; k < num_colors; ++k) {
                    solver.add_short_clause(-x_vc(j, k), -x_vc(i, k));
                }
            }
            for (std::size_t k = 0; k < num_colors; ++k) {
                clause_buffer[k] = x_vc(i, k);
            }
            solver.add_clause(clause_buffer.begin(), clause_buffer.end());
        }
        auto result = solver.solve(time_limit);
        if (!result)
            return std::nullopt;
        if (!*result) {
            m_recorder->store_event(
                "NEW_LOWER_BOUND",
                {{"prev", best_known_lower_bound}, {"new", num_colors + 1}},
                "prev", "new");
            best_known_lower_bound = num_colors + 1;
            return false;
        }
        auto model = solver.get_model();
        const auto& solution = model.raw();
        m_recorder->store_event(
            "NEW_COLORING", {{"prev", best_num_colors}, {"new", num_colors}},
            "prev", "new");
        for (std::size_t k = 0, offs = 0; k < num_colors; ++k) {
            for (std::size_t v = 0; v < n; ++v, ++offs) {
                if (solution[offs]) {
                    best_coloring[v] = k;
                }
            }
        }
        best_num_colors = p_compact_coloring(best_coloring, num_colors);
        return true;
    }

    bool optimize_coloring(
        double time_limit = std::numeric_limits<double>::infinity()) {
        if (time_limit <= 0)
            return false;
        auto begin = Clock::now();
        while (best_known_lower_bound < best_num_colors) {
            std::size_t mid = (best_known_lower_bound + best_num_colors) / 2;
            double trem = time_limit - seconds_between(begin, Clock::now());
            m_recorder->store_event("BEGIN_SAT_COLORING",
                                    {{"best_lb", best_known_lower_bound},
                                     {"best_coloring", best_num_colors},
                                     {"query", mid},
                                     {"n", m_subgraph->n()}},
                                    "best_lb", "best_coloring", "query", "n");
            auto result = run_with_num_colors(mid, trem);
            if (!result) {
                m_recorder->store_event("SAT_COLORING_TIMEOUT");
                return false;
            }
            m_recorder->store_event("DONE_SAT_COLORING",
                                    {{"best_lb", best_known_lower_bound},
                                     {"best_coloring", best_num_colors}},
                                    "best_lb", "best_coloring");
        }
        m_recorder->store_event("OPTIMUM_SUBGRAPH_COLORING",
                                {{"best_coloring", best_num_colors}},
                                "best_coloring");
        return true;
    }

    std::pair<std::vector<SharedDBPropagator>, std::vector<Vertex>>
    coloring_to_classes() {
        auto& prop = m_subgraph->get_propagator();
        prop.reset_or_throw();
        return sammy::coloring_to_classes(prop, m_subgraph->vertex_set(),
                                          best_coloring, best_num_colors);
    }

  private:
    std::size_t p_compact_coloring(std::vector<std::size_t>& coloring,
                                   std::size_t num_colors) {
        const auto n = coloring.size();
        std::vector<std::vector<std::size_t>> color_classes(num_colors);
        for (std::size_t i = 0; i < n; ++i) {
            color_classes[coloring[i]].push_back(i);
        }
        std::size_t num_used = 0;
        std::vector<DynamicBitset> vertex_compatibility(num_colors,
                                                        DynamicBitset(n, true));
        for (std::size_t k = num_colors - 1; k < num_colors; --k) {
            const auto& cclass = color_classes[k];
            for (std::size_t v : cclass) {
                for (std::size_t cand = 0; cand < num_colors; ++cand) {
                    if (vertex_compatibility[cand][v]) {
                        if (cand >= num_used)
                            num_used = cand + 1;
                        vertex_compatibility[cand] -= m_subgraph->matrix_row(v);
                        coloring[v] = cand;
                        break;
                    }
                }
            }
        }
        return num_used;
    }

    UniverseSubgraph* m_subgraph;
    EventRecorder* m_recorder;
    std::size_t best_known_lower_bound;
    std::vector<Vertex> best_clique;
    std::vector<std::size_t> clique_indices;
    std::vector<std::size_t> best_coloring;
    std::size_t best_num_colors;
};

} // namespace sammy

#endif
==> ./memory_usage.h <==
#ifndef BARRAGE_MEMORY_USAGE_H_INCLUDED_
#define BARRAGE_MEMORY_USAGE_H_INCLUDED_

#include <cstddef>
#include <cstdint>
#include <cstdlib>

namespace sammy {

/**
 * @brief Get the current peak resident set size (RSS) in bytes,
 * i.e., the maximum RSS since the start of the process, so far.
 */
inline std::size_t current_peak_rss();

} // namespace sammy

#if defined(__APPLE__) && defined(__MACH__)

#include <sys/resource.h>

std::size_t sammy::current_peak_rss() {
    struct rusage u;
    if (getrusage(RUSAGE_SELF, &u) < 0) {
        std::abort();
    }
    return u.ru_maxrss;
}

#elif __has_include(<sys/resource.h>)

#include <sys/resource.h>

std::size_t sammy::current_peak_rss() {
    struct rusage u;
    if (getrusage(RUSAGE_SELF, &u) < 0) {
        std::abort();
    }
    return u.ru_maxrss * 1024;
}

#elif __has_include(<psapi.h>) && __has_include(<windows.h>) && defined(WIN32)

#define WIN32_LEAN_AND_MEAN
#include <psapi.h>
#include <windows.h>

std::size_t sammy::current_peak_rss() {
    PROCESS_MEMORY_COUNTERS pmc;
    if (!GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
        std::abort();
    }
    return pmc.PeakWorkingSetSize;
}

#endif

#endif
==> ./incremental_sat_lns.h <==
#ifndef SAMMY_INCREMENTAL_SAT_LNS_H_INCLUDED_
#define SAMMY_INCREMENTAL_SAT_LNS_H_INCLUDED_

#include "barrage.h"
#include "barrage_lns_subproblem.h"
#include "clause_db.h"
#include "lazy_g2_adjacency_matrix.h"
#include "output.h"
#include "shared_db_propagator.h"
#include <algorithm>
#include <variant>

namespace sammy {

/**
 * Generalized subproblem solver interface:
 *  - creation from LNSSubproblem, Propagator and EventRecorder*,
 *  - interruptible solve method,
 *  - solution extraction method,
 *  - abort method,
 *  - MES extraction method.
 */

template <typename IncrementalSatSolver>
class FixedMESIncrementalSATImprovementSolver {
  public:
    using SatSolver = IncrementalSatSolver;
    using Lit = typename SatSolver::Lit;
    using SLit = sammy::Lit;
    using LitOrVal = std::variant<Lit, bool>;

    static std::string name() {
        return std::string("IncrementalSAT<") + SatSolver::name() + ">";
    }

  private:
    // basic incremental SAT-solver backend
    SatSolver m_solver;
    // subproblem information
    LNSSubproblem m_subproblem;
    // event recorder
    EventRecorder* m_recorder;
    // worker id for events
    std::size_t m_worker_id;
    // propagator for the SAT model
    SharedDBPropagator m_propagator;
    // flag to set when, during construction, we find infeasibility
    bool m_infeasible_by_construction;
    // contiguous array of all second literals of vertices
    std::vector<SLit> m_literal_partners_array;
    // bitmap that tracks, for each entry in m_literal_partners_array, whether
    // the corresponding vertex is covered.
    DynamicBitset m_covered_literal_partners;
    // tracks which literal partner entries are explicitly covered by
    // clauses/variables
    DynamicBitset m_explicitized_literal_partners;
    // how many vertices are explicitly covered?
    std::size_t m_num_explicitly_covered = 0;
    // for each first literal of a vertex, a [begin_index, end_index) of its
    // second literals as they appear in m_literal_partners_array
    std::vector<std::pair<std::size_t, std::size_t>> m_literal_partners_of;
    // (outer) index k: configuration index [0, allowed_configurations());
    // inner index: feature index [0, propagator().db().num_vars())
    std::vector<std::vector<LitOrVal>> m_config_vars;
    // current solution: for each configuration, assignment of variables
    std::vector<DynamicBitset> m_current_solution;
    // buffer for solver literals/clauses
    std::vector<Lit> m_buffer;
    // buffer for selecting uncovered vertices to be explicitized
    std::vector<std::pair<Vertex, std::size_t>> m_selected_uncovered;
    // set by abort method in addition to aborting the SAT solver itself
    std::atomic<bool> m_abort{false};

  public:
    FixedMESIncrementalSATImprovementSolver(PortfolioSolver* portfolio,
                                            LNSSubproblem&& subproblem,
                                            SharedDBPropagator prop,
                                            EventRecorder* recorder,
                                            std::size_t worker_id)
        : m_solver(), m_subproblem(std::move(subproblem)), m_recorder(recorder),
          m_worker_id(worker_id), m_propagator(std::move(prop)),
          m_infeasible_by_construction(allowed_configurations() < mes_size()) {
        try {
            if (m_infeasible_by_construction)
                return;
            p_prepare_vertex_set();
            throw_if_interrupted();
            p_make_literal_partners();
            throw_if_interrupted();
            p_make_config_vars();
            if (m_infeasible_by_construction)
                return;
            throw_if_interrupted();
            p_make_vertical_clauses();
            throw_if_interrupted();
        } catch (...) {
            subproblem = std::move(m_subproblem);
            throw;
        }
    }

    std::size_t allowed_configurations() const noexcept {
        return m_subproblem.removed_configurations.size() - 1;
    }

    std::size_t mes_size() const noexcept {
        return m_subproblem.mutually_exclusive_set.size();
    }

    std::size_t num_uncovered() const noexcept {
        return m_subproblem.uncovered_universe.size();
    }

    const std::vector<Vertex>& mes_vertices() const noexcept {
        return m_subproblem.mutually_exclusive_set;
    }

    std::optional<bool> solve() {
        if (m_infeasible_by_construction) {
            p_store_event("INCREMENTAL_SAT_INFEASIBLE_BY_CONSTRUCTION");
            return false;
        }
        for (;;) {
            OutputObject event_data{
                {"universe_size", num_uncovered()},
                {"num_removed", m_subproblem.removed_configurations.size()},
                {"mes_size", mes_size()},
                {"num_explicitly_covered", m_num_explicitly_covered}};
            p_store_event("INCREMENTAL_SAT_LNS_BEGIN_SAT_SOLVE", event_data,
                          "universe_size", "num_removed", "mes_size",
                          "num_explicitly_covered");
            std::optional<bool> res = p_solve_iteration();
            if (!res) {
                p_store_event("INCREMENTAL_SAT_LNS_SOLUTION_ABORTED",
                              event_data, "universe_size", "num_removed",
                              "mes_size", "num_explicitly_covered");
                return std::nullopt;
            }
            if (!*res) {
                p_store_event("INCREMENTAL_SAT_LNS_SOLUTION_WAS_OPTIMAL",
                              event_data, "universe_size", "num_removed",
                              "mes_size", "num_explicitly_covered");
                return false;
            }
            p_store_event("INCREMENTAL_SAT_LNS_SATISFIABLE", event_data,
                          "universe_size", "num_removed", "mes_size",
                          "num_explicitly_covered");
            p_extract_solution();
            p_identify_covered();
            p_sample_uncovered();
            if (m_selected_uncovered.empty()) {
                p_store_event("INCREMENTAL_SAT_LNS_SOLUTION_IMPROVED",
                              event_data, "universe_size", "num_removed",
                              "mes_size", "num_explicitly_covered");
                return true;
            }
            for (auto entry : m_selected_uncovered) {
                p_make_explicit_unfixed(entry.second, entry.first);
            }
            if (m_num_explicitly_covered >
                0.33 * m_subproblem.uncovered_universe.size())
            {
                p_make_all_explicit();
            } else if (m_selected_uncovered.size() <
                       0.025 * m_subproblem.uncovered_universe.size())
            {
                p_sample_more_explicit();
            }
            if (m_infeasible_by_construction) {
                p_store_event("INCREMENTAL_SAT_LNS_SOLUTION_WAS_OPTIMAL",
                              event_data, "universe_size", "num_removed",
                              "mes_size", "num_explicitly_covered");
                return false;
            }
        }
    }

    LNSSubproblem move_out_subproblem() noexcept {
        return std::move(m_subproblem);
    }

    const std::vector<DynamicBitset>& get_solution() const {
        return m_current_solution;
    }

    void abort() {
        m_abort.store(true);
        m_solver.terminate();
    }

  private:
    void p_prepare_vertex_set() {
        if (m_subproblem.uncovered_universe.empty()) {
            throw std::logic_error(
                "Improvement algorithm run without uncovered vertices!");
        }
        auto canonical = [](Vertex v) -> Vertex {
            return {std::min(v.first, v.second), std::max(v.first, v.second)};
        };
        std::transform(m_subproblem.uncovered_universe.begin(),
                       m_subproblem.uncovered_universe.end(),
                       m_subproblem.uncovered_universe.begin(), canonical);
        std::sort(m_subproblem.uncovered_universe.begin(),
                  m_subproblem.uncovered_universe.end());
        std::transform(m_subproblem.mutually_exclusive_set.begin(),
                       m_subproblem.mutually_exclusive_set.end(),
                       m_subproblem.mutually_exclusive_set.begin(), canonical);
        std::sort(m_subproblem.mutually_exclusive_set.begin(),
                  m_subproblem.mutually_exclusive_set.end());
    }

    void p_make_literal_partners() {
        const auto& vertices = m_subproblem.uncovered_universe;
        const std::size_t nconcrete = m_subproblem.num_concrete;
        const std::size_t nconclit = 2 * nconcrete;
        m_literal_partners_array.reserve(vertices.size());
        m_literal_partners_of.reserve(nconclit);
        SLit last_first = 0;
        std::size_t current_begin = 0, current_size = 0;
        for (Vertex v : vertices) {
            if (v.first != last_first) {
                m_literal_partners_of.emplace_back(current_begin, current_size);
                while (++last_first < v.first) {
                    m_literal_partners_of.emplace_back(current_size,
                                                       current_size);
                }
                current_begin = current_size;
            }
            m_literal_partners_array.push_back(v.second);
            ++current_size;
        }
        m_literal_partners_of.emplace_back(current_begin, current_size);
        while (++last_first < nconclit) {
            m_literal_partners_of.emplace_back(current_size, current_size);
        }
        m_explicitized_literal_partners.assign(m_literal_partners_array.size(),
                                               false);
        m_covered_literal_partners.assign(m_literal_partners_array.size(),
                                          false);
    }

    std::size_t p_find_index(Vertex v) {
        auto p_of = m_literal_partners_of[v.first];
        auto p_of_begin = p_of.first;
        auto p_of_end = p_of.second;
        auto it = m_literal_partners_array.begin() + p_of_begin;
        auto end = m_literal_partners_array.begin() + p_of_end;
        std::size_t result(std::lower_bound(it, end, v.second) -
                           m_literal_partners_array.begin());
        assert(result < p_of_end);
        assert(m_literal_partners_array[result] == v.second);
        return result;
    }

    void p_make_config_vars() {
        const auto allowed = allowed_configurations();
        m_config_vars.reserve(allowed);
        for (Vertex clique_vertex : mes_vertices()) {
            p_create_config_with_clique(clique_vertex);
            std::size_t index = p_find_index(clique_vertex);
            m_explicitized_literal_partners[index].set();
            ++m_num_explicitly_covered;
        }
        while (m_config_vars.size() < allowed) {
            p_create_config_with_clique({NIL, NIL});
        }
    }

    void p_create_config_with_clique(Vertex clique_vertex) {
        const auto nall = m_propagator.db().num_vars();
        m_config_vars.emplace_back();
        auto& config = m_config_vars.back();
        config.reserve(nall);
        m_propagator.reset_to_zero();
        if (clique_vertex.first != clique_vertex.second &&
            push_vertex(m_propagator, clique_vertex) < 0)
        {
            throw std::logic_error("Infeasible interaction in MES!");
        }
        for (SLit var = 0, n = nall; var < n; ++var) {
            SLit pos = lit::positive_lit(var);
            if (m_propagator.is_true(pos)) {
                // true
                config.emplace_back(std::in_place_type<bool>, true);
            } else if (m_propagator.is_false(pos)) {
                // false
                config.emplace_back(std::in_place_type<bool>, false);
            } else {
                // open
                config.emplace_back(std::in_place_type<Lit>,
                                    m_solver.new_var());
            }
        }
        m_propagator.reset_to_zero();
        p_config_insert_binaries();
        p_config_insert_long_clauses();
        throw_if_interrupted();
    }

    LitOrVal p_config_lit_for(const std::vector<LitOrVal>& config,
                              SLit sl) const {
        auto var = lit::var(sl);
        bool is_pos = !lit::negative(sl);
        auto& entry = config[var];
        return std::visit(
            overloaded{
                [&](bool b) -> LitOrVal {
                    return LitOrVal{std::in_place_type<bool>, is_pos ? b : !b};
                },
                [&](Lit l) -> LitOrVal {
                    return LitOrVal{std::in_place_type<Lit>, is_pos ? l : -l};
                }},
            entry);
    }

    void p_config_insert_binaries() {
        const auto& cdb = m_propagator.db();
        const auto& config = m_config_vars.back();
        for (auto bin_clause : cdb.binary_clauses()) {
            SLit l1 = bin_clause.first, l2 = bin_clause.second;
            LitOrVal v1 = p_config_lit_for(config, l1);
            LitOrVal v2 = p_config_lit_for(config, l2);
            if (std::holds_alternative<bool>(v1) ||
                std::holds_alternative<bool>(v2))
            {
                if (std::holds_alternative<bool>(v1) &&
                    std::holds_alternative<bool>(v2))
                {
                    if (!*std::get_if<bool>(&v1) && !*std::get_if<bool>(&v2)) {
                        // both false, construction is infeasible;
                        // this should not actually happen.
                        m_infeasible_by_construction = true;
                        break;
                    }
                }
                // otherwise, one literal is true (either by propagation or
                // directly); need not add clause.
                continue;
            } else {
                m_solver.add_short_clause(*std::get_if<Lit>(&v1),
                                          *std::get_if<Lit>(&v2));
            }
        }
    }

    void p_config_insert_long_clauses() {
        const auto& cdb = m_propagator.db();
        const auto& config = m_config_vars.back();
        for (CRef clause_ref = 1, n = cdb.literal_db_size(); clause_ref < n;
             clause_ref = cdb.next_clause(clause_ref))
        {
            m_buffer.clear();
            bool got_true = false;
            for (SLit l : cdb.lits_of(clause_ref)) {
                LitOrVal v = p_config_lit_for(config, l);
                std::visit(overloaded{[&](bool b) {
                                          if (b) {
                                              got_true = true;
                                          }
                                      },
                                      [&](Lit l) { m_buffer.push_back(l); }},
                           v);
                if (got_true)
                    break;
            }
            if (got_true)
                continue;
            if (m_buffer.empty()) {
                // empty clause, construction is infeasible.
                m_infeasible_by_construction = true;
                break;
            }
            if (m_buffer.size() > 1) {
                m_solver.add_clause(m_buffer.begin(), m_buffer.end());
            }
        }
    }

    void p_make_vertical_clauses() {
        HashSet<SLit> occurring;
        for (Vertex v : m_subproblem.uncovered_universe) {
            occurring.insert(v.first);
            occurring.insert(v.second);
        }
        const auto allowed = allowed_configurations();
        for (SLit l : occurring) {
            m_buffer.clear();
            bool got_true = false;
            for (std::size_t config_index = 0; config_index < allowed;
                 ++config_index)
            {
                LitOrVal lv = p_config_lit_for(m_config_vars[config_index], l);
                std::visit(overloaded{[&](bool b) {
                                          if (b) {
                                              got_true = true;
                                          }
                                      },
                                      [&](Lit l) { m_buffer.push_back(l); }},
                           lv);
                if (got_true)
                    break;
            }
            if (got_true)
                continue;
            if (m_buffer.empty()) {
                m_infeasible_by_construction = true;
                break;
            }
            m_solver.add_clause(m_buffer.begin(), m_buffer.end());
        }
    }

    void p_extract_solution() {
        const Var nall(m_propagator.db().num_vars());
        const auto& model_map = m_solver.get_model();
        if (m_current_solution.empty()) {
            m_current_solution.reserve(allowed_configurations());
            for (std::size_t i = 0; i < allowed_configurations(); ++i) {
                m_current_solution.emplace_back(nall, false);
            }
        }
        for (std::size_t i = 0; i < allowed_configurations(); ++i) {
            auto& current_config = m_config_vars[i];
            auto& current_solution_config = m_current_solution[i];
            for (Var j = 0; j != nall; ++j) {
                std::visit(
                    overloaded{[&](bool b) { current_solution_config[j] = b; },
                               [&](Lit l) {
                                   current_solution_config[j] = model_map[l];
                               }},
                    current_config[j]);
            }
        }
    }

    void p_identify_covered() {
        const auto allowed = allowed_configurations();
        m_covered_literal_partners = m_explicitized_literal_partners;
        for (SLit i = 0, n = m_literal_partners_of.size(); i < n; ++i) {
            for (std::size_t config_index = 0; config_index < allowed;
                 ++config_index)
            {
                const auto& solution = m_current_solution[config_index];
                if (!lit::is_true_in(i, solution)) {
                    continue;
                }
                auto indices = m_literal_partners_of[i];
                for (std::size_t p = indices.first; p < indices.second; ++p) {
                    if (m_covered_literal_partners[p])
                        continue;
                    SLit partner = m_literal_partners_array[p];
                    if (lit::is_true_in(partner, solution)) {
                        m_covered_literal_partners[p].set();
                    }
                }
            }
        }
    }

    void p_make_all_explicit() {
        for (SLit i = 0, n = m_literal_partners_of.size(); i < n; ++i) {
            auto indices = m_literal_partners_of[i];
            for (std::size_t p = indices.first; p < indices.second; ++p) {
                if (m_explicitized_literal_partners[p])
                    continue;
                SLit partner = m_literal_partners_array[p];
                p_make_explicit_possibly_fixed(p, {i, partner});
            }
        }
    }

    void p_sample_more_explicit() {
        std::geometric_distribution<std::size_t> dist(0.025);
        auto& rng = sammy::rng();
        std::size_t skip = dist(rng);
        for (SLit i = 0, n = m_literal_partners_of.size(); i < n; ++i) {
            auto indices = m_literal_partners_of[i];
            for (std::size_t p = indices.first; p < indices.second; ++p) {
                if (skip > 0) {
                    --skip;
                    continue;
                }
                skip = dist(rng);
                if (m_explicitized_literal_partners[p])
                    continue;
                SLit partner = m_literal_partners_array[p];
                p_make_explicit_possibly_fixed(p, {i, partner});
            }
        }
    }

    void p_sample_uncovered() {
        std::size_t goal =
            (std::max)(m_num_explicitly_covered + 1, std::size_t(500));
        m_selected_uncovered.clear();
        for (SLit i = 0, n = m_literal_partners_of.size(); i < n; ++i) {
            std::size_t q = m_literal_partners_of[i].second;
            for (std::size_t p = m_literal_partners_of[i].first; p < q; ++p) {
                if (m_covered_literal_partners[p])
                    continue;
                m_selected_uncovered.push_back(
                    {{i, m_literal_partners_array[p]}, p});
            }
        }
        if (m_selected_uncovered.empty()) {
            return;
        }
        p_store_event("DETECTED_UNCOVERED_INTERACTIONS",
                      {{"explicitly_covered", m_num_explicitly_covered},
                       {"num_uncovered", m_selected_uncovered.size()}},
                      "explicitly_covered", "num_uncovered");
        if (goal >= 0.75 * m_selected_uncovered.size())
            return;
        std::shuffle(m_selected_uncovered.begin(), m_selected_uncovered.end(),
                     sammy::rng());
        m_selected_uncovered.resize(goal);
    }

    std::optional<bool> p_solve_iteration() {
        if (m_abort.load())
            return std::nullopt;
        auto res = m_solver.solve();
        if (m_abort.load())
            return std::nullopt;
        return res;
    }

    /**
     * Called to explicitly require a currently
     * uncovered interaction to be covered.
     */
    void p_make_explicit_unfixed(std::size_t index, Vertex vertex) {
        const auto visit_combination = overloaded{
            [&](bool b1, bool b2) {
                if (!b1 || !b2)
                    return;
                throw std::logic_error(
                    "Asked to explicitize already fixed interaction!");
            },
            [&](bool b1, Lit l2) {
                if (!b1)
                    return;
                m_buffer.push_back(l2);
            },
            [&](Lit l1, bool b2) {
                if (!b2)
                    return;
                m_buffer.push_back(l1);
            },
            [&](Lit l1, Lit l2) {
                Lit coverage_var = m_solver.new_var();
                m_solver.add_short_clause(coverage_var, -l1, -l2);
                m_solver.add_short_clause(-coverage_var, l1);
                m_solver.add_short_clause(-coverage_var, l2);
                m_buffer.push_back(coverage_var);
            }};

        ++m_num_explicitly_covered;
        m_explicitized_literal_partners[index].set();
        m_buffer.clear();
        for (const auto& config : m_config_vars) {
            LitOrVal lv1 = p_config_lit_for(config, vertex.first);
            LitOrVal lv2 = p_config_lit_for(config, vertex.second);
            std::visit(visit_combination, lv1, lv2);
        }
        if (m_buffer.empty()) {
            p_store_event("EMPTY_COVERAGE_CLAUSE");
            m_infeasible_by_construction = true;
        } else {
            m_solver.add_clause(m_buffer.begin(), m_buffer.end());
            m_buffer.clear();
        }
    }

    void p_make_explicit_possibly_fixed(std::size_t index, Vertex vertex) {
        for (const auto& config : m_config_vars) {
            LitOrVal lv1 = p_config_lit_for(config, vertex.first);
            LitOrVal lv2 = p_config_lit_for(config, vertex.second);
            if (std::holds_alternative<bool>(lv1) &&
                std::holds_alternative<bool>(lv2))
            {
                bool b1 = *std::get_if<bool>(&lv1);
                bool b2 = *std::get_if<bool>(&lv2);
                if (b1 && b2) {
                    ++m_num_explicitly_covered;
                    m_explicitized_literal_partners[index].set();
                    return;
                }
            }
        }
        p_make_explicit_unfixed(index, vertex);
    }

    template <typename... Args>
    void p_store_event(const char* name, OutputObject data, Args&&... keys) {
        if (!m_recorder)
            return;
        data["worker_id"] = m_worker_id;
        m_recorder->store_event(name, std::move(data),
                                std::forward<Args>(keys)...);
    }

    void p_store_event(const char* name) {
        p_store_event(name, {{"worker_id", m_worker_id}});
    }
};

} // namespace sammy

#endif
==> ./parallel_bit_filter.h <==
#ifndef SAMMY_PARALLEL_BIT_FILTER_H_INCLUDED_
#define SAMMY_PARALLEL_BIT_FILTER_H_INCLUDED_

#include "dynamic_bitset.h"
#include "thread_group.h"
#include <vector>

namespace sammy {

/**
 * Buffer for parallel bitset operations.
 */
class BitsetOperationsBuffer {
  public:
    explicit BitsetOperationsBuffer(ThreadGroup<void>* pool)
        : thread_pool(pool), bitsets(pool->num_threads() + 1, Bitset{}),
          used_threads(pool->num_threads() + 1, 0) {}

    void prepare() { used_threads.assign(thread_pool->num_threads() + 1, 0); }

    template <typename Combiner>
    Bitset& initialize(std::size_t index, const Bitset& from_bs,
                       Combiner&& combiner) {
        Bitset& s = bitsets[index];
        if (used_threads[index]) {
            // re-initialization (if one thread got two tasks)
            combiner(s, from_bs);
        } else {
            // first time this thread is used this run
            used_threads[index] = 1;
            s = from_bs;
        }
        return s;
    }

    template <typename Reducer>
    void reduce(Bitset& initial_and_result, Reducer&& reducer) {
        for (std::size_t i = 0, nt = bitsets.size(); i < nt; ++i) {
            if (!used_threads[i])
                continue;
            std::forward<Reducer>(reducer)(initial_and_result, bitsets[i]);
        }
    }

    ThreadGroup<void>& thread_group() { return *thread_pool; }

  private:
    ThreadGroup<void>* thread_pool;
    std::vector<Bitset> bitsets;
    std::vector<std::size_t> used_threads;
};

/**
 * Given an existing thread group, a range of bitsets y_1, ..., y_n,
 * and an initial set x, compute x - y_1 - ... - y_n, possibly in parallel
 * using the thread group of op_buffer, and write the result to
 * initial_and_result.
 */
template <typename BitsetsIterator>
inline void bitwise_filter(BitsetOperationsBuffer& op_buffer,
                           DynamicBitset& initial_and_result,
                           BitsetsIterator begin, BitsetsIterator end) {
    auto& tgroup = op_buffer.thread_group();
    auto nt = tgroup.num_threads() + 1;
    std::size_t num_sets = std::distance(begin, end);
    std::size_t total_size = initial_and_result.size() * num_sets / 8;
    if (nt <= 1 || num_sets < 4 * nt || total_size <= 1024 * 1024) {
        std::for_each(begin, end, [&](const DynamicBitset& bs) {
            initial_and_result -= bs;
        });
    } else {
        op_buffer.prepare();
        tgroup.context_function_parallel_foreach_iterator(
            begin, end,
            [&](std::size_t thread_index,
                const BitsetsIterator& first_iter) -> Bitset& {
                return op_buffer.initialize(
                    thread_index, *first_iter,
                    [](Bitset& bout, const Bitset& bin) { bout |= bin; });
            },
            [&](Bitset& output, std::size_t, const BitsetsIterator& iter) {
                output |= *iter;
            },
            [](Bitset&, std::size_t) {});
        op_buffer.reduce(
            initial_and_result,
            [](Bitset& in_out, const Bitset& in) { in_out -= in; });
    }
}

/**
 * Given an existing thread group, a range of bitsets y_1, ..., y_n,
 * and an initial set x, compute x - y_1 - ... - y_n, possibly in parallel
 * using the thread group tgroup, and write the result to initial_and_result.
 */
template <typename BitsetsIterator>
inline void bitwise_and(BitsetOperationsBuffer& op_buffer,
                        DynamicBitset& initial_and_result,
                        BitsetsIterator begin, BitsetsIterator end) {
    auto& tgroup = op_buffer.thread_group();
    auto nt = tgroup.num_threads() + 1;
    std::size_t num_sets = std::distance(begin, end);
    std::size_t total_size = initial_and_result.size() * num_sets / 8;
    if (nt <= 1 || num_sets < 4 * nt || total_size <= 1024 * 1024) {
        std::for_each(begin, end, [&](const DynamicBitset& bs) {
            initial_and_result &= bs;
        });
    } else {
        op_buffer.prepare();
        tgroup.context_function_parallel_foreach_iterator(
            begin, end,
            [&](std::size_t thread_index,
                const BitsetsIterator& first_iter) -> Bitset& {
                return op_buffer.initialize(
                    thread_index, *first_iter,
                    [](Bitset& bout, const Bitset& bin) { bout &= bin; });
            },
            [&](Bitset& output, std::size_t, const BitsetsIterator& iter) {
                output &= *iter;
            },
            [](Bitset&, std::size_t) {});
        op_buffer.reduce(
            initial_and_result,
            [](Bitset& in_out, const Bitset& in) { in_out &= in; });
    }
}

} // namespace sammy

#endif
==> ./initial_phase_result.h <==
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
==> ./clique_storage.h <==
#ifndef SAMMY_CLIQUE_STORAGE_H_INCLUDED_
#define SAMMY_CLIQUE_STORAGE_H_INCLUDED_

#include "algorithm_ex.h"
#include "literals.h"
#include "range.h"
#include <boost/iterator/iterator_facade.hpp>
#include <memory>
#include <mutex>

namespace sammy {

class CliqueStorage {
  public:
    /**
     * Create a clique storage that can take
     * cliques with a combined number of up to vertex_limit vertices.
     */
    explicit CliqueStorage(std::size_t vertex_limit = 10'000'000)
        : m_vertex_list(new Vertex[p_make_limit(vertex_limit)]),
          m_index_list(new std::size_t[vertex_limit / 2 + 1]),
          m_last_used(new std::size_t[vertex_limit / 2]),
          m_vertex_limit(vertex_limit) {
        m_index_list[0] = 0;
    }

    /**
     * A read-only view of the contents of the clique storage
     * at the point the view was acquired.
     * This allows us to iterate through the cliques without
     * holding a lock while other threads may be adding new cliques
     * (or even purging/moving around cliques when the storage is full).
     * Behaves like a container of CliqueView (IteratorRange<const Vertex*>).
     * The vertices are stored in a contiguous array.
     */
    class StorageView {
        std::shared_ptr<Vertex[]> vertex_list;
        std::shared_ptr<std::size_t[]> clique_index_list;
        std::shared_ptr<std::size_t[]> last_used;
        std::size_t clique_count;

        friend class CliqueStorage;

        StorageView(std::shared_ptr<Vertex[]> v,
                    std::shared_ptr<std::size_t[]> i,
                    std::shared_ptr<std::size_t[]> u, std::size_t cc) noexcept
            : vertex_list(std::move(v)), clique_index_list(std::move(i)),
              last_used(std::move(u)), clique_count(cc) {}

      public:
        using CliqueView = IteratorRange<const Vertex*>;

        CliqueView operator[](std::size_t index) const noexcept {
            assert(index < clique_count);
            std::size_t ibeg = clique_index_list[index];
            std::size_t iend = clique_index_list[index + 1];
            return CliqueView{vertex_list.get() + ibeg,
                              vertex_list.get() + iend};
        }

        bool empty() const noexcept { return clique_count == 0; }

        std::size_t size() const noexcept { return clique_count; }

        class Iterator
            : public boost::iterator_facade<Iterator, CliqueView,
                                            std::random_access_iterator_tag,
                                            CliqueView, std::ptrdiff_t> {
          public:
            Iterator() = default;

            explicit Iterator(const Vertex* vbase,
                              const std::size_t* index_ptr) noexcept
                : vbase(vbase), index_ptr(index_ptr) {}

          private:
            friend boost::iterator_core_access;
            friend class CliqueStorage;

            CliqueView dereference() const noexcept {
                std::size_t ibeg = index_ptr[0];
                std::size_t iend = index_ptr[1];
                return CliqueView{vbase + ibeg, vbase + iend};
            }

            bool equal(const Iterator& other) const noexcept {
                return index_ptr == other.index_ptr;
            }

            void increment() noexcept { ++index_ptr; }

            void decrement() noexcept { --index_ptr; }

            void advance(std::ptrdiff_t diff) noexcept { index_ptr += diff; }

            std::ptrdiff_t distance_to(const Iterator& other) const noexcept {
                return other.index_ptr - index_ptr;
            }

            const Vertex* vbase;
            const std::size_t* index_ptr;
        };

        Iterator begin() const noexcept {
            return Iterator{vertex_list.get(), clique_index_list.get()};
        }

        Iterator end() const noexcept {
            return Iterator{vertex_list.get(),
                            clique_index_list.get() + clique_count};
        }
    };

    /**
     * Add a new clique to the CliqueStorage.
     * Two things can happen:
     *  * The CliqueStorage has space for the new clique.
     *    In that case, the clique is simply added.
     *  * The CliqueStorage is full.
     *    In that case, the clique storage is reduced by
     *    throwing out a quarter of all cliques
     *    until there is space for the new clique,
     *    using a LRU pattern to decide which cliques are purged.
     */
    template <typename Iterator>
    void push_clique(Iterator begin, Iterator end) {
        std::size_t cs = std::distance(begin, end);
        if (cs <= 1)
            return;

        std::unique_lock l{m_lock};
        std::size_t* il = m_index_list.get();
        Vertex* vl = m_vertex_list.get();
        std::size_t current_size = il[m_clique_count];
        std::size_t new_end = current_size + cs;
        if (new_end >= m_vertex_limit) {
            if (2 * cs > m_vertex_limit) {
                m_vertex_limit = 2 * cs;
            }
            p_make_space(cs);
        }
        il[++m_clique_count] = new_end;
        std::copy(begin, end, &vl[current_size]);
        m_last_used[m_clique_count - 1] = m_timestamp++;
    }

    /**
     * Notify the CliqueStorage that a certain clique was used.
     */
    void used_clique(StorageView& view, std::size_t index) {
        std::unique_lock l{m_lock};
        view.last_used[index] = m_timestamp++;
    }

    /**
     * Notify the CliqueStorage that a certain clique was used.
     */
    void used_clique(StorageView& view, StorageView::Iterator iterator) {
        auto index =
            std::size_t(iterator.index_ptr - view.clique_index_list.get());
        used_clique(view, index);
    }

    /**
     * Obtain a view of the clique list.
     * The part of the clique list that is returned is immutable
     * but will stay valid and does not require holding a lock
     * while walking the list of cliques (which may take long).
     */
    StorageView obtain_view() const {
        std::unique_lock l{m_lock};
        return StorageView{m_vertex_list, m_index_list, m_last_used,
                           m_clique_count};
    }

  private:
    static std::size_t& p_make_limit(std::size_t& x) noexcept {
        if (x % 2)
            x += 1;
        if (x < 1000)
            x = 1000;
        return x;
    }

    std::size_t p_clique_size(std::size_t index) {
        return m_index_list[index + 1] - m_index_list[index];
    }

    void p_make_space(std::size_t required_space) {
        auto compare_lru = [&](std::size_t i1, std::size_t i2) {
            return m_last_used[i1] < m_last_used[i2];
        };
        std::vector<std::size_t> indices = vector(range(m_clique_count));
        std::size_t remove_cliques =
            (std::min)(1 + m_clique_count / 4, m_clique_count - 1);
        auto removed_end = indices.begin() + remove_cliques;
        for (;;) {
            std::nth_element(indices.begin(), removed_end, indices.end(),
                             compare_lru);
            std::size_t new_free_space =
                m_vertex_limit - m_index_list[m_clique_count];
            std::for_each(indices.begin(), removed_end, [&](std::size_t i) {
                new_free_space += p_clique_size(i);
            });
            if (new_free_space >= required_space)
                break;
            remove_cliques = (std::min)(2 * remove_cliques, m_clique_count - 1);
            removed_end = indices.begin() + remove_cliques;
        }
        std::shared_ptr<Vertex[]> new_vertex_list(new Vertex[m_vertex_limit]);
        std::shared_ptr<std::size_t[]> new_index_list(
            new std::size_t[m_vertex_limit / 2 + 1]);
        std::shared_ptr<std::size_t[]> new_last_used(
            new std::size_t[m_vertex_limit / 2]);
        Vertex* new_out = new_vertex_list.get();
        std::size_t out_i = 0;
        new_index_list[0] = 0;
        auto copy_clique = [&](std::size_t old_index) {
            std::size_t old_begin = m_index_list[old_index];
            std::size_t old_end = m_index_list[old_index + 1];
            new_out = std::copy(&m_vertex_list[old_begin],
                                &m_vertex_list[old_end], new_out);
            new_index_list[out_i + 1] =
                std::size_t(new_out - new_vertex_list.get());
            new_last_used[out_i] = m_last_used[old_index];
            ++out_i;
        };
        std::for_each(removed_end, indices.end(), copy_clique);
        m_clique_count = out_i;
        m_vertex_list = std::move(new_vertex_list);
        m_index_list = std::move(new_index_list);
        m_last_used = std::move(new_last_used);
    }

    mutable std::mutex m_lock;
    std::shared_ptr<Vertex[]> m_vertex_list;
    std::shared_ptr<std::size_t[]> m_index_list;
    std::shared_ptr<std::size_t[]> m_last_used;
    std::size_t m_clique_count = 0;
    std::size_t m_vertex_limit;
    std::size_t m_timestamp = 0;
};

} // namespace sammy

#endif
==> ./external_sat_solver.h <==
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
                assert(index < raw.size());
                return !m_raw[index];
            } else {
                auto index = l - 1;
                assert(index < raw.size());
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
==> ./instance.h <==
#ifndef SAMMY_INSTANCE_H_INCLUDED_
#define SAMMY_INSTANCE_H_INCLUDED_

namespace sammy {}

#endif
==> ./eliminate_subsumed.h <==
#ifndef SAMMY_ELIMINATE_SUBSUMED_H_INCLUDED_
#define SAMMY_ELIMINATE_SUBSUMED_H_INCLUDED_

#include "error.h"
#include "literals.h"
#include "simplification_stats.h"
#include "stamp_set.h"

namespace sammy {

template <typename ClauseType> class SubsumptionChecker {
  public:
    SubsumptionChecker(std::vector<ClauseType>& clauses, Var n_all)
        : m_nv(n_all), m_nl(2 * n_all), m_clauses(clauses), m_in_clause(m_nl),
          m_watching_clauses(m_nl), stats(&stats_buffer) {
        p_init_watches();
    }

    void set_stats(SimplificationStats* stats) noexcept { this->stats = stats; }

    void remove_subsumed() {
        for (CRef c = 0, n = m_clauses.size(); c < n; ++c) {
            p_empty_if_subsumed(c);
        }
        auto deleted_begin =
            std::remove_if(m_clauses.begin(), m_clauses.end(),
                           [](const ClauseType& cl) { return cl.empty(); });
        stats->clauses_subsumed += std::size_t(m_clauses.end() - deleted_begin);
        m_clauses.erase(deleted_begin, m_clauses.end());
    }

  private:
    bool p_walk_watch_list(CRef index, Lit l) {
        auto& watch_list = m_watching_clauses[l];
        auto end = watch_list.end();
        auto out = watch_list.begin();
        bool subsumed = false;
        for (auto in = watch_list.begin(); in != end; ++in) {
            CRef cother = *in;
            // we cannot subsume ourself. stay in the watch list.
            if (cother == index) {
                *out++ = cother;
                continue;
            }
            const ClauseType& other_lits = m_clauses[cother];
            // subsumed clauses do not participate in subsumption anymore;
            // they are dropped from watch lists without replacement when we
            // encounter them here.
            if (other_lits.empty()) {
                continue;
            }
            // find replacement watch (must not be in the current clause).
            auto replacement =
                std::find_if(other_lits.begin(), other_lits.end(),
                             [&](Lit l) { return !m_in_clause.count(l); });
            if (replacement == other_lits.end()) {
                // cother subsumes us.
                subsumed = true;
                // copy remaining watching clauses.
                out = std::copy(in, end, out);
                break;
            } else {
                // cother does not subsume us.
                m_watching_clauses[*replacement].push_back(cother);
            }
        }
        // trim watch list
        watch_list.erase(out, end);
        return subsumed;
    }

    void p_empty_if_subsumed(CRef index) {
        ClauseType& clause = m_clauses[index];
        m_in_clause.assign(clause.begin(), clause.end());
        for (Lit l : clause) {
            if (p_walk_watch_list(index, l)) {
                clause.clear();
                return;
            }
        }
    }

    void p_init_watches() {
        for (std::size_t ci = 0, cn = m_clauses.size(); ci < cn; ++ci) {
            const auto& cl = m_clauses[ci];
            m_watching_clauses[cl[0]].push_back(CRef(ci));
        }
    }

    Var m_nv;
    Lit m_nl;
    std::vector<ClauseType>& m_clauses;
    StampSet<Lit, std::uint16_t> m_in_clause;
    std::vector<std::vector<CRef>> m_watching_clauses;
    SimplificationStats stats_buffer;
    SimplificationStats* stats;
};

template <typename ClauseType>
inline void eliminate_subsumed(std::vector<ClauseType>& clauses, Var n_all,
                               SimplificationStats* stats = nullptr) {
    SubsumptionChecker<ClauseType> subsumption_checker{clauses, n_all};
    if (stats)
        subsumption_checker.set_stats(stats);
    subsumption_checker.remove_subsumed();
}

} // namespace sammy

#endif
==> ./best_k.h <==
#ifndef SAMMY_BEST_K_INCLUDED_
#define SAMMY_BEST_K_INCLUDED_

#include <algorithm>
#include <utility>
#include <vector>

namespace sammy {

/**
 * A container for keeping track of the
 * best k elements of a set of elements.
 * Which elements are best is determined
 * using Compare::operator(), defaulting
 * to < on the elements (in which case,
 * the smallest elements are best).
 */
template <typename T, typename Compare = std::less<T>> class BestK {
  public:
    explicit BestK(std::size_t k, Compare&& compare = {})
        : m_k(k), m_compare(std::forward<Compare>(compare)) {
        m_elements.reserve(m_k);
    }

    void push(T&& element) { p_push(std::move(element)); }

    void push(const T& element) { p_push(element); }

    template <typename ComparableToElement>
    bool would_push(ComparableToElement&& score) const noexcept {
        if (m_elements.size() < m_k) {
            return true;
        }
        return !m_compare(m_elements[0], score);
    }

    template <typename... Args> void emplace(Args&&... args) {
        T new_element(std::forward<Args>(args)...);
        push(std::move(new_element));
    }

    void clear() noexcept { m_elements.clear(); }

    const std::vector<T>& elements() const noexcept { return m_elements; }

  private:
    template <typename TT> void p_push(TT&& element) {
        if (m_elements.size() < m_k) {
            m_elements.emplace_back(std::forward<TT>(element));
            if (m_elements.size() == m_k) {
                std::make_heap(m_elements.begin(), m_elements.end(), m_compare);
            }
            return;
        }
        if (!m_compare(element, m_elements[0])) {
            // we don't have to do anything for
            // elements that aren't better than
            // the worst of the best k elements so far
            return;
        }
        m_elements[0] = std::forward<TT>(element);
        p_sift_down();
    }

    void p_sift_down() {
        const std::size_t s = m_k;
        std::size_t ei = 0;
        while (2 * ei + 2 < s) {
            std::size_t c = 2 * ei + 1;
            if (m_compare(m_elements[c], m_elements[c + 1])) {
                // the left child is better,
                // so if any child, the right one moves up
                ++c;
            }
            // if the current element isn't better than its
            // worst child, we are done.
            if (!m_compare(m_elements[ei], m_elements[c]))
                return;
            std::swap(m_elements[ei], m_elements[c]);
            ei = c;
        }
        const std::size_t c = 2 * ei + 1;
        if (c < s) {
            // only one child; if we are better, swap
            if (m_compare(m_elements[ei], m_elements[c])) {
                std::swap(m_elements[ei], m_elements[c]);
            }
        }
    }

    std::size_t m_k;
    std::vector<T> m_elements;
    Compare m_compare;
};

} // namespace sammy

#endif
==> ./clique_or_indset_builder.h <==
#ifndef SAMMY_CLIQUE_OR_INDSET_BUILDER_H_INCLUDED_
#define SAMMY_CLIQUE_OR_INDSET_BUILDER_H_INCLUDED_

#include "dynamic_bitset.h"
#include "literals.h"
#include "parallel_bit_filter.h"
#include <boost/iterator/transform_iterator.hpp>

namespace sammy {

template <typename GraphType, bool BuildingClique>
class CliqueOrIndependentSetBuilder {
  public:
    explicit CliqueOrIndependentSetBuilder(
        GraphType* graph, BitsetOperationsBuffer* parallel_bits)
        : graph(graph), parallel_bits(parallel_bits),
          possible_vertices(graph->n(), true) {}

    void reset_vertices() {
        auto g_n = graph->n();
        if (possible_vertices.size() != g_n) {
            possible_vertices.resize(g_n);
        }
        possible_vertices.set();
        clique_vertices.clear();
    }

    void reset_vertices(const std::vector<Vertex>& initial) {
        reset_vertices();
        if (initial.empty()) {
            return;
        }
        auto transform_to_index = [this](Vertex v) {
            return graph->vertex_index(v);
        };
        std::transform(initial.begin(), initial.end(),
                       std::back_inserter(clique_vertices), transform_to_index);
        for (std::size_t ci : clique_vertices) {
            possible_vertices[ci] = false;
        }
        auto index_to_bitset = [&](std::size_t i) -> const DynamicBitset& {
            return graph->matrix_row(i);
        };
        auto bitset_begin =
            boost::make_transform_iterator<decltype(index_to_bitset)>(
                clique_vertices.begin());
        auto bitset_end =
            boost::make_transform_iterator<decltype(index_to_bitset)>(
                clique_vertices.end());
        if constexpr (BuildingClique) {
            sammy::bitwise_and(*parallel_bits, possible_vertices, bitset_begin,
                               bitset_end);
        } else {
            sammy::bitwise_filter(*parallel_bits, possible_vertices,
                                  bitset_begin, bitset_end);
        }
    }

    template <typename RNG> void randomly_make_maximal(RNG& rng) {
        for (;;) {
            bool success = false;
            std::uniform_int_distribution<std::size_t> indices(0,
                                                               graph->n() - 1);
            for (std::size_t trial_count = 0; trial_count < 100; ++trial_count)
            {
                auto index = indices(rng);
                if (possible_vertices[index]) {
                    add_vertex(index);
                    success = true;
                    break;
                }
            }
            if (!success) {
                possible_buffer.clear();
                std::copy(possible_vertices.ones_begin(),
                          possible_vertices.ones_end(),
                          std::back_inserter(possible_buffer));
                while (!possible_buffer.empty()) {
                    std::uniform_int_distribution<std::size_t> iindex_dist(
                        0, possible_buffer.size() - 1);
                    auto iindex = iindex_dist(rng);
                    auto vindex = possible_buffer[iindex];
                    add_vertex(vindex);
                    possible_buffer.erase(
                        std::remove_if(possible_buffer.begin(),
                                       possible_buffer.end(),
                                       [&](std::size_t v) {
                                           return !possible_vertices[v];
                                       }),
                        possible_buffer.end());
                }
                return;
            }
        }
    }

    void add_vertex(std::size_t index) {
        clique_vertices.push_back(index);
        possible_vertices[index] = false;
        if (BuildingClique) {
            possible_vertices &= graph->matrix_row(index);
        } else {
            possible_vertices -= graph->matrix_row(index);
        }
    }

    void add_vertex(Vertex v) { add_vertex(graph->vertex_index(v)); }

    std::size_t size() const noexcept { return clique_vertices.size(); }

    const std::vector<std::size_t>& get_indices() const noexcept {
        return clique_vertices;
    }

    std::vector<Vertex> get_vertices() const {
        std::vector<Vertex> result;
        result.reserve(clique_vertices.size());
        std::transform(clique_vertices.begin(), clique_vertices.end(),
                       std::back_inserter(result),
                       [&](std::size_t i) { return graph->vertex(i); });
        return result;
    }

    template <typename VertexIterator>
    void greedily_extend_on(VertexIterator begin, VertexIterator end) {
        for (auto v : IteratorRange(begin, end)) {
            auto index = graph->vertex_index(v);
            if (possible_vertices[index]) {
                add_vertex(index);
            }
        }
    }

    template <typename VertexIndexIterator>
    void greedily_extend_on_indices(VertexIndexIterator begin,
                                    VertexIndexIterator end) {
        for (auto index : IteratorRange(begin, end)) {
            if (possible_vertices[index]) {
                add_vertex(index);
            }
        }
    }

    template <typename WeightedIndexIterator>
    double greedily_extend_weighted(WeightedIndexIterator begin,
                                    WeightedIndexIterator end) {
        double total = 0.0;
        for (auto [index, weight] : IteratorRange(begin, end)) {
            if (possible_vertices[index]) {
                add_vertex(index);
                total += weight;
            }
        }
        return total;
    }

  private:
    GraphType* graph;
    BitsetOperationsBuffer* parallel_bits;
    std::vector<std::size_t> clique_vertices;
    DynamicBitset possible_vertices;
    std::vector<std::size_t> possible_buffer;
};

template <typename GraphType>
using CliqueBuilder = CliqueOrIndependentSetBuilder<GraphType, true>;
template <typename GraphType>
using IndependentSetBuilder = CliqueOrIndependentSetBuilder<GraphType, false>;

} // namespace sammy

#endif
