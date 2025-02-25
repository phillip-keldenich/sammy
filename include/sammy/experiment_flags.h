#ifndef SAMMY_EXPERIMENT_FLAGS_H_INCLUDED_
#define SAMMY_EXPERIMENT_FLAGS_H_INCLUDED_

#include "run_initial.h"
#include "output.h"
#include "lower_bound_mip_settings.h"
#include "basic_cuda_iteration.h"
#include <cstdlib>
#include <boost/program_options.hpp>

namespace sammy {

template<typename ValueType>
inline auto value_with_default(ValueType& vref)
{
    ValueType default_value{vref};
    return boost::program_options::value(&vref)->default_value(default_value);
}

template<typename ValueType>
inline auto required_value(ValueType& vref)
{
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
        if(!output_file.empty() && !recognized_output_extension(output_file)) {
            return "Error: Unrecognized output file extension (use .json or .jsonl)!";
        }
        if(!std::filesystem::is_regular_file(input_file)) {
            return "Error: The given input file '" + input_file + "' is not an existing regular file!";
        }
        if(cuda_mode == "off") { set_cuda_mode(CUDA_USAGE_DISABLED); }
        else if(cuda_mode == "forced") { set_cuda_mode(CUDA_USAGE_FORCED); }
        else { set_cuda_mode(CUDA_USAGE_IF_AVAILABLE); }
        return std::nullopt;
    }
};

inline
void add_gurobi_experiment_flags(boost::program_options::options_description& description,
                                 ExperimentFlagsConfig& config)
{
    auto& lcfg = config.lb_mip_config;
    description.add_options()
        ("lb-exactly-solve-full-mip",
         bool_switch(lcfg.solve_exactly_full_mip),
         "Try to solve the clique problem on the universe conflict graph exactly, "
         "using a MIP with additional cutting planes applied at the root node.")
        ("lb-exactly-branch-and-price",
         bool_switch(lcfg.solve_exactly_cut_and_price),
         "Try to solve the clique problem on the universe conflict graph exactly, "
         "using not a MIP but cutting planes and pricing for upper bounds and "
         "'accidental' integrality, rounding and other heuristics for lower bounds.")
        ("lb-mip-stop-cuts-iteration-window",
         value_with_default(lcfg.stop_cuts_iteration_window),
         "Specify the length (in iterations) of the window to consider when deciding "
         "whether to stop adding cutting planes.")
        ("lb-mip-stop-cuts-below-relative-improvement",
         value_with_default(lcfg.stop_cuts_below_relative_improvement),
         "Specify the required relative improvement of the upper bound "
         "needed accross an iteration window below which generating cutting planes "
         "shall be stopped.")
        ("separation-rounds-before-pricing", value_with_default(lcfg.separation_rounds_before_pricing),
         "Specify the number of separation rounds to do between pricing rounds "
         "when faced with a non-optimal-for-subgraph LP relaxation.")
        ("lb-mip-timeout",
         value_with_default(lcfg.total_mip_timeout),
         "Specify the maximum amount of time to spend on lower bound MIP approach.")
        ("lb-mip-gurobi-quiet", bool_switch(lcfg.quiet_gurobi),
         "Tell Gurobi to be quiet while solving the lower bound MIPs/LPs.");
}

inline
void add_sat_coloring_flags(boost::program_options::options_description& description,
                                 ExperimentFlagsConfig& config)
{
    auto& scfg = config.sat_coloring_config;
    description.add_options()
        ("initial-sat-coloring", bool_switch(scfg.initial_sat_coloring),
         "Use SAT coloring for some time to initialize upper and lower bounds.")
        ("initial-sat-coloring-time-limit", value_with_default(scfg.initial_sat_coloring_timeout),
         "Time limit for initial SAT coloring.");
}

inline
void add_experiment_flags(boost::program_options::options_description& description,
                          ExperimentFlagsConfig& config, bool enable_gurobi_clique_flags = false,
                          bool enable_sat_coloring_flags = false) 
{
    auto& icfg = config.initial_heuristic_config;
    description.add_options()
        ("help,h", "produce help output")
        ("dont-simplify", bool_switch(config.dont_simplify), "do not apply simplification")
        ("initial-time-limit", value_with_default(icfg.max_time), 
         "Time limit for the initial heuristic.")
        ("cuda-mode", value_with_default(config.cuda_mode), 
         "select CUDA mode: either 'off', 'force', or 'auto'. "
         "If 'off', no CUDA will be used. If 'force', CUDA will be used; "
         "an error is raised if this cannot be done. "
         "If 'auto' (or any other value), CUDA will be used if possible, falling back to CPU otherwise.")
        ("initial-iteration-goal", value_with_default(icfg.goal_iterations),
         "Iterations of the initial heuristic to run at least (unless the time limit is exceeded).")
        ("initial-goal-time", value_with_default(icfg.min_time),
         "Time (in seconds) to spend on the initial heuristic at least "
         "(unless the time limit is exceeded or optimality is proved).")
        ("initial-random-clique-restarts",
         value_with_default(icfg.random_clique_restarts_per_iteration),
         "Number of repeats (per thread) to run the randomized greedy clique "
         "algorithm for per iteration of the initial heuristic.")
        ("initial-heuristic-progress", bool_switch(config.initial_progress),
         "Print progress information during the initial heuristic")
        ("print-events", bool_switch(config.print_events),
         "Print time-stamped event information")
        ("output-file,o", value_with_default(config.output_file));
    if(enable_gurobi_clique_flags) {
        add_gurobi_experiment_flags(description, config);
    }
    if(enable_sat_coloring_flags) {
        add_sat_coloring_flags(description, config);
    }
}

inline
void add_input_file_flag(boost::program_options::options_description& description,
                         boost::program_options::positional_options_description& positional,
                         ExperimentFlagsConfig& config)
{
    description.add_options()
        ("input-file", boost::program_options::value<std::string>(&config.input_file)->required(),
         "The input file to read.");
    positional.add("input-file", 1);
}

inline
bool parse_cmdline(const boost::program_options::options_description& all_options,
                   const boost::program_options::positional_options_description& positional,
                   boost::program_options::variables_map& vmap, int argc, char** argv)
{
    auto cmdline_parser = boost::program_options::command_line_parser(argc, const_cast<const char* const*>(argv));
    boost::program_options::store(cmdline_parser.options(all_options).positional(positional).run(), vmap);
    vmap.notify();
    return !vmap.count("help");
}

inline void get_config_or_exit(int argc, char** argv, ExperimentFlagsConfig& defaults_and_out,
                               const boost::program_options::options_description& extra_options = {},
                               bool enable_gurobi_clique_flags = false,
                               bool enable_sat_coloring_flags = false)
{
    boost::program_options::options_description visible("Allowed options");
    boost::program_options::options_description hidden("Hidden options");
    boost::program_options::positional_options_description positional;
    visible.add(extra_options);
    add_experiment_flags(visible, defaults_and_out, enable_gurobi_clique_flags, enable_sat_coloring_flags);
    add_input_file_flag(hidden, positional, defaults_and_out);
    boost::program_options::options_description all;
    all.add(visible).add(hidden);
    boost::program_options::variables_map vmap;
    try {
        if(!parse_cmdline(all, positional, vmap, argc, argv)) {
            std::cout << visible << std::endl;
            std::exit(0);
        }
    } catch(const boost::program_options::error& err) {
        int ecode = 1;
        if(!vmap.count("help")) {
            std::cerr << "Invalid command line options!" << std::endl;
            std::cerr << err.what() << std::endl;
            ecode = 0;
        }
        std::cerr << visible << std::endl;
        std::exit(ecode);
    }
    defaults_and_out.process_parsed_args();
}

inline
ExperimentFlagsConfig get_config_or_exit(int argc, char** argv,
                                         const boost::program_options::options_description& extra_options = {},
                                         bool enable_gurobi_clique_flags = false,
                                         bool enable_sat_coloring_flags = false)
{
    ExperimentFlagsConfig config;
    get_config_or_exit(argc, argv, config, extra_options, enable_gurobi_clique_flags, enable_sat_coloring_flags);
    return config;
}

}

#endif
