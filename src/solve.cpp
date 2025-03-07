#include <sammy/barrage.h>
#include <sammy/barrage_worker_cnp.h>
#include <sammy/barrage_worker_exact.h>
#include <sammy/barrage_worker_lns.h>
#include <sammy/cmsat5_solver.h>
#include <sammy/experiment_flags.h>
#include <sammy/external_sat_solver.h>
#include <sammy/incremental_sat_lns.h>
#include <sammy/kissat_solver.h>
#include <sammy/lingeling_solver.h>
#include <sammy/memory_usage.h>
#include <sammy/output.h>
#include <sammy/partial_solution.h>
#include <sammy/run_initial.h>
#include <sammy/sat_lns.h>
#include <sammy/subproblem_solver_with_mes.h>
#include <sammy/thread_clauses.h>

using namespace sammy;
namespace po = boost::program_options;

void run_simplification(EventRecorder& rec, const InputData& original,
                        std::optional<SimplifyDatastructure>& simplifier,
                        std::optional<SimplifiedInstance>& simplified) {
    rec.store_event("BEGIN_SIMPLIFICATION");
    SimplificationStats simp_stats;
    simp_stats.capture_before(original.formula, original.num_concrete);
    simplifier.emplace(original.formula, original.num_concrete);
    simplified = run_simplifier(*simplifier, simp_stats);
    simp_stats.capture_after(simplified->formula, simplified->num_concrete);
    OutputObject simp_data;
    export_simplification(simp_data, *simplifier, *simplified, &simp_stats);
    rec.store_event("DONE_SIMPLIFICATION", simp_data);
}

InitialPhaseResult run_initial_phase(EventRecorder& rec,
                                     const RunInitialConfig& config,
                                     ClauseDB& formula, Var num_concrete) {
    ThreadGroup<void> initial_pool;
    InitialHeuristicData initial{formula, num_concrete, &rec, config,
                                 &initial_pool};

    auto is_optimal = [&]() -> bool {
        std::size_t ub = initial.best_ub_value();
        std::size_t lb = initial.best_lb_value();
        return ub == lb || ub == 1;
    };

    auto begin = Clock::now();
    std::size_t iteration = 1;
    initial.first_solution();
    double secs;
    double secs_last = -100.0;
    while ((secs = seconds_between(begin, Clock::now())) < config.max_time &&
           !is_optimal() &&
           (iteration < config.goal_iterations || secs < config.min_time))
    {
        initial.follow_up_iteration(++iteration);
        if (secs - secs_last > 0.75) {
            rec.store_event("INITIAL_ITERATION_DONE",
                            {{"iteration", iteration}}, "iteration");
            secs_last = secs;
        } else {
            rec.store_event_silent("INITIAL_ITERATION_DONE",
                                   {{"iteration", iteration}});
        }
    }
    std::size_t universe_size = initial.infeasibility_map().count_vertices();
    std::vector<Vertex> coloring_order =
        initial.primal_solver().extract_coloring_order();
    rec.store_event("INITIAL_PHASE_DONE",
                    {{"ub", initial.get_solution().size()},
                     {"lb", initial.get_bound().size()},
                     {"num_interactions", universe_size}},
                    "lb", "ub", "num_interactions");
    InitialPhaseResult result{std::move(initial.infeasibility_map()),
                              initial.get_solution(),
                              initial.get_best_spawners(),
                              initial.get_bound(),
                              initial.get_all_spawners(),
                              std::move(coloring_order),
                              universe_size};
    return result;
}

/**
 * Cut-and-price portfolio element type.
 */
using CNPElement = PortfolioElementWithCore<CutAndPricePortfolioCore>;

/**
 * Exact solver (Sat + Clique + DSatur incremental) portfolio element type.
 */
using ExactElement = PortfolioElementWithCore<CliqueSatDSaturExactSolverCore>;

/**
 * Actual SAT solver to use.
 */
using SatSolver = ExternalNonIncrementalSAT<ExternalSolverType::KISSAT>;
// using SatSolver = CadicalSolver;

/**
 * MES + SAT as LNS elements.
 */
using LNSInner = sammy::FixedMESSATImprovementSolver<SatSolver>;
// using LNSInner =
// sammy::FixedMESIncrementalSATImprovementSolver<CadicalSolver>;
using LNSOuter = sammy::SubproblemSolverWithMES<LNSInner>;
using LNSCore = sammy::SubproblemLNSSolverCore<LNSOuter>;
using LNSElement = PortfolioElementWithCore<LNSCore>;
using OldLNSElement = CliqueSatDSaturLNSElement<CMSAT5Solver>;

class Main {
  public:
    struct Config {
        std::string input_file;
        std::string output_file;
        std::string dump_initial_phase_file;
        std::string subproblem_report_dir;
        bool all_concrete = false;
        bool print_events = false;
        bool print_portfolio_events = false;
        bool print_global_stats = false;
        bool print_initial_progress = false;
        bool dont_simplify = false;
        std::uint64_t implied_reduction_limit = 2'000'000;
        std::uint64_t exact_limit = 100'000;
        std::size_t max_lns_workers = std::thread::hardware_concurrency();
        std::size_t initial_goal_iterations = 10;
        double time_limit = 3600.0;
        double initial_time_limit = 1800.0;
        double initial_min_time = 5.0;

        void fill_options_description(po::options_description& opts,
                                      po::positional_options_description& pos) {
            opts.add_options()("input-file", required_value(input_file),
                               "The input file to read (only positional arg).");
            pos.add("input-file", 1);
            opts.add_options()("help,h", "produce help output")(
                "output-file,o", value_with_default(output_file),
                "The output file to write the results to.")(
                "dump-initial-phase",
                value_with_default(dump_initial_phase_file),
                "If set, dump the result of the initial phase (formula, "
                "universe, ...)"
                " to the given file, then terminate.")(
                "all-concrete", bool_switch(all_concrete),
                "If set, all features are considered concrete instead of the "
                "number "
                "given in the instance.")(
                "print-events", bool_switch(print_events),
                "Print time-stamped event information.")(
                "report-subproblems-to",
                value_with_default(subproblem_report_dir),
                "If set, subproblems are exported to this directory.")(
                "print-portfolio-events", bool_switch(print_portfolio_events),
                "If set, the portfolio solver parts print their events to "
                "stdout.")("print-global-stats",
                           bool_switch(print_global_stats),
                           "Print global statistics after solving.")(
                "print-initial-progress", bool_switch(print_initial_progress),
                "Print progress information during the initial heuristic.")(
                "dont-simplify", bool_switch(dont_simplify),
                "Do not apply simplification.")(
                "implied-reduction-limit",
                value_with_default(implied_reduction_limit),
                "Maximum number of valid interactions to attempt removal of "
                "implied "
                "interactions for.")(
                "exact-limit", value_with_default(exact_limit),
                "Maximum number of non-eliminated interactions to "
                "attempt an exact solve for.")(
                "max-lns-workers", value_with_default(max_lns_workers),
                "Maximum number of LNS workers to run in parallel.")(
                "time-limit", value_with_default(time_limit),
                "Time limit in seconds.")(
                "initial-time-limit", value_with_default(initial_time_limit),
                "Time limit for the initial heuristic.")(
                "initial-min-time", value_with_default(initial_min_time),
                "Time to spend on the initial heuristic at least.")(
                "initial-iteration-goal",
                value_with_default(initial_goal_iterations),
                "Iterations of the initial heuristic to run at least.");
        }

        static Config parse_options(int argc, char** argv) {
            Config config;
            po::options_description opts("Sammy options");
            po::positional_options_description pos;
            config.fill_options_description(opts, pos);
            po::variables_map vmap;
            try {
                if (!parse_cmdline(opts, pos, vmap, argc, argv)) {
                    // help was requested
                    std::cerr << opts << std::endl;
                    std::exit(0);
                }
            } catch (const std::exception& e) {
                // error parsing command line
                if (vmap.count("help")) {
                    std::cerr << opts << std::endl;
                    std::exit(0);
                } else {
                    std::cerr << "Error parsing command line: " << e.what()
                              << std::endl;
                    std::cerr << opts << std::endl;
                    std::exit(1);
                }
            }
            config.initial_time_limit =
                (std::min)(config.initial_time_limit, config.time_limit);
            config.initial_min_time =
                (std::min)(config.initial_min_time, config.time_limit / 3.0);
            return config;
        }
    };

    Main(int argc, char** argv)
        : config(Config::parse_options(argc, argv)), recorder() {
        recorder.set_print_events(config.print_events);
        recorder.store_event("COMMAND_LINE_PARSED");
    }

    void read_input() {
        try {
            std::filesystem::path input_path(config.input_file);
            input = ::read_input(input_path);
            if (config.all_concrete)
                input->num_concrete = input->formula.num_vars();
            recorder.store_event("INPUT_READ",
                                 {{"input_file", config.input_file},
                                  {"num_vars", input->formula.num_vars()},
                                  {"num_clauses", input->formula.num_clauses()},
                                  {"num_concrete", input->num_concrete}},
                                 "input_file", "num_vars", "num_clauses",
                                 "num_concrete");
        } catch (const std::exception& e) {
            std::cerr << "Failed to read input file: " << e.what() << std::endl;
            std::exit(1);
        }
    }

    void simplify_if_configured() {
        actual_clause_db = &input->formula;
        actual_num_concrete = input->num_concrete;
        if (!config.dont_simplify) {
            run_simplification(recorder, *input, simplifier, simplified);
            simplified.emplace(remove_subsumed(*simplified));
            actual_clause_db = &simplified->formula;
            actual_num_concrete = simplified->num_concrete;
        }
    }

    void dump_initial_phase_result() {
        OutputObject data_dump =
            export_initial_phase_result(*actual_clause_db, *initial_result);
        std::filesystem::path dump_path(config.dump_initial_phase_file);
        std::ofstream dump_file(dump_path, std::ios::out | std::ios::trunc);
        dump_file << data_dump;
        std::exit(0);
    }

    bool initial_phase() {
        try {
            RunInitialConfig initial_config{
                config.initial_min_time, config.initial_time_limit,     20, 10,
                !config.dont_simplify,   !config.print_initial_progress};
            initial_result =
                run_initial_phase(recorder, initial_config, *actual_clause_db,
                                  actual_num_concrete);
            if (!config.dump_initial_phase_file.empty()) {
                dump_initial_phase_result();
            }
            return initial_result->best_solution.size() == 1 ||
                   initial_result->best_solution.size() ==
                       initial_result->best_mutually_exclusive.size();
        } catch (const UNSATError&) {
            std::cerr << "The feature model is infeasible (UNSAT)!"
                      << std::endl;
            std::exit(0);
        }
    }

    void universe_reduction() {
        const auto& implied_cache = portfolio->implied_cache();
        std::size_t original_size = implied_cache.original_universe_size();
        if (original_size < config.implied_reduction_limit) {
            portfolio->reduce_universe();
        } else {
            double t = (std::min)(10.0, 0.02 * config.time_limit);
            portfolio->limited_reduce_universe(t);
        }
        portfolio->add_clique(portfolio->get_best_mes());
    }

    void setup_element(PortfolioElement& element) {
        if (config.print_portfolio_events) {
            element.set_recorder_quiet(false);
        }
        element.synchronize_recorder(recorder);
    }

    void add_exact_solver() {
        const auto& implied = portfolio->implied_cache();
        if (implied.get_reduced_universe().size() < config.exact_limit) {
            have_exact = true;
            ExactElement& exact = portfolio->emplace_element<ExactElement>(
                &*portfolio, &CliqueSatDSaturExactSolverCore::factory,
                "Exact Clique & SatDSatur");
            setup_element(exact);
        }
    }

    void add_cnp_element() {
        CNPElement& element = portfolio->emplace_element<CNPElement>(
            &*portfolio, &CutAndPricePortfolioCore::factory, "Cut & Price");
        setup_element(element);
    }

    void add_lns_elements() {
        std::size_t remaining_threads = std::thread::hardware_concurrency() - 1;
        remaining_threads -= std::size_t(have_exact);
        remaining_threads =
            (std::min)(remaining_threads, config.max_lns_workers);
        for (std::size_t i : range(remaining_threads)) {
            static_cast<void>(i); // not used
            setup_element(portfolio->emplace_element<LNSElement>(
                &*portfolio, &LNSCore::factory, "Non-Incremental SAT LNS"));
        }
    }

    std::vector<std::vector<bool>>
    reconstruct_sample(const std::vector<std::vector<bool>>& sample) {
        if (simplifier && simplified) {
            return simplifier->reconstruct_sample(*simplified, sample);
        } else {
            return sample;
        }
    }

    std::vector<Vertex> reconstruct_mes(const std::vector<Vertex>& mes) {
        if (simplifier && simplified) {
            return simplifier->reconstruct_lb(*simplified, mes);
        } else {
            return mes;
        }
    }

    void handle_output_common() {
        recorder.store_event("CONSTRUCTING_OUTPUT_DATA");
        output["experiment_tag"] = "sammy_solve";
        output["instance_name"] = input->name;
        output["initial_min_time"] = config.initial_min_time;
        output["initial_time_limit"] = config.initial_time_limit;
        output["initial_iteration_goal"] = config.initial_goal_iterations;
        output["num_threads"] = std::thread::hardware_concurrency();
    }

    void handle_output_portfolio() {
        handle_output_common();
        output["ran_exact_solver"] = have_exact;
        output["lb"] = portfolio->get_best_lower_bound();
        output["ub"] = portfolio->get_best_solution_size();
        output["optimal"] = portfolio->get_best_lower_bound() ==
                            portfolio->get_best_solution_size();
        auto assignments =
            portfolio->get_best_solution().assignments_as<std::vector<bool>>();
        output["best_solution"] = reconstruct_sample(assignments);
        output["mutually_exclusive_set"] =
            lit::externalize(reconstruct_mes(portfolio->get_best_mes()));
        output["lb_vertices"] = lit::externalize(
            reconstruct_mes(portfolio->get_best_lower_bound_vertex_set()));
        export_events(output, recorder.events());
        if (!config.output_file.empty()) {
            output_data(output, config.output_file);
        }
        if (config.print_global_stats) {
            std::cout << get_global_stats() << std::endl;
        }
    }

    void handle_output_initial_phase_optimal() {
        handle_output_common();
        output["ran_exact_solver"] = false;
        output["lb"] = initial_result->best_solution.size();
        output["ub"] = initial_result->best_solution.size();
        output["optimal"] = true;
        output["best_solution"] =
            reconstruct_sample(initial_result->best_solution);
        output["mutually_exclusive_set"] = lit::externalize(
            reconstruct_mes(initial_result->best_mutually_exclusive));
        output["lb_vertices"] = output["mutually_exclusive_set"];
        export_events(output, recorder.events());
        if (!config.output_file.empty()) {
            output_data(output, config.output_file);
        }
    }

    void main() {
        read_input();
        simplify_if_configured();
        if (initial_phase()) {
            // optimality after initial phase
            recorder.store_event("OPTIMALITY_REACHED");
            handle_output_initial_phase_optimal();
            return;
        }
        // setup clauses
        clause_db_ticket = publish_clauses(*actual_clause_db);
        // setup portfolio
        portfolio.emplace(clause_db_ticket, &recorder,
                          std::move(*initial_result));
        universe_reduction();
        add_exact_solver();
        add_cnp_element();
        add_lns_elements();
        if (!config.subproblem_report_dir.empty()) {
            portfolio->enable_subproblem_reporting(
                config.subproblem_report_dir);
        }
        recorder.store_event("PORTFOLIO_START");
        portfolio->run(config.time_limit - recorder.now());
        portfolio->merge_events();
        recorder.store_event("PORTFOLIO_DONE");
        handle_output_portfolio();
    }

  private:
    // the configuration
    Config config;
    // the event recorder
    EventRecorder recorder;
    // input data, if already constructed
    std::optional<InputData> input;
    // simplification data and result, if we are simplifying
    std::optional<SimplifyDatastructure> simplifier;
    std::optional<SimplifiedInstance> simplified;
    // the clause DB and number of concrete features actually used
    ClauseDB* actual_clause_db;
    Var actual_num_concrete;
    // initial phase result, if we have it (after initial phase,
    // before constructing portfolio)
    std::optional<InitialPhaseResult> initial_result;
    // the clause database ticket, once initial phase is done
    ClausesTicket clause_db_ticket;
    // the portfolio solver, once its constructed
    std::optional<PortfolioSolver> portfolio;
    // whether we have an exact solver
    bool have_exact{false};
    // output data
    OutputObject output;
};

int main(int argc, char** argv) {
    // check if we are running a SAT-only call
    auto sat_only = is_sat_only_call(argc, argv);
    if (sat_only) {
        return sat_only_entry_point(*sat_only);
    }

    // run the main program
    Main main(argc, argv);
    main.main();
    return 0;
}
