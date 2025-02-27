#include <sammy/barrage.h>
#include <sammy/barrage_worker_cnp.h>
#include <sammy/barrage_worker_exact.h>
#include <sammy/barrage_worker_lns.h>
#include <sammy/experiment_flags.h>
#include <sammy/incremental_sat_lns.h>
#include <sammy/kissat_solver.h>
#include <sammy/lingeling_solver.h>
#include <sammy/output.h>
#include <sammy/partial_solution.h>
#include <sammy/run_initial.h>
#include <sammy/sat_lns.h>
#include <sammy/subproblem_solver_with_mes.h>
#include <sammy/thread_clauses.h>

using namespace sammy;

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
                                     SimplifiedInstance& instance,
                                     const ExperimentFlagsConfig& config) {
    const auto& iconf = config.initial_heuristic_config;
    ThreadGroup<void> initial_pool;
    InitialHeuristicData initial{instance.formula, instance.num_concrete, &rec,
                                 config.initial_heuristic_config,
                                 &initial_pool};
    auto begin = Clock::now();
    std::size_t iteration = 1;
    initial.first_solution();
    double secs;
    double secs_last = -100.0;
    while ((secs = seconds_between(begin, Clock::now())) < iconf.max_time &&
           initial.best_lb_value() != initial.best_ub_value() &&
           (iteration < iconf.goal_iterations || secs < iconf.min_time))
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
 * MES + Incremental SAT as LNS elements.
 */
using LNSInner = sammy::FixedMESSATImprovementSolver<CadicalSolver>;
// using LNSInner =
// sammy::FixedMESIncrementalSATImprovementSolver<CadicalSolver>;
using LNSOuter = sammy::SubproblemSolverWithMES<LNSInner>;
using LNSCore = sammy::SubproblemLNSSolverCore<LNSOuter>;
using LNSElement = PortfolioElementWithCore<LNSCore>;
using OldLNSElement = CliqueSatDSaturLNSElement<CMSAT5Solver>;

int main(int argc, char** argv) {
    // set defaults for initial phase
    ExperimentFlagsConfig config;
    config.initial_heuristic_config.goal_iterations = 10;
    config.initial_heuristic_config.min_time = 5.0;
    config.initial_heuristic_config.max_time = 1800.0;
    double time_limit = 3600.0;
    std::uint64_t implied_reduction_limit = 2'000'000;
    std::uint64_t exact_limit = 100'000;
    bool all_concrete = false;
    bool print_portfolio_events = false;
    bool old_lns = false;
    std::string subproblem_report_dir;
    std::string dump_initial_phase_file;
    std::size_t max_lns_workers = std::thread::hardware_concurrency();
    boost::program_options::options_description extra_options(
        "Barrage options");
    extra_options.add_options()("time-limit", value_with_default(time_limit),
                                "Time limit for the portfolio solver.")(
        "implied-reduction-limit", value_with_default(implied_reduction_limit),
        "Maximum number of valid interactions to attempt removal of implied "
        "interactions for.")("exact-limit", value_with_default(exact_limit),
                             "Maximum number of non-eliminated interactions to "
                             "attempt an exact solve for.")(
        "all-concrete", bool_switch(all_concrete),
        "If set, all features are considered concrete instead of the number "
        "given in the instance.")(
        "report-subproblems-to", value_with_default(subproblem_report_dir),
        "If set, subproblems are exported to this directory.")(
        "print-portfolio-events", bool_switch(print_portfolio_events),
        "If set, the portfolio solver parts print their events to stdout.")(
        "dump-initial-phase", value_with_default(dump_initial_phase_file),
        "If set, dump the result of the initial phase (formula, universe, ...)"
        " to the given file, then terminate.")(
        "old-lns", bool_switch(old_lns),
        "use the old LNS worker (CliqueSatDSatur)")(
        "max-lns-workers", value_with_default(max_lns_workers),
        "Maximum number of LNS workers to run.");
    get_config_or_exit(argc, argv, config, extra_options, true, false);

    config.initial_heuristic_config.max_time =
        std::min(config.initial_heuristic_config.max_time, time_limit);
    config.initial_heuristic_config.min_time =
        std::min(config.initial_heuristic_config.min_time, time_limit / 3.0);

    // read input and start recording events/timings
    std::filesystem::path input_path(config.input_file);
    EventRecorder recorder;
    recorder.set_print_events(config.print_events);
    InputData input = read_input(input_path);
    recorder.store_event("INPUT_READ");
    if (all_concrete)
        input.num_concrete = input.formula.num_vars();

    // run simplification on the given SAT formula
    std::optional<SimplifyDatastructure> simplifier;
    std::optional<SimplifiedInstance> simplified;
    run_simplification(recorder, input, simplifier, simplified);

    // run the initial heuristics multiple times
    InitialPhaseResult initial_result =
        run_initial_phase(recorder, *simplified, config);
    // remove any subsumed clauses (the initial heuristics learn clauses!)
    simplified.emplace(remove_subsumed(*simplified));

    // 'publish' the clause database so each
    // worker thread gets its own local copy
    // that it can freely change
    auto clause_db_ticket = publish_clauses(simplified->formula);

    // dump the initial phase result if requested
    if (!dump_initial_phase_file.empty()) {
        OutputObject data_dump =
            export_initial_phase_result(simplified->formula, initial_result);
        std::filesystem::path dump_path(dump_initial_phase_file);
        std::ofstream dump_file(dump_path, std::ios::out | std::ios::trunc);
        dump_file << data_dump;
        return 0;
    }

    // setup & run portfolio
    PortfolioSolver solver{clause_db_ticket, &recorder,
                           std::move(initial_result)};
    bool have_exact = false;
    if (initial_result.universe_size < implied_reduction_limit) {
        solver.reduce_universe();
    } else {
        solver.limited_reduce_universe((std::min)(10.0, 0.02 * time_limit));
    }
    solver.add_clique(solver.get_best_mes());

    auto setup_element = [&](PortfolioElement& element) {
        if (print_portfolio_events) {
            element.set_recorder_quiet(false);
        }
        element.synchronize_recorder(recorder);
    };
    if (solver.implied_cache().have_reduced_universe() &&
        solver.implied_cache().get_reduced_universe().size() < exact_limit)
    {
        have_exact = true;
        setup_element(solver.emplace_element<ExactElement>(
            &solver, &CliqueSatDSaturExactSolverCore::factory,
            "Exact Clique & SatDSatur"));
    }
    setup_element(solver.emplace_element<CNPElement>(
        &solver, &CutAndPricePortfolioCore::factory, "Cut & Price"));
    std::size_t remaining_threads =
        std::thread::hardware_concurrency() - 1 - std::size_t(have_exact);
    if (remaining_threads > max_lns_workers) {
        remaining_threads = max_lns_workers;
    }
    for (std::size_t i = 0; i < remaining_threads; ++i) {
        if (old_lns) {
            setup_element(solver.emplace_element<OldLNSElement>(
                &solver, solver.get_best_solution_size(), i));
        } else {
            setup_element(solver.emplace_element<LNSElement>(
                &solver, &LNSCore::factory, "Non-Incremental SAT LNS"));
        }
    }
    if (!subproblem_report_dir.empty()) {
        solver.enable_subproblem_reporting(subproblem_report_dir);
    }
    recorder.store_event("PORTFOLIO_START");
    solver.run(time_limit - recorder.now());
    solver.merge_events();
    recorder.store_event("PORTFOLIO_DONE");

    // output data/statistics
    OutputObject output;
    output["experiment_tag"] = "run_barrage";
    output["instance_name"] = input.name;
    output["initial_min_time"] = config.initial_heuristic_config.min_time;
    output["initial_time_limit"] = config.initial_heuristic_config.max_time;
    output["initial_iteration_goal"] =
        config.initial_heuristic_config.goal_iterations;
    output["num_threads"] = std::thread::hardware_concurrency();
    output["ran_exact_solver"] = have_exact;
    output["lb"] = solver.get_best_lower_bound();
    output["ub"] = solver.get_best_solution_size();
    output["optimal"] =
        (solver.get_best_lower_bound() == solver.get_best_solution_size());
    auto assignments =
        solver.get_best_solution().assignments_as<std::vector<bool>>();
    output["best_solution"] =
        simplifier->reconstruct_sample(*simplified, assignments);
    output["mutually_exclusive_set"] = lit::externalize(
        simplifier->reconstruct_lb(*simplified, solver.get_best_mes()));
    output["lb_vertices"] = lit::externalize(simplifier->reconstruct_lb(
        *simplified, solver.get_best_lower_bound_vertex_set()));
    export_events(output, recorder.events());
    if (!config.output_file.empty())
        output_data(output, config.output_file);
    return 0;
}
