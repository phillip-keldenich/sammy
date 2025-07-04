#include "io.h"
#include <functional>
#include <sammy/barrage.h>
#include <sammy/cadical_solver.h>
#include <sammy/clique_sat_dsatur.h>
#include <sammy/cmsat5_solver.h>
#include <sammy/incremental_sat_lns.h>
#include <sammy/kissat_solver.h>
#include <sammy/lingeling_solver.h>
#include <sammy/output.h>
#include <sammy/sat_dsatur.h>
#include <sammy/sat_lns.h>
#include <sammy/subproblem_solver_with_mes.h>

using namespace sammy;

namespace po = boost::program_options;

/**
 * Result from solving a subproblem.
 * Contains the outcome of the solve, the configurations
 * found, and the time taken to build the solver and solve the subproblem.
 */
struct SolveResult {
    std::string outcome;
    std::vector<DynamicBitset> configurations;
    double build_time;
    double solve_time;

    std::optional<std::vector<std::vector<bool>>> get_assignments() const {
        if (outcome != "IMPROVED_SOLUTION") {
            return std::nullopt;
        }
        std::vector<std::vector<bool>> assignments;
        std::transform(configurations.begin(), configurations.end(),
                       std::back_inserter(assignments),
                       [](const DynamicBitset& config) {
                           return static_cast<std::vector<bool>>(config);
                       });
        return assignments;
    }
};

/**
 * Run a solver with a timeout outside of the usual
 * portfolio structure that would enforce the timeout.
 */
template <typename TypeWithSolveAndAbort>
auto timeout_wrapped_solve(TypeWithSolveAndAbort& solver, double time_limit) {
    std::mutex m;
    std::condition_variable cv;
    bool done = false;
    auto timeout_watcher_routine = [&]() {
        std::unique_lock l{m};
        if (!cv.wait_for(l, std::chrono::duration<double>(time_limit),
                         [&]() { return done; }))
        {
            solver.abort();
        }
    };
    std::thread timeout_watcher(timeout_watcher_routine);
    auto res = solver.solve();
    {
        std::unique_lock l{m};
        done = true;
        cv.notify_one();
    }
    timeout_watcher.join();
    return res;
}

/**
 * Create and run a common interface subproblem solver.
 */
template <typename SolverType>
std::pair<LNSSubproblem, SolveResult>
make_and_run_common_interface_solver(LNSSubproblem&& subproblem,
                                     ClauseDB* clauses, EventRecorder* recorder,
                                     double time_limit) {
    auto before_make = std::chrono::steady_clock::now();
    SolverType solver(nullptr, std::move(subproblem),
                      SharedDBPropagator(clauses), recorder, 0);
    auto after_make = std::chrono::steady_clock::now();
    auto res = timeout_wrapped_solve(solver, time_limit);
    auto after_solve = std::chrono::steady_clock::now();
    SolveResult result;
    result.build_time = seconds_between(before_make, after_make);
    result.solve_time = seconds_between(after_make, after_solve);
    if (!res) {
        result.outcome = "TIMEOUT";
    } else if (!*res) {
        result.outcome = "SOLUTION_WAS_OPTIMAL";
    } else {
        result.outcome = "IMPROVED_SOLUTION";
        result.configurations = solver.get_solution();
    }
    return {solver.move_out_subproblem(), std::move(result)};
}

/**
 * Write the results of the solve to a JSON file.
 * The JSON file will contain the outcome of the solve,
 * the improved configurations (if any), the build time, and the solve time,
 * as well as any events recorded during the solve.
 */
void write_result(SolveResult result, const EventRecorder& recorder,
                  const std::string& output_file_path) {
    nlohmann::json output_json;
    output_json["outcome"] = result.outcome;
    auto assignments = result.get_assignments();
    if (assignments) {
        output_json["assignments"] = *result.get_assignments();
    } else {
        output_json["assignments"] = nlohmann::json{};
    }
    output_json["build_time"] = result.build_time;
    output_json["solve_time"] = result.solve_time;
    export_events(output_json, recorder);
    write_json_path(output_file_path, output_json);
}

/**
 * Run a specific LNS subproblem solver on a
 * given subproblem (given as *amended* subproblem file),
 * produced by amend_subproblem_with_default_mes, which
 * contains a problem, subproblem and a default MES.
 */
int main(int argc, char** argv) {
    double solve_time_limit = 1800.0; //< time limit for the subproblem solver
    std::string solver_name;          //< name of the solver to run
    std::string amended_subproblem;   //< amended subproblem file
    std::string output_file_path;     //< path to the output JSON file

    // command line options
    boost::program_options::options_description desc("Options");
    desc.add_options()("amended-subproblem-file",
                       required_value(amended_subproblem),
                       "Path to the MES entry point file.")(
        "solver-name", required_value(solver_name),
        "Name of the solver to run.")("solve-time-limit",
                                      value_with_default(solve_time_limit),
                                      "Time limit for the solver.")(
        "output-file,o", value_with_default(output_file_path),
        "Path to the output JSON file.")("help,h", "Print this help message.");

    // parse command line
    try {
        po::variables_map vmap;
        auto cmdline_parser = boost::program_options::command_line_parser(
            argc, const_cast<const char* const*>(argv));
        po::store(cmdline_parser.options(desc).run(), vmap);
        vmap.notify();
        if (vmap.count("help")) {
            std::cerr << desc << std::endl;
            return 0;
        }
    } catch (const std::exception& e) {
        std::cerr << desc << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }

    /**
     * Parse input (problem and subproblem), and prepare the input.
     */
    nlohmann::json input = read_json_path(amended_subproblem);
    const nlohmann::json& amendment = input.at("amendment");
    ProblemInput problem_input =
        interpret_problem_json(amendment.at("problem_input"));
    const nlohmann::json& json_mes = amendment.at("mes");
    auto mes = lit::internalize(json_mes.get<std::vector<ExternalVertex>>());
    SubproblemInput subproblem_input =
        interpret_subproblem_json(std::move(input));

    /**
     * Setup event recorder and subproblem.
     */
    LNSSubproblem subproblem{
        subproblem_input.uncovered_vertices, mes,
        subproblem_input.covering_assignments,
        static_cast<Lit>(problem_input.initial_phase.inf_map.get_n_concrete())};

    /**
     * Handle cases usually handled elsewhere
     * in the LNS portfolio where the MES precludes any improvement.
     */
    if (subproblem.mutually_exclusive_set.size() >=
        subproblem_input.covering_assignments.size())
    {
        EventRecorder recorder;
        recorder.set_print_events(false);
        recorder.store_event("MES_PRECLUDES_IMPROVEMENT");
        write_result({"SOLUTION_WAS_OPTIMAL", {}, 0.0, 0.0}, recorder,
                     output_file_path);
        return 0;
    }

    /**
     * Make and run a common interface solver
     * (new SatDSatur, incremental SAT, or fixed SAT solver).
     */
    auto run_solver_type = [&](auto solver_nullptr) {
        using SolverType = std::remove_pointer_t<decltype(solver_nullptr)>;
        EventRecorder recorder;
        recorder.set_print_events(false);
        auto subproblem_and_result =
            make_and_run_common_interface_solver<SolverType>(
                std::move(subproblem), &problem_input.clauses, &recorder,
                solve_time_limit);
        subproblem = std::move(subproblem_and_result.first);
        write_result(std::move(subproblem_and_result.second), recorder,
                     output_file_path);
    };

    // run a CNPSATDSatur-based solver
    auto run_cnpsatdsatur = [&](auto solver_nullptr) {
        using SatSolverType = std::remove_pointer_t<decltype(solver_nullptr)>;
        using SolverType = CliqueSatDSaturSolver<SatSolverType>;
        auto before_make = std::chrono::steady_clock::now();
        EventRecorder recorder;
        SolverType solver{
            subproblem_input.uncovered_vertices,
            &problem_input.initial_phase.inf_map,
            problem_input.clauses,
            mes,
            problem_input.initial_phase.best_mutually_exclusive.size(),
            problem_input.initial_phase.best_mutually_exclusive.size(),
            subproblem_input.covering_assignments};
        solver.set_event_recorder(&recorder);
        auto after_make = std::chrono::steady_clock::now();
        auto res = timeout_wrapped_solve(solver, solve_time_limit);
        auto after_solve = std::chrono::steady_clock::now();
        std::vector<DynamicBitset> solution;
        if (res == SolverType::SolveResult::IMPROVED_SOLUTION) {
            solution = solver.get_partial_solution().assignments();
        }
        SolveResult solve_result{
            (res == SolverType::SolveResult::ABORTED ? "TIMEOUT"
             : res == SolverType::SolveResult::SOLUTION_WAS_OPTIMAL
                 ? "SOLUTION_WAS_OPTIMAL"
                 : "IMPROVED_SOLUTION"),
            std::move(solution), seconds_between(before_make, after_make),
            seconds_between(after_make, after_solve)};
        write_result(std::move(solve_result), recorder, output_file_path);
    };

    // solver name table and definitions
#define FIXED_SOLVER(SolverType)                                               \
    [&]() {                                                                    \
        using Solver = FixedMESSATImprovementSolver<SolverType>;               \
        run_solver_type(static_cast<Solver*>(nullptr));                        \
    }

#define INCRE_SOLVER(SolverType)                                               \
    [&]() {                                                                    \
        using Solver = FixedMESIncrementalSATImprovementSolver<SolverType>;    \
        run_solver_type(static_cast<Solver*>(nullptr));                        \
    }

#define CNPSD_SOLVER(SolverType)                                               \
    [&]() { run_cnpsatdsatur(static_cast<SolverType*>(nullptr)); }

#define NSDS_SOLVER(SolverType)                                                \
    [&]() {                                                                    \
        using Solver = FixedMESSatDSaturSolver<SolverType>;                    \
        run_solver_type(static_cast<Solver*>(nullptr));                        \
    }

    std::unordered_map<std::string, std::function<void()>> solvers{
        {"fixed_sat[kissat]", FIXED_SOLVER(KissatSolver)},
        {"fixed_sat[cadical]", FIXED_SOLVER(CadicalSolver)},
        {"fixed_sat[cryptominisat]", FIXED_SOLVER(CMSAT5Solver)},
        {"fixed_sat[lingeling]", FIXED_SOLVER(LingelingSolver)},
        {"incremental_sat[cadical]", INCRE_SOLVER(CadicalSolver)},
        {"incremental_sat[cryptominisat]", INCRE_SOLVER(CMSAT5Solver)},
        {"incremental_sat[lingeling]", INCRE_SOLVER(LingelingSolver)},
        {"satdsatur[cadical]", CNPSD_SOLVER(CadicalSolver)},
        {"satdsatur[cryptominisat]", CNPSD_SOLVER(CMSAT5Solver)},
        {"satdsatur[lingeling]", CNPSD_SOLVER(LingelingSolver)},
        {"newsatdsatur[cadical]", NSDS_SOLVER(CadicalSolver)},
        {"newsatdsatur[cryptominisat]", NSDS_SOLVER(CMSAT5Solver)},
        {"newsatdsatur[lingeling]", NSDS_SOLVER(LingelingSolver)}};
#undef NSDS_SOLVER
#undef CNPSD_SOLVER
#undef INCRE_SOLVER
#undef FIXED_SOLVER

    // look up and run the solver
    if (solvers.count(solver_name) == 0) {
        std::cerr << "Unknown solver name: " << solver_name << "\n";
        std::cerr << "Available solvers:\n";
        for (const auto& e : solvers) {
            std::cerr << "    " << e.first << "\n";
        }
        return 3;
    }
    auto& solver_func = solvers.at(solver_name);
    solver_func();
    return 0;
}
