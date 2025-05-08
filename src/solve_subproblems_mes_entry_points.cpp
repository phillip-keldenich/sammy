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
#include <sammy/sat_lns.h>
#include <sammy/subproblem_solver_with_mes.h>

using namespace sammy;

namespace po = boost::program_options;

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

int main(int argc, char** argv) {
    double solve_time_limit = 1800.0;
    std::string solver_name;
    std::string mes_entry_point_file;
    std::string output_file_path;
    std::size_t num_runs = 3;
    std::size_t target_mes_size = 0;

    // command line options
    boost::program_options::options_description desc("Options");
    desc.add_options()("entry-point-file", required_value(mes_entry_point_file),
                       "Path to the MES entry point file.")(
        "solver-name", required_value(solver_name),
        "Name of the solver to run.")("mes-size",
                                      required_value(target_mes_size),
                                      "The MES size to run the solver for.")(
        "num-runs", value_with_default(num_runs),
        "Number of times to run the solver.")(
        "solve-time-limit", value_with_default(solve_time_limit),
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

    // read and parse input
    nlohmann::json input = read_json_path(mes_entry_point_file);
    ProblemInput problem_input =
        interpret_problem_json(input.at("problem_input"));
    SubproblemInput subproblem_input =
        interpret_subproblem_json(input.at("subproblem_input"));

    // initialize with first MES
    std::vector<std::vector<Vertex>> meses;
    auto internalize = [](const auto& j) {
        return lit::internalize(j.template get<std::vector<ExternalVertex>>());
    };
    meses.push_back(internalize(input.at("initial_mes")));

    // MES history: each entry is a new, better MES found
    const auto& mes_history_data = input.at("mes_history");
    std::transform(mes_history_data.begin(), mes_history_data.end(),
                   std::back_inserter(meses), [&](const auto& entry) {
                       return internalize(entry.at("mes"));
                   });

    // find the target MES
    auto mes_found =
        std::find_if(meses.begin(), meses.end(), [&](const auto& mes) {
            return mes.size() == target_mes_size;
        });
    if (mes_found == meses.end()) {
        std::cerr << "Could not find a MES of the given size to solve!\n";
        return 2;
    }

    // the actual MES to work on
    std::vector<Vertex> mes = *mes_found;
    // and a subproblem based on it
    LNSSubproblem subproblem{
        subproblem_input.uncovered_vertices, mes,
        subproblem_input.covering_assignments,
        static_cast<Lit>(problem_input.initial_phase.inf_map.get_n_concrete())};

    // keys in each output row:
    // build_time, configurations, events, outcome, solve_time
    OutputObject output_rows;

    // add another row to output_rows
    auto add_outcome_row = [&](SolveResult&& result,
                               const EventRecorder& recorder) {
        OutputObject entry;
        auto assignments = result.get_assignments();
        entry["outcome"] = std::move(result.outcome);
        entry["build_time"] = result.build_time;
        entry["solve_time"] = result.solve_time;
        if (assignments) {
            entry["configurations"] = std::move(*assignments);
        } else {
            entry["configurations"] = OutputObject{};
        }
        export_events(entry, recorder.events());
        output_rows.push_back(std::move(entry));
    };

    // run a SAT-based solver (either incremental or fixed)
    auto run_solver_type = [&](auto solver_nullptr) {
        using SolverType = std::remove_pointer_t<decltype(solver_nullptr)>;
        EventRecorder recorder;
        auto subproblem_and_result =
            make_and_run_common_interface_solver<SolverType>(
                std::move(subproblem), &problem_input.clauses, &recorder,
                solve_time_limit);
        subproblem = std::move(subproblem_and_result.first);
        add_outcome_row(std::move(subproblem_and_result.second), recorder);
    };

    // run a CNPSATDSatur-based solver
    auto run_cnpsatdsatur = [&](auto solver_nullptr) {
        using SatSolverType = std::remove_pointer_t<decltype(solver_nullptr)>;
        using SolverType = CliqueSatDSaturSolver<SatSolverType>;
        auto before_make = std::chrono::steady_clock::now();
        SolverType solver{
            subproblem_input.uncovered_vertices,
            &problem_input.initial_phase.inf_map,
            problem_input.clauses,
            mes,
            problem_input.initial_phase.best_mutually_exclusive.size(),
            problem_input.initial_phase.best_mutually_exclusive.size(),
            subproblem_input.covering_assignments,
            false};
        EventRecorder recorder;
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
        add_outcome_row(std::move(solve_result), recorder);
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
        {"satdsatur[lingeling]", CNPSD_SOLVER(LingelingSolver)}};
#undef CNPSD_SOLVER
#undef INCRE_SOLVER
#undef FIXED_SOLVER

    // look up and run the solver
    if (solvers.count(solver_name) == 0) {
        std::cerr << "Unknown solver name: " << solver_name << "\n";
        return 3;
    }
    auto& solver_func = solvers.at(solver_name);
    for (std::size_t i = 0; i < num_runs; ++i) {
        solver_func();
    }

    // write output
    if (!output_file_path.empty()) {
        write_json_path(output_file_path, output_rows);
    } else {
        std::cout << output_rows.dump(2) << std::endl;
    }
    return 0;
}
