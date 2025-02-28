#include "io.h"
#include <sammy/barrage.h>
#include <sammy/output.h>
#include <sammy/subproblem_solver_with_mes.h>
#include <sammy/sat_lns.h>
#include <sammy/incremental_sat_lns.h>
#include <sammy/clique_sat_dsatur.h>
#include <sammy/kissat_solver.h>
#include <sammy/cadical_solver.h>
#include <sammy/cmsat5_solver.h>
#include <sammy/lingeling_solver.h>
#include <functional>

using namespace sammy;

using SatSolver = KissatSolver;

namespace po = boost::program_options;

struct SolveResult {
    std::string outcome;
    std::vector<DynamicBitset> configurations;
    double build_time;
    double solve_time;

    std::optional<std::vector<std::vector<bool>>> get_assignments() const {
        if(outcome != "IMPROVED_SOLUTION") {
            return std::nullopt;
        }
        std::vector<std::vector<bool>> assignments;
        std::transform(configurations.begin(), configurations.end(), std::back_inserter(assignments),
            [](const DynamicBitset& config) { return static_cast<std::vector<bool>>(config); }
        );
        return assignments;
    }
};

template<typename TypeWithSolveAndAbort>
auto timeout_wrapped_solve(TypeWithSolveAndAbort& solver, double time_limit) {
    std::mutex m;
    std::condition_variable cv;
    bool done = false;
    auto timeout_watcher_routine = [&]() {
        std::unique_lock l{m};
        if (!cv.wait_for(l, std::chrono::duration<double>(time_limit), [&]() { return done; })) {
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

template<typename SolverType>
std::pair<LNSSubproblem, SolveResult> 
    make_and_run_common_interface_solver(LNSSubproblem&& subproblem, ClauseDB* clauses, 
                                         EventRecorder* recorder, double time_limit) 
{
    auto before_make = std::chrono::steady_clock::now();
    SolverType solver(nullptr, std::move(subproblem), SharedDBPropagator(clauses), recorder, 0);
    auto after_make = std::chrono::steady_clock::now();
    auto res = timeout_wrapped_solve(solver, time_limit);
    auto after_solve = std::chrono::steady_clock::now();
    SolveResult result;
    result.build_time = seconds_between(before_make, after_make);
    result.solve_time = seconds_between(after_make, after_solve);
    if(!res) {
        result.outcome = "TIMEOUT";
    }
    if(!*res) {
        result.outcome = "SOLUTION_WAS_OPTIMAL";
    }
    if(*res) {
        result.outcome = "IMPROVED_SOLUTION";
        result.configurations = solver.get_solution();
    }
    return {solver.move_out_subproblem(), std::move(result)};
}

int main(int argc, char** argv) {
    std::string mes_entry_point_file;
    std::string output_file_path;
    double solve_time_limit = 1800.0;

    // command line options
    boost::program_options::options_description desc("Options");
    desc.add_options()("entry-point-file", required_value(mes_entry_point_file),
                       "Path to the MES entry point file.")
                      ("solve-time-limit", value_with_default(solve_time_limit),
                       "Time limit for the solver.")
                      ("output-file,o", value_with_default(output_file_path),
                       "Path to the output JSON file.")
                      ("help,h", "Print this help message.");
    
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
    } catch(const std::exception& e) {
        std::cerr << desc << std::endl;
        std::cerr << e.what() << std::endl;
        return 1;
    }
    
    // read MES entry point file
    std::ifstream input_file(mes_entry_point_file.c_str(), std::ios::in);
    auto input = nlohmann::json::parse(input_file);

    ProblemInput problem_input = interpret_problem_json(input.at("problem_input"));
    SubproblemInput subproblem_input = interpret_subproblem_json(input.at("subproblem_input"));
    OutputObject mes_history_data = input.at("mes_history");
    OutputObject solver_run_data;
    OutputObject overall_output{
        {"formula_file", input.at("formula_file")},
        {"subproblem_file", input.at("subproblem_file")},
        {"entry_point_file", mes_entry_point_file},
        {"solve_time_limit", solve_time_limit},
        {"solver_runs", OutputObject{}}
    };

    auto add_solver_outcome = [&] (const char* solver_name, SolveResult&& result, const EventRecorder& recorder) {
        OutputObject entry;
        auto assignments = result.get_assignments();
        entry["outcome"] = std::move(result.outcome);
        entry["build_time"] = result.build_time;
        entry["solve_time"] = result.solve_time;
        if(assignments) {
            entry["configurations"] = std::move(*assignments);
        } else {
            entry["configurations"] = OutputObject{};
        }
        export_events(entry, recorder.events());
        solver_run_data[solver_name].push_back(std::move(entry));
    };

    // MES history: each entry is a new, better MES found
    for(const auto& entry : mes_history_data) {
        std::vector<Vertex> mes = lit::internalize(entry.at("mes").get<std::vector<ExternalVertex>>());
        LNSSubproblem subproblem{
            subproblem_input.uncovered_vertices,
            mes,
            subproblem_input.covering_assignments,
            static_cast<Lit>(problem_input.initial_phase.inf_map.get_n_concrete())
        };

        auto run_solver_type = [&] (const char* name, auto solver_nullptr) {
            using SolverType = std::remove_pointer_t<decltype(solver_nullptr)>;
            EventRecorder recorder;
            auto subproblem_and_result = make_and_run_common_interface_solver<SolverType>(
                std::move(subproblem),
                &problem_input.clauses,
                &recorder, 
                solve_time_limit
            );
            subproblem = std::move(subproblem_and_result.first);
            std::cerr << "\t - " << subproblem_and_result.second.outcome << std::endl;
            add_solver_outcome(name, std::move(subproblem_and_result.second), recorder);
        };

        auto run_cnpsatdsatur = [&] (const char* name, auto solver_nullptr) {
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
                false
            };
            EventRecorder recorder;
            solver.set_event_recorder(&recorder);
            auto after_make = std::chrono::steady_clock::now();
            auto res = timeout_wrapped_solve(solver, solve_time_limit);
            auto after_solve = std::chrono::steady_clock::now();
            std::vector<DynamicBitset> solution;
            if(res == SolverType::SolveResult::IMPROVED_SOLUTION) {
                solution = solver.get_partial_solution().assignments();
            }
            SolveResult solve_result{
                (res == SolverType::SolveResult::ABORTED ? "TIMEOUT" :
                 res == SolverType::SolveResult::SOLUTION_WAS_OPTIMAL ? "SOLUTION_WAS_OPTIMAL" :
                       "IMPROVED_SOLUTION"),
                std::move(solution),
                seconds_between(before_make, after_make),
                seconds_between(after_make, after_solve)
            };
            add_solver_outcome(name, std::move(solve_result), recorder);
        };
        
        for(std::size_t i = 0; i < 3; ++i) {
            // non-incremental SAT solver backed
            run_solver_type("fixed_sat[kissat]", static_cast<FixedMESSATImprovementSolver<KissatSolver>*>(nullptr));
            run_solver_type("fixed_sat[cadical]", static_cast<FixedMESSATImprovementSolver<CadicalSolver>*>(nullptr));
            run_solver_type("fixed_sat[cryptominisat]", static_cast<FixedMESSATImprovementSolver<CMSAT5Solver>*>(nullptr));
            run_solver_type("fixed_sat[lingeling]", static_cast<FixedMESSATImprovementSolver<LingelingSolver>*>(nullptr));
            
            // incremental SAT solver backed
            run_solver_type("incremental_sat[cadical]", static_cast<FixedMESIncrementalSATImprovementSolver<CadicalSolver>*>(nullptr));
            run_solver_type("incremental_sat[cryptominisat]", static_cast<FixedMESIncrementalSATImprovementSolver<CMSAT5Solver>*>(nullptr));
            run_solver_type("incremental_sat[lingeling]", static_cast<FixedMESIncrementalSATImprovementSolver<LingelingSolver>*>(nullptr));
        
            // SATDSatur (different interface)
            run_cnpsatdsatur("satdsatur[cadical]", static_cast<CadicalSolver*>(nullptr));
            run_cnpsatdsatur("satdsatur[cryptominisat]", static_cast<CMSAT5Solver*>(nullptr));
            run_cnpsatdsatur("satdsatur[lingeling]", static_cast<LingelingSolver*>(nullptr));
        }
        overall_output.at("solver_runs").push_back(std::move(solver_run_data));
        solver_run_data = OutputObject{};
    }

    if(!output_file_path.empty()) {
        std::ofstream output;
        output.exceptions(std::ios::failbit | std::ios::badbit);
        output.open(output_file_path, std::ios::out | std::ios::trunc);
        output << overall_output.dump(2);
    }
    return 0;
}
