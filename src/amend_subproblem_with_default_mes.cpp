#include "io.h"
#include <sammy/barrage.h>
#include <sammy/experiment_flags.h>
#include <sammy/output.h>
#include <sammy/subproblem_solver_with_mes.h>

using namespace sammy;
namespace po = boost::program_options;

int main(int argc, char** argv) {
    std::string formula_file;
    std::string subproblem_file;
    std::string output_file;

    /**
     * Parse command line arguments.
     */
    boost::program_options::options_description desc("Barrage options");
    desc.add_options()("universe-and-clauses", required_value(formula_file),
                       "Path to the universe-and-clauses file.")(
        "subproblem", required_value(subproblem_file),
        "Path to the subproblem file.")("output", required_value(output_file),
                                        "Path to the output file.")(
        "help", "Print this help message.");
    po::variables_map vmap;
    auto cmdline_parser = boost::program_options::command_line_parser(
        argc, const_cast<const char* const*>(argv));
    boost::program_options::store(cmdline_parser.options(desc).run(), vmap);
    vmap.notify();
    if (vmap.count("help")) {
        std::cerr << desc << std::endl;
        return 0;
    }

    /**
     * Parse input (problem and subproblem), and prepare the input.
     */
    ProblemInput problem_input = load_problem_input(formula_file);
    SubproblemInput subproblem_input = load_subproblem_input(subproblem_file);
    auto ticket = publish_clauses(problem_input.clauses);
    EventRecorder recorder;
    recorder.set_print_events(false);
    PortfolioSolver portfolio{ticket, &recorder,
                              std::move(problem_input.initial_phase)};
    if (portfolio.get_universe_size() < 2'000'000) {
        portfolio.reduce_universe();
    } else {
        portfolio.limited_reduce_universe(20.0);
    }

    /**
     * Construct a subproblem object from the input.
     */
    LNSSubproblem subproblem{
        std::move(subproblem_input.uncovered_vertices),
        std::move(subproblem_input.best_local_mes),
        std::move(subproblem_input.covering_assignments),
        static_cast<Lit>(portfolio.get_infeasibility_map().get_n_concrete())};
    const auto& ic = portfolio.implied_cache();
    ic.remove_implied(subproblem.uncovered_universe);
    ic.replace_implied(subproblem.mutually_exclusive_set);

    /**
     * Create solver and run it using the
     * default MES refinement settings.
     */
    SubproblemMESSolver subproblem_mes_solver{
        &portfolio, std::move(subproblem),
        SharedDBPropagator(&portfolio.get_clauses()), &recorder, 0};
    subproblem_mes_solver.solve();

    /**
     * Store the result.
     */
    nlohmann::json amendment;
    recorder.store_event("STORING_AMENDMENT");
    export_events(amendment, recorder);
    amendment["mes"] = lit::externalize(subproblem_mes_solver.mes_vertices());
    amendment["mes_size"] = amendment["mes"].size();
    amendment["problem_input"] = problem_input.raw_input;
    subproblem_input.raw_input["amendment"] = std::move(amendment);
    write_json_path(output_file, subproblem_input.raw_input);
    return 0;
}
