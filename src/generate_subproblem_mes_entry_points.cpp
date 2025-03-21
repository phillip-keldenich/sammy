#include "io.h"
#include <sammy/barrage.h>
#include <sammy/experiment_flags.h>
#include <sammy/output.h>
#include <sammy/subproblem_solver_with_mes.h>

using namespace sammy;
namespace po = boost::program_options;

class MESSolverObserver {
  public:
    explicit MESSolverObserver(std::size_t removed_configs,
                               std::vector<Vertex> initial_mes)
        : m_removed_configs{removed_configs},
          m_initial_mes{std::move(initial_mes)},
          m_start_time{std::chrono::steady_clock::now()} {}

    bool operator()(SubproblemMESSolver& solver, bool aborted,
                    bool proved_mes_optimal, bool failed) {
        double t =
            seconds_between(m_start_time, std::chrono::steady_clock::now());
        if (m_init_time < 0) {
            m_init_time = t;
        }
        if (aborted) {
            m_aborted = true;
            m_total_time = t;
            m_mes_outcome = "ABORTED";
            return false;
        }
        if (proved_mes_optimal) {
            m_proved_optimal = true;
            if (last_mes_size() == m_removed_configs) {
                m_established_improvement_impossible = true;
            }
            m_total_time = t;
            p_handle_regular_iteration(solver, t);
            m_mes_outcome = "MES_OPTIMAL";
            return false;
        }
        if (failed) {
            m_failed = true;
            m_total_time = t;
            m_mes_outcome = "FAILED";
            return false;
        }
        p_handle_regular_iteration(solver, t);
        return true;
    }

    std::size_t last_mes_size() const {
        if (m_mes_history.empty()) {
            return m_initial_mes.size();
        }
        return m_mes_history.back().mes.size();
    }

    OutputObject get_data() const {
        OutputObject iteration_data;
        for (const auto& e : m_relaxation_history) {
            iteration_data.push_back(e.export_data());
        }
        OutputObject mes_data;
        for (const auto& e : m_mes_history) {
            mes_data.push_back(e.export_data());
        }
        return OutputObject{{"iterations", std::move(iteration_data)},
                            {"mes_history", std::move(mes_data)},
                            {"removed_configs", m_removed_configs},
                            {"initial_mes", lit::externalize(m_initial_mes)},
                            {"initial_time", m_init_time},
                            {"total_time", m_total_time},
                            {"aborted", m_aborted},
                            {"proved_optimal", m_proved_optimal},
                            {"failed", m_failed},
                            {"established_improvement_impossible",
                             m_established_improvement_impossible},
                            {"mes_outcome", m_mes_outcome}};
    }

  private:
    void p_handle_regular_iteration(SubproblemMESSolver& solver, double time) {
        ++m_iteration;
        const auto& mes = solver.mes_vertices();
        std::size_t mes_size = mes.size();
        std::size_t last_mes = last_mes_size();
        if (mes_size > last_mes) {
            m_mes_history.push_back({m_iteration, mes, time});
        }
        m_relaxation_history.push_back(
            {m_iteration, time, last_mes, mes_size,
             solver.last_change_was_cutting(), solver.last_change_was_pricing(),
             solver.num_lp_variables(), solver.num_lp_constraints(),
             solver.num_graph_vertices(), solver.num_lp_solve_calls(),
             solver.num_lp_nonzeros(), solver.get_relaxation_value(),
             solver.get_subgraph_ub(), solver.fractional_mes_support_size()});
    }

    struct RelaxationHistoryEntry {
        std::size_t iteration;
        double time_after_iteration;
        std::size_t mes_size_before_iteration;
        std::size_t mes_size_after_iteration;
        bool change_was_cut;
        bool change_was_price;
        std::size_t num_variables;
        std::size_t num_constraints;
        std::size_t num_graph_vertices;
        std::size_t total_solve_calls;
        std::size_t total_nonzeros;
        double relaxation_value_after;
        std::size_t subgraph_ub;
        std::size_t fractional_mes_support_size;

        OutputObject export_data() const {
            return OutputObject{
                {"iteration", iteration},
                {"time_after_iteration", time_after_iteration},
                {"mes_size_before_iteration", mes_size_before_iteration},
                {"mes_size_after_iteration", mes_size_after_iteration},
                {"change_was_cut", change_was_cut},
                {"change_was_price", change_was_price},
                {"num_variables", num_variables},
                {"num_constraints", num_constraints},
                {"num_graph_vertices", num_graph_vertices},
                {"total_solve_calls", total_solve_calls},
                {"total_nonzeros", total_nonzeros},
                {"relaxation_value_after", relaxation_value_after},
                {"subgraph_ub", subgraph_ub},
                {"fractional_mes_support_size", fractional_mes_support_size}};
        }
    };
    struct MESHistoryEntry {
        std::size_t iteration;
        std::vector<Vertex> mes;
        double discovery_time;

        OutputObject export_data() const {
            return OutputObject{{"iteration", iteration},
                                {"mes", lit::externalize(mes)},
                                {"discovery_time", discovery_time}};
        }
    };

    std::size_t m_removed_configs;
    std::vector<RelaxationHistoryEntry> m_relaxation_history;
    std::vector<MESHistoryEntry> m_mes_history;
    std::vector<Vertex> m_initial_mes;
    std::chrono::steady_clock::time_point m_start_time;
    std::string m_mes_outcome;
    double m_init_time{-1};
    double m_total_time{-1};
    std::size_t m_iteration = std::numeric_limits<std::size_t>::max();
    bool m_aborted{false};
    bool m_proved_optimal{false};
    bool m_failed{false};
    bool m_established_improvement_impossible{false};
};

/**
 * Generate, for a given subproblem,
 * a trace of the MES solver, including
 * the MES sizes and changes to it,
 * the behavior of the relaxation, and so on.
 */
int main(int argc, char** argv) {
    std::string formula_file;
    std::string subproblem_file;
    std::string output_file;
    double time_limit = 3600.0;

    boost::program_options::options_description desc("Barrage options");
    desc.add_options()("formula", required_value(formula_file),
                       "Path to the formula file.")(
        "subproblem", required_value(subproblem_file),
        "Path to the subproblem file.")("output", required_value(output_file),
                                        "Path to the output file.")(
        "time-limit", value_with_default(time_limit),
        "Time limit for the solver.")("help", "Print this help message.");

    po::variables_map vmap;
    auto cmdline_parser = boost::program_options::command_line_parser(
        argc, const_cast<const char* const*>(argv));
    boost::program_options::store(cmdline_parser.options(desc).run(), vmap);
    vmap.notify();
    if (vmap.count("help")) {
        std::cerr << desc << std::endl;
        return 0;
    }

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
        portfolio.limited_reduce_universe((std::min)(10.0, 0.02 * time_limit));
    }
    LNSSubproblem subproblem{
        std::move(subproblem_input.uncovered_vertices),
        std::move(subproblem_input.best_local_mes),
        std::move(subproblem_input.covering_assignments),
        static_cast<Lit>(portfolio.get_infeasibility_map().get_n_concrete())};
    portfolio.implied_cache().remove_implied(subproblem.uncovered_universe);
    portfolio.implied_cache().replace_implied(
        subproblem.mutually_exclusive_set);
    MESSolverObserver observer(subproblem.removed_configurations.size(),
                               subproblem.mutually_exclusive_set);
    SubproblemMESSolver subproblem_mes_solver{
        &portfolio, std::move(subproblem),
        SharedDBPropagator(&portfolio.get_clauses()), &recorder, 0};
    auto time_before = std::chrono::steady_clock::now();
    auto deadline = time_before + std::chrono::duration<double>(time_limit);
    std::mutex m;
    std::condition_variable cv;
    bool done = false;
    auto timeout_watcher_routine = [&]() {
        std::unique_lock l{m};
        if (!cv.wait_until(l, deadline, [&]() { return done; })) {
            subproblem_mes_solver.abort();
        }
    };
    std::thread timeout_watcher(timeout_watcher_routine);
    subproblem_mes_solver.extended_solve(observer);
    {
        std::unique_lock l{m};
        done = true;
        cv.notify_one();
    }
    timeout_watcher.join();

    OutputObject overall_output_data = observer.get_data();
    overall_output_data["formula_file"] = formula_file;
    overall_output_data["subproblem_file"] = subproblem_file;
    overall_output_data["problem_input"] = std::move(problem_input.raw_input);
    overall_output_data["subproblem_input"] =
        std::move(subproblem_input.raw_input);
    write_json_path(output_file, overall_output_data);
    return 0;
}
