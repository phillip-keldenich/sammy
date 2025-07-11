#include <sammy/cadical_solver.h>
#include <sammy/io.h>
#include <thread>
#include <mutex>
#include <limits>
#include <atomic>

using namespace sammy;

/**
 * Instance data and checks of feasibility, counting of interactions, etc.
 */
struct InstanceData {
    InstanceData(const nlohmann::json& json)
        : clauses(json.at("cnf_clauses")
                      .get<std::vector<std::vector<std::int32_t>>>()),
          num_variables(json.at("num_variables").get<std::int32_t>()),
          num_concrete(json.at("num_concrete_features").get<std::int32_t>()) {}

    std::optional<std::string>
    check_assignment(const std::vector<bool>& assignment) const {
        for (const auto& clause : clauses) {
            bool satisfied = false;
            for (auto lit : clause) {
                if (lit < 0) {
                    lit = -lit - 1;
                    if (!assignment[lit]) {
                        satisfied = true;
                        break;
                    }
                } else {
                    lit = lit - 1;
                    if (assignment[lit]) {
                        satisfied = true;
                        break;
                    }
                }
            }
            if (!satisfied) {
                return "Clause not satisfied!";
            }
        }
        return std::nullopt;
    }

    template <typename Assignments>
    std::optional<std::string>
    check_assignments(const Assignments& assignments) const {
        for (const auto& assignment : assignments) {
            auto result = check_assignment(assignment);
            if (result) {
                return result;
            }
        }
        return std::nullopt;
    }

    using Interaction = std::pair<std::int32_t, std::int32_t>;
    using MES = std::vector<Interaction>;
    using InteractionSet = HashSet<Interaction>;

    template<typename Assignments> 
    InteractionSet get_interaction_set(const Assignments& asmts) const {
        std::size_t num_concrete = std::size_t(this->num_concrete);
        std::size_t num_threads = std::thread::hardware_concurrency();
        std::atomic<std::size_t> current_var1(0);
        std::vector<std::vector<Interaction>> thread_interactions(num_threads);
        std::vector<std::thread> threads(num_threads);
        for(std::size_t i = 0; i < num_threads; ++i) {
            threads[i] = std::thread([&] (std::size_t thread_index) {
                auto& interactions = thread_interactions[thread_index];
                std::vector<bool> positive_partners(2 * num_concrete, false);
                std::vector<bool> negative_partners(2 * num_concrete, false);
                for(;;) {
                    std::size_t var1 = current_var1++;
                    if(var1 >= num_concrete) {
                        return; // no more variables to process
                    }
                    positive_partners.assign(2 * num_concrete, false);
                    negative_partners.assign(2 * num_concrete, false);
                    for(const auto& assignment : asmts) {
                        auto& out_partners = assignment[var1] ? 
                            positive_partners : negative_partners;
                        for(std::size_t var2 = var1 + 1; var2 < num_concrete; ++var2) {
                            std::size_t op_index = 2 * var2 + 1 - assignment[var2];
                            out_partners[op_index] = true;
                        }
                    }
                    std::int32_t plit1 = std::int32_t(var1 + 1);
                    std::int32_t nlit1 = -plit1;
                    for(std::size_t var2 = var1 + 1; var2 < num_concrete; ++var2) {
                        std::int32_t plit2 = std::int32_t(var2 + 1);
                        std::int32_t nlit2 = -plit2;
                        if(positive_partners[2 * var2]) {
                            interactions.emplace_back(plit1, plit2);
                        }
                        if(positive_partners[2 * var2 + 1]) {
                            interactions.emplace_back(nlit2, plit1);
                        }
                        if(negative_partners[2 * var2]) {
                            interactions.emplace_back(nlit1, plit2);
                        }
                        if(negative_partners[2 * var2 + 1]) {
                            interactions.emplace_back(nlit2, nlit1);
                        }
                    }
                }
            }, i);
        }
        std::for_each(threads.begin(), threads.end(), [] (auto& t) {t.join();});
        std::size_t num_interactions = 0;
        for(const auto& interactions : thread_interactions) {
            num_interactions += interactions.size();
        }
        InteractionSet interaction_set;
        interaction_set.reserve(num_interactions);
        for(const auto& interactions : thread_interactions) {
            for(const auto interaction : interactions) {
                interaction_set.emplace(interaction);
            }
        }
        return interaction_set;
    }

    template <typename Assignments>
    std::size_t count_covered_interactions(const Assignments& asmts,
                                           const MES& mes) const {
        InteractionSet covered_interactions = get_interaction_set(asmts);
        for (auto v : mes) {
            v = std::make_pair((std::min)(v.first, v.second),
                               (std::max)(v.first, v.second));
            if (!covered_interactions.count(v)) {
                throw std::runtime_error(
                    "Mutually exclusive set contains uncovered interaction: " +
                    std::to_string(v.first) + ", " + std::to_string(v.second));
            }
        }
        return covered_interactions.size();
    }

    std::vector<std::vector<std::int32_t>> clauses;
    std::int32_t num_variables;
    std::int32_t num_concrete;
};

/**
 * Class to check the mutual exclusivity of a set of interactions.
 */
class MESExclusivityChecker {
  public:
    using MES = std::vector<std::pair<std::int32_t, std::int32_t>>;

    MESExclusivityChecker(const InstanceData& instance, const MES& mes)
        : instance(&instance), solver(), mes(mes),
          locally_excluded(mes.size(), false) 
    {
        p_create_main_vars();
        p_add_main_clauses();
        p_add_mes_clauses();
        // incremental mode needed for Automotive02 checks
        //p_add_two_interaction_clauses(); // for non-incremental mode
    }

    /**
     * Check the mutual exclusivity of the interactions,
     * using parallel threads and the incremental mode.
     */
    void check_exclusivity() {
        std::size_t nthreads = std::thread::hardware_concurrency() / 2;
        if(nthreads <= 2) nthreads = 2;

        std::vector<std::unique_ptr<MESExclusivityChecker>> mes_checkers;
        std::mutex glob_excl_mutex;
        std::vector<bool> globally_excluded(mes.size(), false);
        std::size_t current_unchecked = 0;
        std::optional<std::string> error;
        mes_checkers.reserve(nthreads);
        for(std::size_t i = 0; i < nthreads; ++i) {
            mes_checkers.emplace_back(std::make_unique<MESExclusivityChecker>(*instance, mes));
        }
        std::vector<std::thread> threads(nthreads);
        for(std::size_t i = 0; i < nthreads; ++i) {
            threads[i] = std::thread([&] (std::size_t thread_index) {
                auto& mes_checker = *mes_checkers[thread_index];
                for(;;) {
                    std::size_t our_check;
                    {
                        std::unique_lock<std::mutex> lock(glob_excl_mutex);
                        if(error || current_unchecked >= mes.size()) {
                            return; // stop if an error has occurred
                        }
                        our_check = current_unchecked++;
                        mes_checker.exclude_interactions(globally_excluded);
                    }
                    auto cres = mes_checker
                            .check_exclusivity_of_interaction(our_check);
                    if(!cres) {
                        // we should have an error stored
                        return;
                    }
                    {
                        std::unique_lock<std::mutex> lock(glob_excl_mutex);
                        if(cres->empty()) {
                            // no error
                            globally_excluded[our_check] = true;
                        } else {
                            // error: store it and stop everyone
                            if(!error) {
                                error = *cres;
                                for(std::size_t i = 0; i < mes_checkers.size(); ++i) {
                                    mes_checkers[i]->solver.terminate();
                                }
                            }
                            return;
                        }
                    }
                }
            }, i);
        }
        std::for_each(threads.begin(), threads.end(),
                      [] (std::thread& t) { t.join(); });
        if(error) {
            throw std::runtime_error(*error);
        }
    }

    std::optional<std::string> 
        check_exclusivity_of_interaction(std::size_t interaction_index) 
    {
        std::vector<CadicalSolver::Lit> assumption;
        assumption.push_back(mes_vars[interaction_index]);
        std::vector<CadicalSolver::Lit> clause;
        for(std::size_t i = 0; i < mes.size(); ++i) {
            if(i == interaction_index) {
                clause.push_back(-mes_vars[i]);
            } else {
                clause.push_back(mes_vars[i]);
            }
        }
        solver.add_clause(clause.begin(), clause.end());
        auto result = solver.solve(assumption);
        if(!result) {
            return std::nullopt; // timeout means another solver found an error
        }
        if(!*result) {
            return std::string{}; // no error, interaction is exclusive
        }
        std::vector<std::pair<std::int32_t, std::int32_t>> simultaneous;
        auto model = solver.get_model();
        for (std::size_t i = 0; i < mes.size(); ++i) {
            if (model[mes_vars[i]]) {
                simultaneous.push_back(mes[i]);
            }
        }
        std::ostringstream msg;
        msg << "Mutually exclusive set contains simultaneously feasible "
               "interactions: ";
        for (const auto& interaction : simultaneous) {
            msg << "(" << interaction.first << ", " << interaction.second
                << ") ";
        }
        msg << std::endl;
        return msg.str();
    }

    void exclude_interactions(const std::vector<bool>& excludable) {
        for(std::size_t i = 0; i < mes.size(); ++i) {
            if(excludable[i]) {
                 if(!locally_excluded[i]) {
                    locally_excluded[i] = true;
                    solver.fix(-mes_vars[i]);
                 }
            }
        }
    }

    /* // nonincremental version
    void check_exclusivity() {
        auto result = solver.solve();
        if (!result) {
            throw std::runtime_error("Solver timed out?");
        }
        if (!*result) {
            return;
        }
        auto model = solver.get_model();
        MES simultaneous;
        for (std::size_t i = 0; i < mes.size(); ++i) {
            if (model[mes_vars[i]]) {
                simultaneous.push_back(mes[i]);
            }
        }
        std::ostringstream msg;
        msg << "Mutually exclusive set contains simultaneously feasible "
               "interactions: ";
        for (const auto& interaction : simultaneous) {
            msg << "(" << interaction.first << ", " << interaction.second
                << ") ";
        }
        throw std::runtime_error(msg.str());
    }*/

  private:
    using SatVar = CadicalSolver::Lit;

    /**
     * Add clauses to enforce that two of the MES interactions
     * are simultaneously true; non-incremental version.
     */
    void p_add_two_interaction_clauses() {
        for (std::size_t mes_index = 0; mes_index < mes.size(); ++mes_index) {
            for (std::size_t mi = 0; mi != mes.size(); ++mi) {
                if (mi == mes_index)
                    continue;
                solver.add_literal(mes_vars[mi]);
            }
            solver.finish_clause();
        }
    }

    /**
     * Create main variables for each variable in the instance.
     */
    void p_create_main_vars() {
        for (std::int32_t i = 0; i < instance->num_variables; ++i) {
            main_vars.push_back(solver.new_var());
        }
    }

    /**
     * Add clauses for the main variables.
     */
    void p_add_main_clauses() {
        for (const auto& clause : instance->clauses) {
            for (const auto lit : clause) {
                SatVar var =
                    (lit < 0) ? -main_vars[-(lit + 1)] : main_vars[lit - 1];
                solver.add_literal(var);
            }
            solver.finish_clause();
        }
    }

    /**
     * Create mes variables and ensure that
     * mes_vars[i] is true iff the interaction
     * is true in main_vars.
     */
    void p_add_mes_clauses() {
        for (const auto& interaction : mes) {
            mes_vars.push_back(solver.new_var());
            SatVar mes_var = mes_vars.back();
            SatVar lit1, lit2;
            if (interaction.first < 0) {
                lit1 = -main_vars[-(interaction.first + 1)];
            } else {
                lit1 = main_vars[interaction.first - 1];
            }
            if (interaction.second < 0) {
                lit2 = -main_vars[-(interaction.second + 1)];
            } else {
                lit2 = main_vars[interaction.second - 1];
            }
            solver.add_short_clause(-mes_var, lit1);
            solver.add_short_clause(-mes_var, lit2);
            solver.add_short_clause(mes_var, -lit1, -lit2);
        }
    }

    const InstanceData* instance;
    CadicalSolver solver;
    std::vector<SatVar> main_vars;
    std::vector<SatVar> mes_vars;
    MES mes;
    std::vector<bool> locally_excluded; // for incremental mode
};


int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <solution_file> <instance_file>"
                  << std::endl;
        return EXIT_FAILURE;
    }

    // read solution and instance
    InstanceData instance(read_json_path(argv[2]));
    using SolType = std::vector<std::vector<bool>>;
    using MesType = std::vector<std::pair<std::int32_t, std::int32_t>>;
    SolType solution;
    MesType mes;

    {
        auto solution_data = read_json_path(argv[1]);
        solution = solution_data.at("best_solution").get<SolType>();
        mes = solution_data.at("mutually_exclusive_set").get<MesType>();
        
        // validate the optimality flag
        if (solution_data.at("optimal").get<bool>()) {
            if(solution_data.at("lb").get<std::size_t>() !=
            solution_data.at("ub").get<std::size_t>()) 
            {
                std::cerr << "Error: Solution is marked as optimal, but lower "
                            "bound does not match upper bound!"
                        << std::endl;
                return EXIT_FAILURE;
            }
        } else {
            if (solution_data.at("lb").get<std::size_t>() ==
                solution_data.at("ub").get<std::size_t>()) {
                std::cerr << "Error: Solution is not optimal, but lower bound "
                            "matches upper bound!"
                        << std::endl;
                return EXIT_FAILURE;
            }
        }

        // validate the bounds
        if(solution_data.at("lb").get<std::size_t>() >
        solution_data.at("ub").get<std::size_t>()) 
        {
            std::cerr << "Error: Lower bound is greater than upper bound!" 
                    << std::endl;
            return EXIT_FAILURE;
        }

        // validate the solution size against the upper bound
        if (solution.size() != solution_data.at("ub").get<std::size_t>()) {
            std::cerr << "Error: Solution size does "
                        "not match upper bound!"
                    << std::endl;
            return EXIT_FAILURE;
        }
    }

    // validate that mes entries are on concrete features
    for(auto mes_entry : mes) {
        std::int32_t var1 = (std::abs)(mes_entry.first);
        std::int32_t var2 = (std::abs)(mes_entry.second);
        if(var1 <= 0 || var1 > instance.num_concrete ||
           var2 <= 0 || var2 > instance.num_concrete) 
        {
            std::cerr << "Error: MES entry (" << mes_entry.first << ", "
                      << mes_entry.second
                      << ") is not on concrete features!" << std::endl;
            std::cerr << "Concrete features: " << instance.num_concrete
                      << std::endl;
            return EXIT_FAILURE;
        }
    }

    // validate the size of individual configurations
    auto check_config_size = [&](const std::vector<bool>& config) {
        return config.size() == std::size_t(instance.num_variables);
    };
    if (!std::all_of(solution.begin(), solution.end(), check_config_size)) {
        std::cerr << "Error: Not all configurations "
                     "have the correct size!"
                  << std::endl;
        return EXIT_FAILURE;
    }

    // check feasibility of the solution
    auto feas_error = instance.check_assignments(solution);
    if (feas_error) {
        std::cerr << "Error: Solution is not feasible! " << *feas_error
                  << std::endl;
        return EXIT_FAILURE;
    }

    // count covered interactions and check MES vertices are covered
    std::size_t covered;
    try {
        covered = instance.count_covered_interactions(solution, mes);
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    // check the mutual exclusivity
    MESExclusivityChecker mes_checker(instance, mes);
    try {
        mes_checker.check_exclusivity();
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    nlohmann::json output{
        {"covered_interactions", covered},
        {"exclusivity_checked", true},
        {"feasibility_checked", true},
        {"solution_size", solution.size()},
        {"mes_size", mes.size()}
    };
    std::cout << output.dump(4) << std::endl;
    return 0;
}
