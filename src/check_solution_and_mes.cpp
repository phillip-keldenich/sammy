#include <sammy/cadical_solver.h>
#include <sammy/io.h>

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

    using MES = std::vector<std::pair<std::int32_t, std::int32_t>>;

    template <typename Assignments>
    std::size_t count_covered_interactions(const Assignments& asmts,
                                           const MES& mes) const {
        HashSet<std::pair<std::int32_t, std::int32_t>> covered_interactions;
        for (const auto& assignment : asmts) {
            for (std::int32_t var = 1; var <= num_concrete; ++var) {
                std::int32_t lit1 = assignment[var - 1] ? var : -var;
                for (std::int32_t var2 = var + 1; var2 <= num_concrete; ++var2)
                {
                    std::int32_t lit2 = assignment[var2 - 1] ? var2 : -var2;
                    std::int32_t l1 = (std::min)(lit1, lit2);
                    std::int32_t l2 = (std::max)(lit1, lit2);
                    covered_interactions.emplace(l1, l2);
                }
            }
        }
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
        : instance(&instance), solver(), mes(mes) {
        p_create_main_vars();
        p_add_main_clauses();
        p_add_mes_clauses();
        p_add_two_interaction_clauses();
    }

    /**
     * Check the mutual exclusivity of the interactions.
     */
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
    }

  private:
    using SatVar = CadicalSolver::Lit;

    /**
     * Add clauses to enforce that two of the MES interactions
     * are simultaneously true.
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
