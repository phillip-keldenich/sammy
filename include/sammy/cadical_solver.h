#ifndef SAMMY_CADICAL_BINDINGS_H_INCLUDED_
#define SAMMY_CADICAL_BINDINGS_H_INCLUDED_

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <vector>

namespace CaDiCaL {

/**
 * Forward declaration of the actual solver type.
 */
class Solver;

} // namespace CaDiCaL

namespace sammy {

/**
 * Bindings to make CaDiCaL usable in our framework,
 * and to avoid symbols seeping into the rest of our
 * projects.
 */
class CadicalSolver {
  public:
    using Lit = int;

    /**
     * Add a new variable.
     * Reusable is ignored in this solver.
     */
    Lit new_var(bool reusable = true);

    /**
     * Get the number of variables.
     */
    Lit num_vars() const noexcept;

    /**
     * Add a short clause.
     */
    template <typename... Lits> void add_short_clause(Lits... lits) {
        add_literals(lits...);
        finish_clause();
    }

    /**
     * Add a variadic number of literals to the current clause.
     */
    template <typename... Lits> void add_literals(Lits... lits) {
        (add_literal(lits), ...);
    }

    /**
     * Finish the current clause.
     */
    void finish_clause() { add_literal(0); }

    /**
     * Add a clause from a sequence of literals.
     */
    template <typename LitIterator,
              std::enable_if_t<!std::is_integral_v<LitIterator>, int> = 0>
    void add_clause(LitIterator begin, LitIterator end) {
        std::for_each(begin, end, [this](Lit l) { add_literal(l); });
        finish_clause();
    }

    /**
     * Add a single literal to the current clause.
     */
    void add_literal(Lit l);

    /**
     * Create an empty solver.
     */
    CadicalSolver();

    /**
     * Destroy the solver.
     */
    ~CadicalSolver();

    /**
     * Reset the solver.
     */
    void reset();

    /**
     * Fix a variable to a given value.
     */
    void fix(Lit l) { add_short_clause(l); }

    /**
     * Asynchronously terminate the solver.
     */
    void terminate();

    /**
     * Reset the termination flag.
     */
    void reset_terminate();

    /**
     * Solve the current formula.
     */
    std::optional<bool>
    solve(const std::vector<Lit>& assumptions = {},
          double time_limit = std::numeric_limits<double>::infinity());

    /**
     * Store a model returned by the solver;
     * allows querying the truth values of literals.
     */
    class ModelMap {
      public:
        ModelMap() = default;

        bool operator[](Lit l) const noexcept {
            if (l < 0) {
                return !model_map[-l];
            } else {
                return model_map[l];
            }
        }

        std::vector<bool> model_map;
    };

    /**
     * After a successful solve, get the model.
     */
    ModelMap get_model() const;

    /**
     * Get the solver name.
     */
    static const char* name() noexcept { return "CaDiCaL"; }

  private:
    struct Terminator;
    std::unique_ptr<CaDiCaL::Solver> m_solver;
    std::unique_ptr<Terminator> m_terminator;
    Lit m_num_vars = 0;
};

} // namespace sammy

#endif
