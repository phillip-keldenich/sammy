#ifndef SAMMY_CMSAT5_SOLVER_H_INCLUDED_
#define SAMMY_CMSAT5_SOLVER_H_INCLUDED_

#include <vector>
#include <optional>
#include <limits>
#include <atomic>
#include <algorithm>
#include <cmath>
#include <climits>
#include <cfloat>
#include <cassert>
#include <cstdint>

/**
 * Like all the other solvers, compiler-firewall the
 * internals of cryptominisat5 to prevent any header content
 * from cryptominisat seeping into the rest of the codebase.
 * The code in cmsat5_binding and cryptominisat itself are the only
 * places that see cryptominisat5 headers or any other internals.
 */

namespace sammy {

class CMSAT5Solver {
  private:
    class Impl;

  public:
    class Lit {
      public:
        Lit() noexcept = default;
        Lit operator-() const noexcept;
        bool operator==(Lit other) const noexcept { return m_lit == other.m_lit; }
        bool operator!=(Lit other) const noexcept { return m_lit != other.m_lit; }
        bool operator< (Lit other) const noexcept { return m_lit < other.m_lit; }
        bool operator> (Lit other) const noexcept { return m_lit > other.m_lit; }
        bool operator<=(Lit other) const noexcept { return m_lit <= other.m_lit; }
        bool operator>=(Lit other) const noexcept { return m_lit >= other.m_lit; }

      private:
        friend class CMSAT5Solver;
        friend class CMSAT5Solver::Impl;

        explicit Lit(std::uint32_t x) noexcept : m_lit(x) {}

        std::uint32_t m_lit;
    };

    class ModelMap {
      public:
        bool operator[](Lit l) const noexcept;

        std::vector<bool> &raw() noexcept { return m_model; }
        const std::vector<bool> &raw() const noexcept { return m_model; }

      private:
        friend class CMSAT5Solver;
        friend class CMSAT5Solver::Impl;

        ModelMap(std::vector<bool> model) noexcept :
            m_model(std::move(model))
        {}

        std::vector<bool> m_model;
    };

    CMSAT5Solver();
    ~CMSAT5Solver();

    // not copyable
    CMSAT5Solver(const CMSAT5Solver&) = delete;
    CMSAT5Solver& operator=(const CMSAT5Solver&) = delete;

    // nothrow movable
    CMSAT5Solver(CMSAT5Solver&& o) noexcept :
        m_impl(o.m_impl),
        m_clause_buffer(std::move(o.m_clause_buffer))
    {
        o.m_impl = nullptr;
    }
    CMSAT5Solver &operator=(CMSAT5Solver&& o) noexcept {
        std::swap(m_impl, o.m_impl);
        std::swap(m_clause_buffer, o.m_clause_buffer);
        return *this;
    }

    /**
     * Get the number of variables.
     */
    std::size_t num_vars() const noexcept;

    /**
     * With cryptominisat, all variables are reusable.
     */
    Lit new_var(bool reusable = true);

    /**
     * Create a number of new variables, returning the first one.
     */
    Lit new_vars(std::size_t num_vars);

    /**
     * Add a short clause to the solver, consisting of the given literals.
     */
    template<typename... Lits>
    void add_short_clause(Lits&&... lits) {
        p_add_to_clause(std::forward<Lits>(lits)...);
        finish_clause();
    }

    /**
     * Add literals to the clause that is currently being built.
     */
    template<typename... Lits>
    void add_literals(Lits&&... lits) {
        p_add_to_clause(std::forward<Lits>(lits)...);
    }

    /**
     * Add a single literal to the clause that is currently being built.
     */
    void add_literal(Lit l) {
        m_clause_buffer.push_back(l);
    }

    /**
     * Add the clause that is currently being built to the solver.
     */
    void finish_clause();

    /**
     * Add a clause from a range of Lits.
     */
    template<typename LitIterator>
    void add_clause(LitIterator begin, LitIterator end) {
        std::copy(begin, end, std::back_inserter(m_clause_buffer));
        finish_clause();
    }

    /**
     * Trigger the solver to run.
     */
    std::optional<bool> solve(const std::vector<Lit>& assumptions = {},
                              double time_limit = std::numeric_limits<double>::infinity());

    /**
     * Get the model; throws if the solver has not returned SAT (i.e., {true}).
     */
    ModelMap get_model() const;

    /**
     * Interrupt the solver.
     */
    void terminate();

    /**
     * Reset the termination flag.
     */
    void reset_terminate();

    /**
     * Get the termination flag status.
     */
    bool should_terminate() const noexcept;

    /**
     * Fix the given literal to true.
     */
    void fix(Lit l);

    /**
     * Reset the solver, removing all variables and clauses.
     */
    void reset();

    /**
     * Get the name of the solver.
     */
    static const char* name() noexcept {
        return "cryptominisat5";
    }

    /**
     * Translate a DIMACS integer to a Lit.
     */
    Lit lit_from_dimacs_int(std::int32_t l) const;

  private:
    template<typename... Lits>
    void p_add_to_clause(Lits&&... lits) {
        (add_literal(std::forward<Lits>(lits)), ...);
    }

    Impl* m_impl;
    std::vector<Lit> m_clause_buffer;
};

}

#endif
