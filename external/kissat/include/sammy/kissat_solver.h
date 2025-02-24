#ifndef SAMMY_KISSAT_SOLVER_H_INCLUDED_
#define SAMMY_KISSAT_SOLVER_H_INCLUDED_

#include <algorithm>
#include <vector>
#include <optional>
#include <thread>
#include <limits>
#include <future>
#include <cmath>
#include <climits>
#include <cfloat>

namespace sammy {

/*
BaseSolver:
    BaseSolver()
    type Lit (supports operator-, operator==/!=; may be a primitive type such as int)
    Lit new_var()
    void add_short_clause(Lit...)
    void add_clause(LitIterator, LitIterator)
    void add_literals(Lit...)
    void add_literal(Lit)
    void finish_clause()
    std::optional<bool> solve(time_limit) // must only be called once!
    SomeMapType(Lit->bool) get_model()
    void terminate()
*/

/**
 * C++ binding for kissat to keep all
 * kissat headers from seeping into
 * the remaining project.
 */
class KissatSolver {
  public:
    using Lit = int;

    class ModelMap {
      public:
        bool operator[](Lit l) const {
            bool r = model[std::abs(l) - 1];
            return l < 0 ? !r : r;
        }

        const std::vector<bool>& raw() const {
            return model;
        }
        
        std::vector<bool>& raw() {
            return model;
        }

      private:
        std::vector<bool> model;

        explicit ModelMap(std::vector<bool> model) noexcept :
            model(std::move(model))
        {}

        friend class KissatSolver;
    };

    KissatSolver();
    ~KissatSolver();

    /**
     * Prepare the model so it has at least num_vars variables.
     */
    void reserve(Lit num_vars);

    /**
     * Solve without time limit (but still might
     * be aborted or hit another limit).
     * Returns nullopt if aborted, {true} if sat and
     * {false} if unsat.
     */
    std::optional<bool> solve();

    /**
     * Solve with time limit in seconds (but still might
     * be aborted earlier or hit another limit).
     * Returns nullopt if aborted, {true} if sat and
     * {false} if unsat.
     */
    std::optional<bool> solve(double time_limit) {
        if(time_limit <= 0) return std::nullopt;
        if(!std::isfinite(time_limit)) return solve();
        std::promise<std::optional<bool>> res_promise;
        std::future<std::optional<bool>> result = res_promise.get_future();
        auto smain = [&] () { res_promise.set_value(solve()); };
        std::thread t{smain};
        auto status = result.wait_for(std::chrono::duration<double>(time_limit));
        if(status == std::future_status::timeout) {
            this->terminate(); // abort solution process on timeout
        }
        std::optional<bool> r = result.get();
        t.join();
        return r;
    }

    /**
     * Add a clause from a sequence of literals.
     */
    template<typename InputIterator,
             std::enable_if_t<!std::is_integral_v<InputIterator>,int> = 0>
    void add_clause(InputIterator begin, InputIterator end)
    {
        std::for_each(begin, end, [&] (Lit l) {add_literal(l);});
        finish_clause();
    }

    /**
     * Add a small clause from a given variadic number of literals.
     */
    template<typename... Args>
    void add_short_clause(Lit l1, Args&&... lits) {
        add_literal(l1);
        add_short_clause(std::forward<Args>(lits)...);
    }

    /**
     * Add a small clause from a given variadic number of literals.
     */
    void add_short_clause() {
        finish_clause();
    }

    /**
     * Add a variadic number of literals to the current clause.
     */
    template<typename... Args>
    void add_literals(Lit l1, Args&&... args) {
        add_literal(l1);
        add_literals(std::forward<Args>(args)...);
    }
    void add_literals() noexcept {}

    /**
     * Add a single literal to the current clause.
     */
    void add_literal(Lit lit);

    /**
     * Finish the current clause.
     */
    void finish_clause();
    
    /**
     * Create a new variable.
     */
    Lit new_var() noexcept {
        return ++m_num_vars;
    }

    /**
     * Get the number of variables.
     */
    Lit num_vars() const noexcept {
        return m_num_vars;
    }

    /**
     * Set limits on the number of conflicts/decisions.
     */
    void set_conflict_limit(unsigned limit);
    void set_decision_limit(unsigned limit);

    /**
     * Abort solution; usually, will result in a
     * nullopt result from solve.
     */
    void terminate();

    /**
     * Get the model.
     * Throws an exception if not already solved.
     * Undefined behavior if solve() returned nullopt or false.
     */
    ModelMap get_model() const;

    static const char* name() noexcept {
        return "kissat";
    }

  private:
    // solver handle (made opaque as void pointer)
    void *solver;
    // number of variables
    Lit m_num_vars;
    // take note if this is already solved
    bool already_solved = false;
};

}

#endif
