#ifndef SAMMY_LINGELING_SOLVER_H_INCLUDED_
#define SAMMY_LINGELING_SOLVER_H_INCLUDED_

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cfloat>
#include <climits>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

namespace sammy {

/*
IncrementalSolver:
    IncrementalSolver()
    type Lit (supports operator-, operator==/!=; may be a primitive type such as
int) Lit new_var(bool reusable = true) -> add a new variable. if reusable=false,
cannot be used in new clauses/assumptions after the next solve is run. Lit
num_vars() const noexcept void add_short_clause(Lit...) void
add_literals(Lit...) void add_literal(Lit l) void finish_clause() void
add_clause(LitIterator begin, LitIterator end) SomeMapType(Lit->bool)
get_model() std::optional<bool> solve(const std::vector<Lit>& assumptions = {},
time_limit=infinity) void terminate() void reset_terminate() void fix(Lit l) //
fix a currently-reusable variable to make the given literal true;
                    // the variable becomes non-reusable in the process
    void reset() // reset the solver, removing all variables and clauses.
*/

class LingelingSolver {
  public:
    using Lit = int;

    class ModelMap {
      public:
        bool operator[](Lit l) const {
            bool r = model[std::abs(l) - 1];
            return l < 0 ? !r : r;
        }

        const std::vector<bool>& raw() const { return model; }

        std::vector<bool>& raw() { return model; }

      private:
        std::vector<bool> model;

        explicit ModelMap(std::vector<bool> model) noexcept
            : model(std::move(model)) {}

        friend class LingelingSolver;
    };

    LingelingSolver();
    LingelingSolver(const LingelingSolver&) = delete;
    LingelingSolver& operator=(const LingelingSolver&) = delete;

    LingelingSolver(LingelingSolver&& o) noexcept
        : solver(o.solver), m_num_vars(o.m_num_vars), m_terminate_flag(false) {
        o.solver = nullptr;
    }

    LingelingSolver& operator=(LingelingSolver&& o) noexcept {
        std::swap(solver, o.solver);
        std::swap(m_num_vars, o.m_num_vars);
        m_terminate_flag.store(o.m_terminate_flag.load());
        return *this;
    }

    ~LingelingSolver();

    Lit num_vars() const noexcept { return m_num_vars; }
    Lit new_var(bool reusable = true);

    template <typename Arg1, typename... Args>
    void add_literals(Arg1&& a1, Args&&... args) {
        add_literal(std::forward<Arg1>(a1));
        add_literals(std::forward<Args>(args)...);
    }
    void add_literals() noexcept {}
    void add_literal(Lit l);
    void finish_clause() { add_literal(0); }

    template <typename... Args> void add_short_clause(Args&&... args) {
        add_literals(std::forward<Args>(args)..., 0);
    }

    template <typename LitIterator>
    void add_clause(LitIterator begin, LitIterator end) {
        assert((std::none_of(begin, end, [](Lit l) { return l == 0; })));
        std::for_each(begin, end, [&](Lit l) { add_literal(l); });
        finish_clause();
    }

    ModelMap get_model() const;
    std::optional<bool>
    solve(const std::vector<Lit>& assumptions = {},
          double time_limit = std::numeric_limits<double>::infinity());
    void terminate();

    void fix(Lit l);

    bool should_terminate() const noexcept { return m_terminate_flag.load(); }

    void reset();

    static const char* name() noexcept { return "Lingeling"; }

  private:
    void* solver;
    Lit m_num_vars;
    std::atomic<bool> m_terminate_flag;
};

} // namespace sammy

#endif
