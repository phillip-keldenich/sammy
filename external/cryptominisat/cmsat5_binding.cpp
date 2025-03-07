#include <sammy/cmsat5_solver.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include "cryptominisat.h"

namespace sammy {

class CMSAT5Solver::Impl {
  public:
    Impl() : m_solver(nullptr, &m_interrupt_flag) {
        m_solver.set_num_threads(1);
    }

    void add_clause(const std::vector<Lit>& clause) {
        m_translated_buffer.clear();
        std::transform(clause.begin(), clause.end(), std::back_inserter(m_translated_buffer), 
                       [] (Lit l) { return CMSat::Lit::toLit(l.m_lit); });
        m_solver.add_clause(m_translated_buffer);
    }

    std::size_t num_vars() const noexcept {
        return m_solver.nVars();
    }

    Lit new_var() {
        std::uint32_t var_index = m_solver.nVars();
        m_solver.new_var();
        return Lit(CMSat::Lit(var_index, false).toInt());
    }

    Lit new_vars(std::size_t new_vars) {
        std::uint32_t var_index = m_solver.nVars();
        m_solver.new_vars(new_vars);
        return Lit(CMSat::Lit(var_index, false).toInt());
    }
    
    std::optional<bool> solve(const std::vector<CMSat::Lit>& assumptions,
                              double time_limit)
    {
        if(m_interrupt_flag.load()) {
            m_interrupt_flag.store(false);
            return std::nullopt;
        }
        if(!std::isfinite(time_limit)) {
            return p_solve(assumptions);
        }
        std::mutex mtx;
        std::condition_variable cond;
        bool done = false;
        std::thread watcher([&] () {
            std::unique_lock l{mtx};
            std::chrono::duration<double> tlim(time_limit);
            if(!cond.wait_for(l, tlim, [&] () { return done; })) {
                m_interrupt_flag.store(true);
            }
        });
        auto result = p_solve(assumptions);
        {
            std::unique_lock l{mtx};
            done = true;
        }
        cond.notify_one();
        watcher.join();
        return result;
    }

    ModelMap get_model() const {
        const std::vector<CMSat::lbool>& cmm = m_solver.get_model();
        const std::size_t n = m_solver.nVars();
        std::vector<bool> result(n, false);
        for(std::size_t i = 0; i < n; ++i) {
            if(cmm[i] == CMSat::l_True) {
                result[i] = true;
            }
        }
        return ModelMap(std::move(result));
    }

    void terminate() {
        m_interrupt_flag.store(true);
    }

    void reset_terminate() {
        m_interrupt_flag.store(false);
    }

    bool should_terminate() const noexcept {
        return m_interrupt_flag.load();
    }

  private:
    std::optional<bool> p_solve(const std::vector<CMSat::Lit>& assumptions) {
        const std::vector<CMSat::Lit> *assumptions_ptr = assumptions.empty() ? nullptr : &assumptions;
        CMSat::lbool res = m_solver.solve(assumptions_ptr);
        if(res == CMSat::l_Undef) {
            return std::nullopt;
        }
        m_interrupt_flag.store(false);
        return res == CMSat::l_True;
    }

    std::atomic<bool> m_interrupt_flag{false};
    CMSat::SATSolver m_solver;
    std::vector<CMSat::Lit> m_translated_buffer;
};

void CMSAT5Solver::finish_clause() {
    assert(m_impl);
    m_impl->add_clause(m_clause_buffer);
    m_clause_buffer.clear();
}

void CMSAT5Solver::reset() {
    m_clause_buffer.clear();
    if(m_impl) {
        delete m_impl;
        m_impl = nullptr;
    }
    m_impl = new Impl();
}

CMSAT5Solver::CMSAT5Solver() :
    m_impl(new Impl())
{}

CMSAT5Solver::~CMSAT5Solver() {
    if(m_impl) {
        delete m_impl;
        m_impl = nullptr;
    }
}

std::size_t CMSAT5Solver::num_vars() const noexcept {
    assert(m_impl);
    return m_impl->num_vars();
}

auto CMSAT5Solver::new_var(bool /*reusable*/) -> Lit {
    assert(m_impl);
    return m_impl->new_var();
}

void CMSAT5Solver::fix(Lit l) {
    add_short_clause(l);
}

CMSAT5Solver::Lit CMSAT5Solver::Lit::operator-() const noexcept {
    CMSat::Lit lit = CMSat::Lit::toLit(m_lit);
    CMSat::Lit negated = ~lit;
    return Lit(negated.toInt());
}

std::optional<bool> 
CMSAT5Solver::solve(const std::vector<Lit>& assumptions,
                    double time_limit)
{
    if(time_limit <= 0) return std::nullopt;
    if(assumptions.empty()) {
        return m_impl->solve({}, time_limit);
    }
    std::vector<CMSat::Lit> translated_assumptions;
    std::transform(assumptions.begin(), assumptions.end(), std::back_inserter(translated_assumptions),
                   [] (Lit l) { return CMSat::Lit::toLit(l.m_lit); });
    return m_impl->solve(translated_assumptions, time_limit);
}

bool CMSAT5Solver::ModelMap::operator[](Lit l) const noexcept {
    CMSat::Lit lit = CMSat::Lit::toLit(l.m_lit);
    std::uint32_t var = lit.var();
    bool negated = lit.sign();
    return m_model[var] ? !negated : negated;
}

CMSAT5Solver::ModelMap CMSAT5Solver::get_model() const {
    assert(m_impl);
    return m_impl->get_model();
}

void CMSAT5Solver::terminate() {
    assert(m_impl);
    m_impl->terminate();
}

void CMSAT5Solver::reset_terminate() {
    assert(m_impl);
    m_impl->reset_terminate();
}

bool CMSAT5Solver::should_terminate() const noexcept {
    assert(m_impl);
    return m_impl->should_terminate();
}

CMSAT5Solver::Lit CMSAT5Solver::lit_from_dimacs_int(std::int32_t l) const {
    assert(l != 0);
    std::uint32_t var_index;
    if(l < 0) {
        var_index = -(l + 1);
    } else {
        var_index = l - 1;
    }
    return Lit(CMSat::Lit(var_index, l < 0).toInt());
}

auto CMSAT5Solver::new_vars(std::size_t num_vars) -> Lit {
    assert(m_impl);
    return m_impl->new_vars(num_vars);
}

}
