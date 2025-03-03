#include <sammy/cadical_solver.h>
#include <cadical.hpp>
#include <future>
#include <thread>
#include <mutex>
#include <iostream>

namespace sammy {

struct CadicalSolver::Terminator : CaDiCaL::Terminator {
    std::atomic<bool> interrupt_flag;

    Terminator() noexcept : interrupt_flag(false) {}

    bool terminate() override {
        return interrupt_flag.load(std::memory_order_relaxed);
    }

    void set_terminate() noexcept {
        interrupt_flag.store(true, std::memory_order_relaxed);
    }

    void reset_terminate() noexcept {
        interrupt_flag.store(false, std::memory_order_relaxed);
    }
};

CadicalSolver::CadicalSolver() {
    m_solver = std::make_unique<CaDiCaL::Solver>();
    m_terminator = std::make_unique<Terminator>();
    m_num_vars = 0;
    m_solver->connect_terminator(m_terminator.get());
}

CadicalSolver::~CadicalSolver() {
    m_solver->disconnect_terminator();
}

void CadicalSolver::reset() {
    m_solver->disconnect_terminator();
    m_solver = std::make_unique<CaDiCaL::Solver>();
    m_solver->connect_terminator(m_terminator.get());
    m_num_vars = 0;
}

auto CadicalSolver::new_var(bool /*reusable*/) -> Lit {
    if(m_num_vars == 0) {
        m_num_vars = m_solver->vars();
    }
    m_num_vars += 1;
    return m_num_vars;
}

auto CadicalSolver::num_vars() const noexcept -> Lit {
    return m_solver->vars();
}

auto CadicalSolver::get_model() const -> ModelMap {
    ModelMap map;
    map.model_map = std::vector<bool>(m_solver->vars() + 1, false);
    for(Lit l = 1, nvars = num_vars(); l <= nvars; ++l) {
        if(m_solver->val(l) > 0) {
            map.model_map[l] = true;
        }
    }
    return map;
}

void CadicalSolver::terminate() {
    m_terminator->set_terminate();
}
    
void CadicalSolver::reset_terminate() {
    m_terminator->reset_terminate();
}

static void cdc_add_assumptions(CaDiCaL::Solver& solver, const std::vector<int>& assumptions) {
    for(int l : assumptions) {
        solver.assume(l);
    }
}

std::optional<bool> 
CadicalSolver::solve(const std::vector<Lit>& assumptions, double time_limit)
{
    if(time_limit <= 0) return std::nullopt;
    cdc_add_assumptions(*m_solver, assumptions);
    reset_terminate();
    if(!std::isfinite(time_limit)) {
        int solve_result = m_solver->solve();
        if(solve_result == 10) {
            return true;
        } else if(solve_result == 20) {
            return false;
        } else {
            return std::nullopt;
        }
    }
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

void CadicalSolver::add_literal(Lit l) {
    m_solver->add(l);
}

} // namespace sammy
