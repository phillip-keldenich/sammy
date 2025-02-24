#include <sammy/lingeling_solver.h>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <iostream>

extern "C" {
#include "lglib.h"
}

namespace sammy {

static int lingeling_terminate_callback(void* obj) {
    return static_cast<LingelingSolver*>(obj)->should_terminate();
}

static void lingeling_init(void*& solver, LingelingSolver* this_) {
    LGL* s = lglinit();
    if(!s) throw std::bad_alloc();
    lglseterm(s, &lingeling_terminate_callback, static_cast<void*>(this_));
    solver = static_cast<void*>(s);
}

LingelingSolver::LingelingSolver() :
    solver(nullptr),
    m_num_vars(0),
    m_terminate_flag(false)
{
    lingeling_init(solver, this);
}

void LingelingSolver::reset() {
    m_terminate_flag.store(false);
    m_num_vars = 0;
    if(solver) {
        lglrelease(static_cast<LGL*>(solver));
        solver = nullptr;
    }
    lingeling_init(solver, this);
}

LingelingSolver::~LingelingSolver() {
    if(solver) {
        lglrelease(static_cast<LGL*>(solver));
        solver = nullptr;
    }
}

auto LingelingSolver::new_var(bool reusable) -> Lit {
    LGL* s = static_cast<LGL*>(solver);
    Lit nv = ++m_num_vars;
    if(reusable) {
        lglfreeze(s, nv);
    }
    return nv;
}

void LingelingSolver::add_literal(Lit l) {
    m_num_vars = (std::max)((std::abs)(l), m_num_vars);
    lgladd(static_cast<LGL*>(solver), l);
}

void LingelingSolver::fix(Lit l) {
    LGL* s = static_cast<LGL*>(solver);
    lgladd(s, l);
    lgladd(s, 0);
    lglmelt(s, l);
}

void LingelingSolver::terminate() {
    m_terminate_flag.store(true);
}

std::optional<bool> LingelingSolver::solve(const std::vector<Lit>& assumptions, double time_limit) {
    if(time_limit <= 0 || should_terminate()) {
        m_terminate_flag.store(false);
        return std::nullopt;
    }
    LGL* s = static_cast<LGL*>(solver);
    auto to_result = [&] (int rcode) -> std::optional<bool> {
        m_terminate_flag.store(false);
        if(rcode == 10) return true;
        if(rcode == 20) return false;
        return std::nullopt;
    };
    std::for_each(assumptions.begin(), assumptions.end(), [&] (Lit l) {lglassume(s, l);});
    if(std::isfinite(time_limit)) {
        std::mutex mtx;
        std::condition_variable cond;
        bool done = false;
        std::thread watcher([&] () {
            std::unique_lock l{mtx};
            std::chrono::duration<double> tlim(time_limit);
            static_cast<void>(cond.wait_for(l, tlim, [&] () {return done;}));
            if(!done) {
                terminate();
            }
        });
        auto result = to_result(lglsat(s));
        {
            std::unique_lock l{mtx};
            done = true;
        }
        cond.notify_one();
        watcher.join();
        return result;
    } else {
        return to_result(lglsat(s));
    }
}

auto LingelingSolver::get_model() const -> ModelMap {
    LGL* s = static_cast<LGL*>(solver);
    std::vector<bool> model_data(std::size_t(m_num_vars), false);
    for(Lit l = 1; l <= m_num_vars; ++l) {
        if(lglderef(s, l) > 0) {
            model_data[l-1] = true;
        }
    }
    return ModelMap{std::move(model_data)};
}

}
