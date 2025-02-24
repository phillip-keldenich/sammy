#include <sammy/kissat_solver.h>

extern "C" {
#include "kissat.h"
}

#include <iostream>

namespace sammy {

KissatSolver::KissatSolver() :
    solver(static_cast<void*>(kissat_init())),
    m_num_vars(0)
{}

KissatSolver::~KissatSolver() {
    kissat_release(static_cast<kissat*>(solver));
}

void KissatSolver::reserve(Lit nv) {
    m_num_vars = (std::max)(nv, m_num_vars);
    kissat_reserve(static_cast<kissat*>(solver), nv);
}

void KissatSolver::add_literal(Lit l) {
    Lit v = (std::abs)(l);
    if(v > m_num_vars) {
        m_num_vars = v;
    }
    kissat_add(static_cast<kissat*>(solver), l);
}

void KissatSolver::finish_clause() {
    kissat_add(static_cast<kissat*>(solver), 0);
}

std::optional<bool> KissatSolver::solve() {
    if(already_solved) 
        throw std::logic_error("Already solved (KISSAT does not support incremental solving)!");
    already_solved = true;
    kissat* s = static_cast<kissat*>(solver);
    int state = kissat_solve(s);
    if(state == 10) {
        // SAT
        /*solution.assign(m_num_vars, false);
        for(Var v = 1; v <= m_num_vars; ++v) {
            Lit l = kissat_value(s, v);
            if(l == v) {
                solution[v-1] = true;
            }
        }*/
        return true;
    } else if(state == 20) {
        // UNSAT
        return false;
    } else {
        // UNKNOWN/LIMIT HIT
        return std::nullopt;
    }
}

auto KissatSolver::get_model() const -> ModelMap {
    if(!already_solved) throw std::logic_error("Called get_model on unsolved formula!");
    kissat* s = static_cast<kissat*>(solver);
    std::vector<bool> solution(m_num_vars, false);
    for(Lit l = 1; l <= m_num_vars; ++l) {
        Lit val = kissat_value(s, l);
        if(val == l) solution[l-1] = true;
    }
    return ModelMap{std::move(solution)};
}

void KissatSolver::set_conflict_limit(unsigned limit) {
    kissat_set_conflict_limit(static_cast<kissat*>(solver), limit);
}

void KissatSolver::set_decision_limit(unsigned limit) {
    kissat_set_decision_limit(static_cast<kissat*>(solver), limit);
}

void KissatSolver::terminate() {
    kissat_terminate(static_cast<kissat*>(solver));
}

}
