#ifndef SAMMY_VERIFY_H_INCLUDED_
#define SAMMY_VERIFY_H_INCLUDED_

#include "clause_db.h"
#include "greedysat.h"
#include "logging.h"
#include "range.h"
#include "shared_db_propagator.h"
#include "vertex_operations.h"
#include <optional>
#include <sstream>
#include <string>

namespace sammy {

template <typename AssignmentType>
inline std::optional<std::string>
solution_has_error(const ClauseDB& clauses, const AssignmentType& solution,
                   Reason* reason = nullptr) {
    if (solution.size() != clauses.num_vars()) {
        return "Error: Solution length (= " + std::to_string(solution.size()) +
               ") does not match number of variables (= " +
               std::to_string(clauses.num_vars()) + ")!";
    }
    auto literal_satisfied = [&](Lit l) -> bool {
        Var v = lit::var(l);
        bool val = solution[v];
        return val ^ lit::negative(l);
    };
    auto clause_safisfied = [&](const Lit* beg, const Lit* end) {
        return std::any_of(beg, end, literal_satisfied);
    };
    auto vio_unary = std::find_if(clauses.unary_literals().begin(),
                                  clauses.unary_literals().end(),
                                  [&](Lit l) { return !literal_satisfied(l); });
    if (vio_unary != clauses.unary_literals().end()) {
        if (reason) {
            *reason = Reason::Unary{*vio_unary};
        }
        return "Error: Unary clause (" +
               std::to_string(lit::externalize(*vio_unary)) + ") is violated!";
    }
    auto vio_binary =
        std::find_if(clauses.binary_clauses().begin(),
                     clauses.binary_clauses().end(), [&](auto cl) {
                         return !literal_satisfied(cl.first) &&
                                !literal_satisfied(cl.second);
                     });
    if (vio_binary != clauses.binary_clauses().end()) {
        if (reason) {
            *reason = Reason::Binary{vio_binary->first, vio_binary->second};
        }
        return "Error: Binary clause (" +
               std::to_string(lit::externalize(vio_binary->first)) + " " +
               std::to_string(lit::externalize(vio_binary->second)) +
               ") is violated!";
    }
    for (CRef c = 1, n = clauses.literal_db_size(); c < n;
         c = clauses.next_clause(c))
    {
        auto lits = clauses.lits_of(c);
        if (!clause_safisfied(lits.begin(), lits.end())) {
            if (reason) {
                *reason = Reason::Clause{clauses.clause_length(c), c};
            }
            std::ostringstream buffer;
            buffer << "Error: Clause (";
            print_clause_external(buffer, lits);
            buffer << ") is violated!";
            return buffer.str();
        }
    }
    return std::nullopt;
}

inline std::string vertex_to_string(Vertex v) {
    return "(" + std::to_string(lit::externalize(v.first)) + ", " +
           std::to_string(lit::externalize(v.second)) + ")";
}

inline std::optional<std::string>
mutually_exclusive_set_has_error(ClauseDB& clauses, Var n_concrete,
                                 const std::vector<Vertex>& vertices_) {
    // normalize vertices and look for invalid ones
    auto normalize = [](Vertex v) -> Vertex {
        return {(std::min)(v.first, v.second), (std::max)(v.first, v.second)};
    };
    auto is_bad_vertex = [](Vertex v) {
        return lit::var(v.first) == lit::var(v.second);
    };
    auto is_oor_vertex = [&](Vertex v) {
        return lit::var(v.second) >= n_concrete;
    };
    std::vector<Vertex> vertices;
    vertices.reserve(vertices_.size());
    std::transform(vertices_.begin(), vertices_.end(),
                   std::back_inserter(vertices), normalize);
    auto bad_vertex =
        std::find_if(vertices.begin(), vertices.end(), is_bad_vertex);
    if (bad_vertex != vertices.end()) {
        return "Error: Vertex " + vertex_to_string(*bad_vertex) +
               " in mutually exclusive set!";
    }
    auto oor_vertex =
        std::find_if(vertices.begin(), vertices.end(), is_oor_vertex);
    if (oor_vertex != vertices.end()) {
        return "Error: Vertex " + vertex_to_string(*oor_vertex) +
               " (n_concrete = " + std::to_string(n_concrete) +
               ") in mutually exclusive set!";
    }

    // check for duplicates
    EdgeSet present_vertices;
    for (Vertex v : vertices) {
        if (!present_vertices.emplace(v).second) {
            return "Error: Duplicate vertex " + vertex_to_string(v) +
                   " in mutually exclusive set!";
        }
    }

    // verify mutual exclusiveness
    SharedDBPropagator propagator{&clauses};
    std::vector<std::pair<Vertex, Vertex>> explicitly_check_vertices;
    for (std::size_t ind_in_vertices = 0, n_vertices = vertices.size();
         ind_in_vertices < n_vertices; ++ind_in_vertices)
    {
        Vertex v = vertices[ind_in_vertices];
        reset_and_push_noresolve(propagator, v);
        for (std::size_t j = ind_in_vertices + 1; j < n_vertices; ++j) {
            Vertex w = vertices[j];
            if (propagator.is_false(w.first) || propagator.is_false(w.second)) {
                continue;
            }
            bool pushed_w1 = false;
            if (propagator.is_open(w.first)) {
                if (!propagator.push_level(w.first) ||
                    propagator.is_false(w.second))
                {
                    propagator.pop_level();
                    continue;
                }
                pushed_w1 = true;
            }
            // propagator has w.first on trail and w.second isn't false
            if (propagator.is_open(w.second)) {
                bool pres = propagator.push_level(w.second);
                propagator.pop_level();
                if (pushed_w1)
                    propagator.pop_level();
                if (!pres) {
                    continue;
                }
            } else if (pushed_w1) {
                propagator.pop_level();
            }
            // now we're back at the level we had after pushing v
            explicitly_check_vertices.emplace_back(v, w);
        }
    }
    for (auto [v, w] : explicitly_check_vertices) {
        Lit literals[4] = {v.first, v.second, w.first, w.second};
        ClauseDB new_clauses{clauses};
        new_clauses.add_clause(&literals[0], &literals[1]);
        new_clauses.add_clause(&literals[1], &literals[2]);
        new_clauses.add_clause(&literals[2], &literals[3]);
        new_clauses.add_clause(&literals[3], literals + 4);
        SharedDBPropagator pp{&new_clauses};
        GreedySAT<PreferFalse> gsolver{&pp, n_concrete, PreferFalse{}};
        try {
            if (gsolver.solve()) {
                return "Error: Vertices " + vertex_to_string(v) + " and " +
                       vertex_to_string(w) + " are not mutually exclusive!\n";
            }
        } catch (UNSATError&) {
            continue;
        }
    }

    return std::nullopt;
}

inline bool verify_solution(const ClauseDB& clauses,
                            const std::vector<bool>& solution) {
    return !solution_has_error(clauses, solution);
}

} // namespace sammy

#endif
