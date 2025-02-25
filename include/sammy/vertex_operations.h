#ifndef SAMMY_VERTEX_OPERATIONS_H_INCLUDED_
#define SAMMY_VERTEX_OPERATIONS_H_INCLUDED_

#include "literals.h"
#include "shared_db_propagator.h"

namespace sammy {

/**
 * @brief Reset a propagator and push a vertex to it.
 *        Does not use conflict resolution and
 *        throws an exception if the vertex cannot
 *        be pushed without conflict.
 *
 * @param propagator
 * @param vertex
 */
inline void reset_and_push_noresolve(SharedDBPropagator& propagator,
                                     Vertex vertex) 
{
    propagator.reset_or_throw();
    Lit lmin = vertex.first;
    Lit lmax = vertex.second;
    if (propagator.is_false(lmin) || propagator.is_false(lmax)) {
        throw std::logic_error("Infeasible vertex pushed!");
    }
    if (propagator.is_open(lmin)) {
        if (!propagator.push_level(lmin) || propagator.is_false(lmax)) {
            throw std::logic_error("Infeasible vertex pushed!");
        }
    }
    if (propagator.is_open(lmax)) {
        if (!propagator.push_level(lmax)) {
            throw std::logic_error("Infeasible vertex pushed!");
        }
    }
}

/**
 * @brief Tests whether the given vertex can be pushed onto the given propagator.
 * Does not learn from conflicts; returns the propagator into the previous state.
 */
inline bool can_push(SharedDBPropagator& propagator, Vertex vertex) {
    Lit lmin = vertex.first;
    Lit lmax = vertex.second;
    if (propagator.is_false(lmin) || propagator.is_false(lmax)) {
        return false;
    }
    int push_count = 0;
    if (propagator.is_open(lmin)) {
        if (!propagator.push_level(lmin)) {
            propagator.pop_level();
            return false;
        }
        if (propagator.is_false(lmax)) {
            propagator.pop_level();
            return false;
        }
        ++push_count;
    }
    bool result = true;
    if (propagator.is_open(lmax)) {
        result = propagator.push_level(lmax);
        ++push_count;
    }
    for (int i = 0; i < push_count; ++i) {
        propagator.pop_level();
    }
    return result;
}

inline int push_vertex(SharedDBPropagator& propagator, Vertex vertex) {
    Lit lmin = vertex.first;
    Lit lmax = vertex.second;
    if (propagator.is_false(lmin) || propagator.is_false(lmax)) {
        return -1;
    }
    int pushed = 0;
    if (propagator.is_open(lmin)) {
        if (!propagator.push_level(lmin) || propagator.is_false(lmax)) {
            propagator.pop_level();
            return -1;
        }
        pushed = 1;
    }
    if(propagator.is_open(lmax)) {
        if(!propagator.push_level(lmax)) {
            propagator.pop_level();
            if(pushed) propagator.pop_level();
            return -1;
        }
        ++pushed;
    }
    return pushed;
}

inline bool push_vertex_pair(SharedDBPropagator& propagator, Vertex v1, Vertex v2) {
    int pc1 = push_vertex(propagator, v1);
    if(pc1 < 0) return false;
    if(push_vertex(propagator, v2) < 0) {
        for(int i = 0; i < pc1; ++i) propagator.pop_level();
        return false;
    }
    return true;
}

}

#endif
