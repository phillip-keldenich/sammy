#ifndef SAMMY_BARRAGE_LNS_SUBPROBLEM_H_INCLUDED_
#define SAMMY_BARRAGE_LNS_SUBPROBLEM_H_INCLUDED_

#include "literals.h"
#include "dynamic_bitset.h"
#include <memory>

namespace sammy {

/**
 * Structure for general representation of a subproblem
 * as encountered by our LNS.
 */
struct LNSSubproblem {
    std::vector<Vertex> uncovered_universe;
    std::vector<Vertex> mutually_exclusive_set;
    std::vector<DynamicBitset> removed_configurations;
    Lit num_concrete;
};

}

#endif
