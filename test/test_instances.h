#ifndef SAMMY_TEST_INSTANCES_H_INCLUDED_
#define SAMMY_TEST_INSTANCES_H_INCLUDED_

#include <vector>
#include <sammy/literals.h>

namespace sammy {

struct VarCount {
    Var n_all;
    Var n_concrete;
};

extern std::vector<ExternalClause> soletta_2017_03_01_15_25_44_clauses();
extern VarCount soletta_2017_03_01_15_25_44_variables();

}

#endif
