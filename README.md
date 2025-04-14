# Write a better readme

* [run_initial.h](./include/sammy/run_initial.h): Contains the code to create the initial solution and determine the coverable interactions.
* [fast_clique.h](./include/sammy/fast_clique.h): Contains the code to quickly find a maximal clique in the graph. This is primarily used to find mutually excluding interactions that can be used as lower bound or for symmetry breaking.
* [literals.h](./include/sammy/literals.h): Basic types for literals, variables, interactions, and more.
* [pair_infeasibility_map.h](./include/sammy/pair_infeasibility_map.h): A container for storing which pairs of interactions are infeasible and which are feasible (and which are unknown).

## Key Files corresponding to Document

### 3.3 Tracking Coverage

* [pair_infeasibility_map.h](./include/sammy/pair_infeasibility_map.h): A container for storing which pairs of interactions are infeasible and which are feasible (and which are unknown).
    * [dynamic_bitset.h](./include/sammy/dynamic_bitset.h): A dynamic bitset implementation as the one from the standard library was not optimal for our needs.


## Probably unused files

* [primal_dual_driver.h](./include/sammy/primal_dual_driver.h)