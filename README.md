# SAMMY - High-Quality Interaction Sampling with Bounds

## Build and Run

```bash
# Make sure you got all submodules
git submodule update --init --recursive
# Install conan
pip install -U conan
# Install Gurobi via local conan recipe
conan create ./cmake/conan/gurobi_public -pr:b default -pr:h default  -s build_type=Release 
--build=missing
# Install dependencies via conanfile.txt
conan install . -pr:b default -pr:h default --build=missing
# Configure the project with CMake
cmake --preset conan-release
# Build the project
cmake --build --preset conan-release
# Run the program
./build/Release/src/sammy_solve
```

If you get issues with building the test cases, you can run the following command to build only the main program:

```bash
cmake --preset conan-release -DSAMMY_SKIP_TESTS=ON && cmake --build --preset conan-release     
```

# Write a better readme

* [run_initial.h](./include/sammy/run_initial.h): Contains the code to create the initial solution and determine the coverable interactions.
* [fast_clique.h](./include/sammy/fast_clique.h): Contains the code to quickly find a maximal clique in the graph. This is primarily used to find mutually excluding interactions that can be used as lower bound or for symmetry breaking.
* [literals.h](./include/sammy/literals.h): Basic types for literals, variables, interactions, and more.
* [pair_infeasibility_map.h](./include/sammy/pair_infeasibility_map.h): A container for storing which pairs of interactions are infeasible and which are feasible (and which are unknown).

## Key Files corresponding to Document

### 3.3 Tracking Coverage

* [pair_infeasibility_map.h](./include/sammy/pair_infeasibility_map.h): A container for storing which pairs of interactions are infeasible and which are feasible (and which are unknown).
    * [dynamic_bitset.h](./include/sammy/dynamic_bitset.h): A dynamic bitset implementation as the one from the standard library was not optimal for our needs.


### Simplification

Direct dependencies from the function implementation:
* [solve.cpp](./src/sammy/solve.cpp): Contains the `run_simplification` function definition.
* [simplification.h](./include/sammy/simplification.h): Contains the `run_simplifier` function that is called by `run_simplification`.

#### Dependencies from `simplification.h`:

* [bounded_variable_elimination.h](./include/sammy/bounded_variable_elimination.h): For bounded variable elimination functionality.
* [compress_binaries.h](./include/sammy/compress_binaries.h): For binary clause compression.
* [detect_equalities.h](./include/sammy/detect_equalities.h): For detecting literal equalities.
* [eliminate_subsumed.h](./include/sammy/eliminate_subsumed.h): For subsumption elimination.
* [simplify_datastructure.h](./include/sammy/simplify_datastructure.h): For the core simplification data structures.
* [vivify.h](./include/sammy/vivify.h): For clause vivification.

## Probably unused files

* [primal_dual_driver.h](./include/sammy/primal_dual_driver.h)
* [instance.h](./include/sammy/instance.h)
* [ortools_solver.h](./include/sammy/ortools_solver.h)
* [problem_input.h](./include/sammy/problem_input.h)
* [sammy.h](./include/sammy/sammy.h)
# Sammy

This repository contains code and data for the paper "TBD".

## Installation

In order to install the package, install all the dependencies with conan:

```bash
conan install . --build=missing
```

This will generate a toolchain file in the `build/` directory. 
You can then use this file to build the project with CMake.
Then, you can build the project with:

```bash
cmake -B build/Release -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build/Release --config Release
```

## Usage
TBD

## Project Structure

The project is organized as follows:

* `src/`: Contains the source code for the project.
* `include/`: Contains the header files for the project.
* `test/`: Contains the test files for the project.
* `external/`: Contains the external libraries used in the project.
* `experiments/`: Contains the scripts and data for the experiments.
