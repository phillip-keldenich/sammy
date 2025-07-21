# Sammy
This repository contains the source code for Sammy,
our software for computing minimum samples with pairwise interaction coverage.
It also contains the data from our experiments as well as scripts to re-run
the experiments (though the experiment scripts are somewhat specialized to our
cluster's scheduler, they should be runnable but will consume very significant
time and memory).
Some python packages have to be installed to run the experiments,
but they are all available on PyPI.

# Setup Instructions
To build Sammy, which is a C++ program, 
we use `CMake` and `conan` (version 2) as a package manager;
you will have to install both of these first (as well as a C++ compiler),
as the following instructions assume they are available in your `PATH`.
`conan` can be installed into a python virtual environment
using `pip` (e.g., `pip install conan`), but it will put C++ packages into
a user-specific directory such as `~/.conan2`, called the conan cache.

## Installing Dependencies
Most of our C++ dependencies are directly available in `conan` (through the
default `conan-center` repository), but `gurobi` is not (the repository
includes our own `gurobi_public` conan package, which must be manually installed
before trying to build Sammy).

### Installing `gurobi_public`
The first step should thus be to install the `gurobi_public` package.
Assuming you start in the directory where this `README.md` is located,
you can run the following command to install it:
```bash
cd custom_deps/gurobi_public
conan create . -s build_type=Release --user sammy --channel stable --build=missing
cd ../..
```

### Build Preparations
After installing `gurobi_public`, you can prepare the build environment
by running the following command in the directory where this 
`README.md` is located.
```bash
conan install . -s build_type=Release --build=missing
```
This will download and build all remaining required dependencies.
It will also create a build directory for `CMake` at `build/Release`.

### CMake Configure Step
You are now ready to run the `CMake` configuration step.
This can be done by running the following command in the directory where this
`README.md` is located:
```bash
cmake -B build/Release -S . -DCMAKE_TOOLCHAIN_FILE=build/Release/generators/conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release
```
This should find all the dependencies and prepare the build system.

### Build Step
Finally, you can build Sammy by running the following command in the directory
where this `README.md` is located.
```bash
cmake --build build/Release --parallel 8
```
This will build Sammy and its required executables, as well as some tests 
and other executables that are used for further experiments.

## Gurobi License
To actually run Sammy (and not just the initial heuristic), you will need
a valid Gurobi license since we use Gurobi as an LP solver in the later stages
of our lower bound approach.
You can obtain a Gurobi license (free for academic users) from the official
[Gurobi website](https://www.gurobi.com/).
To install the Gurobi license, you will usually be instructed to run the command
```bash
grbgetkey <your-license-key>
```
The `grbgetkey` command will be installed to your conan cache 
as part of the `gurobi_public` package, but it is not added to your `PATH`.
To add it to your `PATH` temporarily, you can use the following command
(assuming a unix system; note the double periods and space):
```bash
. ./build/Release/generators/conanrun.sh
```
You should afterwards be able to run `grbgetkey` from the command line.

If you installed the Gurobi license in a non-standard location,
you will need to set the `GUROBI_LICENSE_FILE` environment variable
to point to the license file.

# Running Sammy
After building Sammy, you can solve instances in the `.scm.json` format
(optionally with `.xz`, `.bz2` or `.gz` compression) 
using the `sammy_solve` executable.
For example, in the directory where this `README.md` is located,
try the following command:
```bash
./build/Release/src/sammy_solve sammy_benchmark_instances/soletta_2017-03-09_21-02-40.scm.json.xz --print-events
```
Without the `-o filename.json` option, this will not create an output file,
but it will print event information to the console.

# Instances
The benchmark instances used in our experiments are in the 
`./sammy_benchmark_instances/` directory.
The full instance set is in `./full_instances/`.
These instances are in the `.scm.json` format, which is a JSON-based format
that is used by Sammy to represent instances.
Their counterparts in the original `.xml` or `.dimacs` formats are also
contained in the same directory, with the same name but a different extension.

# Instance Formats
To transform instances from other formats (standard `.xml` or `.dimacs` formats)
into the format accepted by Sammy, you can use the python script
`./scripts/instance_to_generic_json.py`; this script uses the parser from
`samplns` to generate instances in the `.scm.json` format.

# Solution Formats
To transform solutions from the JSON format produced by Sammy,
you can use the script `./scripts/solutions_to_feature_name_mapping.py`.
It will produce a list of mappings, one for each configuration in the sample,
that assigns a Boolean value to each concrete feature.

# Experiment Scripts
There are several experiments in the `./experiments/` directory.
Each experiment has its own subdirectory; each subdirectory contains
a list of files with some name prefixed with some number.
The python scripts should be executed in the order indicated by the numbers
to rerun our experiments.
For instance, `experiments/02_run_default_params` contains the script
`01_run.py`, which will produce outputs in the `02_output` directory,
which can then be checked by running the `03_check_solutions.py` script.

To reproduce the figures and tables in our paper on other computers,
the subdirectory `paper_figures` can be used.
It contains scripts/instructions to extract the data from the experiment
outputs and a notebook to visualize the data and save plots as PDF files.

The subdirectories `06_identify_nontrivial` and 
`07_run_nontrivial_export_subproblems` contain scripts and notebooks
to produce the data for the figures in the appendix.

## Extracting our raw outputs
To emplace the raw output files from our experiments as we ran them,
you can run the script `./experiments/extract_author_experiment_outputs.py`.
This essentially replaces the running of the experiments
`02_run_default_params`, `03_baselines`, `05_run_initial_no_simplification`,
`06_identify_nontrivial`, and `07_run_nontrivial_export_subproblems`,
and places the raw output data as well as summary files into those directories.
Afterwards, the `paper_figures` scripts can be run to produce the figures and
tables in the paper.
Note that this can take quite a while and produces a large number of files.