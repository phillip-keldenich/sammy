# Sammy
This repository contains the source code for Sammy,
our software for computing minimum samples with pairwise interaction coverage.
It also contains the data from our experiments as well as scripts to re-run
the experiments (though the experiment scripts are somewhat specialized to our
cluster's scheduler, they should be runnable but will consume very significant
time and memory on a single machine).
Some python packages have to be installed to run the experiments,
but they are all available on PyPI.

# Long-term availability
There is a long-term available snapshot of the code,
as well as all raw experiment data stored at 
[Zenodo](https://doi.org/10.5281/zenodo.16277857).
This snapshot contains two archives: one contains a snapshot of the source code
and instances and everything that should be necessary to reproduce our experiments
(`sammy_evaluation_package.zip`), while the other contains the source code
and all the raw data resulting from our experiments (`sammy_package.zip`),
which is significantly larger and lacks files we added later to ease reproduction.

> [!NOTE]
> If you downloaded this as a ZIP archive from Zenodo,
> be aware that there is a GitHub repostory containing
> the source code, which may have been updated since the Zenodo snapshot.
> The GitHub repository is at:
> https://github.com/phillip-keldenich/sammy

# Reproduction of Experiments
We offer several possibilities to reproduce our experiments to varying degrees.
There are essentially three levels of reproducibility that should work on sufficiently
recent Linux or MacOS systems.

1. Either using Docker (follow `Simple Docker Usage` below for setup instructions) 
   or using the `Direct Setup Instructions` below, you can build **Sammy**.
   Depending on which option you choose, you can then either execute `run_level1_reproduction_direct.sh` or `run_level1_reproduction_docker.sh` to run **Sammy** five times on each of the
   benchmark instances from the benchmark set (in the `sammy_benchmark_instances/` directory)
   used in our paper.
   This allows you to reproduce most parts of Table E.1 in the full version of the paper.
   Note that the column counting the number of interactions is expected to differ from the paper,
   since it (in your run) counts the interactions after simplification (but before universe reduction),
   unlike the table in the paper.
   To reproduce the number without simplification, edit the script to pass `--dont-simplify` to the invokations of `sammy_solve`.

   Note that either option requires you to have a valid Gurobi license file.
   In case of the Docker option, this needs to be a Docker-compatible license file,
   which has to be placed as `gurobi.lic` in the same directory as this `README.md`.
   In case of the direct setup option, any valid Gurobi license file will do,
   as long as the `GRB_LICENSE_FILE` environment variable is properly set up.

   Either way, running this script will produce an Excel table with the results of the runs on your machine,
   as well as a subdirectory containing the raw output files from each run.
   Note that this will take a few days on a single machine, but the time should still
   be relatively reasonable.

2. To reproduce more experiments from our paper, including the baselines,
   and reproduce most of the figures and tables in the paper,
   after following the `Direct Setup Instructions` below and running `pip install -r requirements.txt`,
   you can run the `run_level2_reproduction.sh` script.
   Note that by default, this will clone the current version of SampLNS from github.
   Should the repository of SampLNS become unavailable in the future,
   we have included a snapshot of the version of SampLNS that we used in the archive,
   named `SampLNS.tar.xz`.
   This will take a significant amount of time on a single machine,
   and relies on enough memory being available; in our experiment
   environment, the scripts distribute the work across 6 identical machines,
   which still may take 1-2 days to complete the experiments.
   The figures that are produced can be found in the `experiments/paper_figures` directory.

3. To reproduce almost all the experiments from our paper, including the 
   identification of non-trivial instances in the full instance set,
   as well as exporting and solving the individual LNS repair subproblems,
   you can run the `run_level3_reproduction.sh` script.
   This will probably be infeasible on a single machine;
   again, in our experiment environment, the scripts distribute the work
   across 6 identical machines, where this can still take 2 weeks.

# Software Versions Used in Experiments
In our experiments, we used the following software versions:
- Ubuntu Linux 24.04.2 LTS
- Compiler: `clang++ 18.1.3 (1ubuntu1)` with flags `-std=c++17 -mavx2 -mpopcnt -mlzcnt` 
on top of the release mode flags introduced by CMake and our `CMakeLists.txt` file,
- `cmake 3.28.3`
- `conan 2.12.1`
- `gurobi 12.0.2` as installed by the bundled `gurobi_public` conan package,
- `Boost 1.88.0` installed via conan (from conan-center),
- `nlohmann_json 3.12.0` installed via conan (from conan-center),
as well as the versions of the external SAT solvers CaDiCaL, Lingeling, Cryptominisat and Kissat,
which are included in this archive.

# Simple Docker Usage
We provide a simplified procedure for running **Sammy** using Docker.
To proceed, ensure that **Docker** is installed on your system and that you have a Docker-compatible Gurobi license file.  

<details>
<summary>Gurobi Academic License</summary>

1. Visit the [Gurobi Academic Program](https://www.gurobi.com/academia/academic-program-and-licenses/) and click **Claim Your Free License Now**.  
2. Sign up with your academic email address.  
3. In the *Licenses* tab, click **Request**.  
4. Select **WLS Academic**, which works in Docker containers but requires access to an academic network (e.g., university VPN).  
   Gurobi describes the WLS Academic license as:  
   - Runs on multiple machines/containers  
   - Requires Gurobi v10+  
   - Needs an internet connection during usage  
   - Academic network required only for generation; usable anywhere afterwards  
   - Valid for 90 days, renewable while eligible  
5. Click **Go to License**.  
6. On the right, click **Open License within the Web License Manager** (icon).  
7. Click **Download**.  
8. Enter any application name (e.g., `sammy`).  
9. Save the `gurobi.lic` file in the same directory as this `README.md`.  

*Note: The process may change; consult the official Gurobi documentation if you encounter issues.*  
</details>

If your problem instances are in FeatureIDE’s XML format or the DIMACS format, convert them to the `.scm.json` format using the script
`scripts/instance_to_generic_json.py`. Place all `.scm.json` files you want to solve in a directory named `instances/` located in the same directory as this `README.md`.
For testing, we recommend to just copy `full_instances/berkeleyDB1.scm.json.xz`.

Also place your Gurobi license file as `gurobi.lic` in the same directory as this `README.md`.

## Building Docker image
You need to build the Docker image yourself by executing the following command in this directory:  

```bash
docker build --platform linux/amd64 -t sammy .
```

Unfortunately, we cannot provide a pre-built Docker image due to licensing restrictions of Gurobi;
we are not allowed to redistribute the Gurobi shared libraries.
This holds for SampLNS as well, which is why we opted to include a snapshot of its
repository for long-term archival purposes.

## Using the Docker image
After building the image, run:

```bash
./solve_instances.sh
```

This script solves all instances in the `instances/` directory and stores the results in the `results/` directory.

You may edit `solve_instances.sh` to adjust the arguments passed to Sammy or to change the locations of the instances, results, or license file.
By default, the script enforces a time limit of 3600 seconds per instance.

> [!NOTE]
>
> To reproduce the basic experiments from our paper (Table E.1), run:
> ```bash
> ./run_level1_reproduction_docker.sh
> ```  
> This script executes each instance five times and stores the results in the `level1_docker_results/` directory, including an Excel summary file.
> Running all experiments may take several days on a single machine.
> Since some benchmark instances are large, a workstation with sufficient memory is required (tested with 96 GiB RAM).
> The column reporting the number of feasible interactions will deviate significantly from the table in the paper,
> since it counts the number of interactions after simplification (but not universe reduction).

> [!WARNING]
>
> Killing the script can be difficult due to Docker's handling of signals.
> `while docker kill $(docker ps -q --filter ancestor=sammy) 2>/dev/null; do sleep 0.5; done` seems to work in most cases.

# Direct Setup Instructions
To build Sammy (which is a C++ program) without using Docker,
offering more flexible access to individual commands,
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
you will need to set the `GRB_LICENSE_FILE` environment variable
to point to the license file.

## Running Sammy
After building Sammy, you can solve instances in the `.scm.json` format
(optionally with `.xz`, `.bz2` or `.gz` compression) 
using the `sammy_solve` executable.
For example, in the directory where this `README.md` is located,
try the following command:
```bash
./build/Release/src/sammy_solve sammy_benchmark_instances/soletta_2017-03-09_21-02-40.scm.json.xz --print-events
```
Without the `-o filename.json[.xz]` option, this will not create an output file,
but it will print event information to the console.

## Instances
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
