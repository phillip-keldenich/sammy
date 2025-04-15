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