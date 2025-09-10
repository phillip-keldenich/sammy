FROM ubuntu:24.04

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3-pip \
    cmake \
    git \
    build-essential \
    python3 \
    python3-venv \
    libboost-all-dev \
    && rm -rf /var/lib/apt/lists/*

# Set workdir
WORKDIR /app

# Copy source code
COPY . /app

# Ensure submodules are initialized
RUN git submodule update --init --recursive

# Create a Python virtual environment for Conan
RUN python3 -m venv /app/conan-venv
ENV PATH="/app/conan-venv/bin:$PATH"

# Install Conan in the virtual environment
RUN /app/conan-venv/bin/pip install --upgrade pip
RUN /app/conan-venv/bin/pip install conan

# Create default Conan profile
RUN conan profile detect --force

# Build Gurobi via local conan recipe
RUN conan create ./custom_deps/gurobi_public -pr:b default -pr:h default -s build_type=Release --build=missing --user sammy --channel stable

# Install dependencies via conanfile.docker.txt
RUN conan install ./conanfile.docker.txt -pr:b default -pr:h default --build=missing

# Configure the project with CMake
RUN cmake --preset conan-release

# Build the project (only 2 threads to keep memory usage low)
RUN cmake --build --preset conan-release --parallel 2

# Run tests (optional, can be commented out)
# RUN ./build/Release/test/sammy_test --success

# Default command: run the main program
ENTRYPOINT ["./build/Release/src/sammy_solve"]