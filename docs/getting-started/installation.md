# Installation

This guide covers the installation of FLINT.

## Prerequisites

Before installing FLINT, ensure you have the required dependencies:

- **CMake** ≥ 3.23
- **Fortran compiler**: Intel (`ifort`/`ifx`) or GNU (`gfortran`)
- **C++ compiler**: Required for Cantera and TecIO (optional)
- **C compiler**: Required for SUNDIALS (optional)
- **Python** ≥ 3.6: For Python post-process (optional)

The following public repositories are included as Git submodules and are automatically configured during installation:

- [OSlo](https://github.com/MarcoGrossi92/OSlo): Provides ODE solvers. SUNDIALS is included
- [ORION](https://github.com/MarcoGrossi92/ORION): Provides the input routines to load the lookup tables in native format

(both hosted under [Marco Grossi](https://github.com/MarcoGrossi92))

## Installation Methods

Choose the installation method that best fits your needs:

=== "Installer"

    Build the Fortran library with Cantera support:

    ```bash
    git clone https://github.com/MarcoGrossi92/FLINT.git
    cd FLINT
    ./install.sh build --compiler=gnu --use-cantera --use-sundials
    ```

=== "Manual CMake"

    Build the Fortran library and command-line converter:

    ```bash
    git clone https://github.com/MarcoGrossi92/FLINT.git
    cd FLINT
    mkdir build
    cd build

    cmake .. \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_Fortran_COMPILER=gfortran \
        -DUSE_CANTERA=ON \
        -DUSE_SUNDIALS=ON

    cmake --build . --parallel
    ```

    This creates:

    - Static library: `lib/libFLINT.a`
    - Test executables: `bin/test/`

## Installer Script

The `install.sh` script provides flexible build configuration:

**Basic Syntax**

```bash
./install.sh [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]
```

**Global Options**

| Option | Description |
|--------|-------------|
| `-v, --verbose` | Enable verbose output for debugging |

## Commands

**build**

Perform a complete build from scratch:

```bash
./install.sh build [OPTIONS]
```

**Options:**

| Option | Values | Description |
|--------|--------|-------------|
| `--compiler` | `gnu`, `intel` | Select compiler suite (default: `gnu`) |
| `--use-cantera` | — | Enable Cantera |
| `--use-sundials` | — | Enable SUNDIALS solvers (via OSlo) |
| `--use-tecio` | — | Enable TecIO for binary Tecplot formats (via ORION) |

!!! warning "SUNDIALS & Cantera"
    To properly compile Cantera SUNDIALS is requested.

**Examples:**

```bash
# Build with GNU compiler
./install.sh build --compiler=gnu

# Build with Intel compiler and all options
./install.sh build --compiler=intel --use-sundials --use-cantera --use-tecio
```

**compile**

Recompile using existing CMake configuration stored in `CMakePresets.json`:

```bash
./install.sh compile
```

Use this after modifying source code to avoid reconfiguring the build system.

## Cantera Support

Cantera is an open-source suite of tools for problems involving chemical kinetics, thermodynamics, and transport processes. It is shipped with a Fortran 90 interface that can be linked to FLINT to fully employ Cantera capabilities.

Cantera will not rebuilt unless its folder `lib/cantera` is cleaned. 

!!! note "C++ Compiler Required"
    Cantera requires a C++ compiler. Ensure `g++` or Intel C++ compiler is available.

!!! warning "Build Cantera"
    Cantera compilation may be a little tricky. If you experiment some problems, try to build Cantera standalone before compiling it within FLINT framework.


## Library Linking (Advanced)

To link FLINT in external Fortran projects:

```bash
# Compile your program
gfortran -I/path/to/FLINT/include \
         -L/path/to/FLINT/lib \
         -lFLINT \
         your_program.f90 -o your_program
```

Or use CMake's `find_package`:

```cmake
find_package(FLINT REQUIRED)
target_link_libraries(your_target FLINT::FLINT)
```

## Next Steps

- **[Quick Start Tutorial](quick-start.md)**: Learn FLINT basics with a hands-on example
- **[User Guide](../user/index.md)**: Explore the converter and API capabilities

---