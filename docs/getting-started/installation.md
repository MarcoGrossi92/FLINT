# Installation

This guide covers the installation of FLINT.

## Prerequisites

Before installing FLINT, ensure you have the required dependencies:

- **CMake** ≥ 3.23
- **Fortran compiler**: Intel (`ifort`/`ifx`) or GNU (`gfortran`)
- **C++ compiler**: Required for Cantera and TecIO (optional)
- **C compiler**: Required for SUNDIALS (optional)
- **Python** ≥ 3.6: For Python post-process (optional)

FLINT uses two public repositories, both hosted under [Marco Grossi](https://github.com/MarcoGrossi92), which CMake
fetches at configure time (see [ORION and OSLO](#orion-and-oslo)):

- [OSLO](https://github.com/MarcoGrossi92/OSLO): Provides ODE solvers. It builds SUNDIALS when `USE_SUNDIALS` is on
- [ORION](https://github.com/MarcoGrossi92/ORION): Provides the input routines to load the lookup tables in native format

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

    - Static library: `libFLINT.a`, in the build folder
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
| `--use-sundials` | — | Enable SUNDIALS solvers (via OSLO) |
| `--use-tecio` | — | Enable TecIO for binary Tecplot formats (via ORION) |
| `--include-orion` | path | Build a local ORION checkout instead of a release |
| `--include-oslo` | path | Build a local OSLO checkout instead of a release |
| `--update-deps` | — | Select again the newest releases of ORION and OSLO |

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

## ORION and OSLO

`CMakeLists.txt` requires a range of releases of each, not a commit: OSLO `1.1` (1.1.0 or a later 1.x) and ORION `2.1`
(2.1.0 or a later 2.x). The configure takes the newest release (a tag `vX.Y.Z`; pre-releases are skipped) that
satisfies the requirement and clones it in the build folder, under `_deps`. A minor or patch release of ORION or OSLO
reaches FLINT with no change to FLINT.

- **Lock file.** The first configure writes the releases in `deps.lock`; the next configures, of any build folder of
  this checkout, use them without the network. FLINT does not commit it: in a larger build tree the top-level project
  chooses the releases, and FLINT only checks that they satisfy its requirements.
- **Newer releases.** `cmake -DSEMVER_UPDATE=ALL <build>` (or `=ORION`, `=OSLO`), or `./install.sh build
  --update-deps`, selects the newest releases again and rewrites `deps.lock`.
- **Local checkout.** `-DFETCHCONTENT_SOURCE_DIR_ORION=/path/to/ORION` (`--include-orion`) builds that folder as it
  is, and the same for OSLO. The configure warns when its version, its newest tag, does not satisfy the requirement.
- **TecIO and SUNDIALS.** ORION and OSLO build them in their copies, so a new build folder builds them again.
  `-DORION_TECIO_DIR=<folder>` and `-DOSLO_SUNDIALS_DIR=<folder>` keep them outside, for every
  build folder; OSLO rebuilds SUNDIALS when its stamp (the SUNDIALS commit, the compilers and the flags) does not match.
  `install.sh` keeps them in `build-deps/`.

The rules are those of `cmake/SemVerDeps.cmake`, described at its top.


## Using FLINT in a CMake Project

FLINT is built inside the tree of the code that uses it, fetched by CMake, not installed apart. The project requires a
range of releases of FLINT with a copy of `cmake/SemVerDeps.cmake`, included by its path:

```cmake
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SemVerDeps.cmake)
semver_dependency(FLINT GIT_REPOSITORY https://github.com/MarcoGrossi92/FLINT.git VERSION 2.0
                  GIT_SUBMODULES "")      # lib/cantera is cloned by FLINT, if USE_CANTERA is on
FetchContent_MakeAvailable(FLINT)
```

One copy of FLINT, ORION and OSLO serves the whole tree. The first `semver_dependency` of a library selects its
release, and the others only check it: the top-level project requires all its libraries before it makes any of them
available, so that it is the one that selects. When a project of the tree needs a release outside the requirement of
another (ORION 1 for a solver, ORION 2 for FLINT), the configure stops and names both.

A project uses FLINT in one of two ways, and every consumer in a build tree must use the same one: two ways in one
tree compile FLINT's sources twice, into the same module files.

=== "Static library"

    The target `FLINT` (alias `FLINT::FLINT`) is visible from every directory of the tree and is compiled once
    however many targets link it. It is built only for the targets that link it. It carries its dependencies (ORION,
    OSLO, BLAS/LAPACK, and Cantera and OpenMP when enabled) and its module folder:

    ```cmake
    target_link_libraries(your_target PUBLIC FLINT)
    ```

    Link it `PUBLIC` when the modules of `your_target` use FLINT's: ifx reads FLINT's module files when it compiles a
    user of those modules. This is the way for a tree where several projects use FLINT.

=== "Source files"

    The list `FLINT_SRC` holds FLINT's source files. The consumer compiles them in one of its own targets, and
    `link_FLINT_dependencies` links what they need:

    ```cmake
    add_library(your_lib ${your_sources} ${FLINT_SRC})
    link_FLINT_dependencies(your_lib)
    ```

    `FLINT_SRC` is set only in the directory that adds FLINT, so only one project of a tree can use it. FLINT's
    directory-level definitions (`CANTERA`) and its OpenMP flag do not reach the consumer's target: the consumer has to
    give them itself.

The top-level project sets the flags before it adds FLINT: `CMAKE_Fortran_FLAGS_<CONFIG>`,
`CMAKE_EXE_LINKER_FLAGS_<CONFIG>` and the options `USE_OPENMP`, `USE_MPI` and `USE_CANTERA` (OSLO's own defaults
turn OpenMP and MPI on). FLINT's objects are compiled with them, the same as the top-level project's own.

!!! warning "Intel IPO"
    With interprocedural optimization (`-ipo`), FLINT's objects hold IR, not machine code. The link of each executable
    then optimizes them together with the code that calls them, and FLINT's small thermodynamic functions are inlined
    into the solver. This needs `-ipo` in the compile and link flags, and an archiver that indexes IR objects
    (`llvm-ar` and `llvm-ranlib` as `CMAKE_AR` and `CMAKE_RANLIB`), all set before FLINT is added. A FLINT compiled
    apart, or with other flags, is left out of that optimization: the solver runs slower, and nothing warns.

FLINT does not install a CMake package: `find_package(FLINT)` is not supported.

## Next Steps

- **[Quick Start Tutorial](quick-start.md)**: Learn FLINT basics with a hands-on example
- **[User Guide](../user/index.md)**: Explore the converter and API capabilities

---