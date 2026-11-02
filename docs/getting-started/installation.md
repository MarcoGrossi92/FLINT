# Installation

This guide covers the installation of FLINT.

## Prerequisites

Before installing FLINT, ensure you have the required dependencies:

- **CMake** ≥ 3.23
- **Git**, and network access at the first configure
- **Fortran compiler**: Intel (`ifort`/`ifx`) or GNU (`gfortran`)
- **C++ compiler**: Required for Cantera and TecIO (optional)
- **C compiler**: Required for SUNDIALS (optional)
- **Python** ≥ 3.10 and **SCons** ≥ 4.0: Required to build Cantera (optional). Cantera also needs the headers of
  Boost 1.70 or newer: those installed on the machine, or else a copy that the configure downloads

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
    ./install.sh build --compilers=gnu --use-cantera --use-sundials
    ```

=== "Manual CMake"

    Build the Fortran library and the test programs, with Cantera support:

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
| `--compilers` | `gnu`, `intel` | Compiler suite: gfortran, gcc, g++ or ifx, icx, icpx (default: those of `FC`, `CC`, `CXX`, or those CMake finds) |
| `--use-cantera` | — | Enable Cantera |
| `--use-sundials` | — | Enable SUNDIALS solvers (via OSLO) |
| `--use-tecio` | — | Enable TecIO for binary Tecplot formats (via ORION) |
| `--include-orion` | path | Build a local ORION checkout instead of a release |
| `--include-oslo` | path | Build a local OSLO checkout instead of a release |
| `--update-deps` | — | Select again the newest releases of ORION and OSLO |

`--use-cantera` turns `--use-sundials` on: Cantera is built on the SUNDIALS of OSLO.

**Examples:**

```bash
# Build with GNU compiler
./install.sh build --compilers=gnu

# Build with Intel compiler and all options
./install.sh build --compilers=intel --use-sundials --use-cantera --use-tecio
```

`install.sh build` writes `CMakePresets.json`; after a change to the sources, `cmake --build build` recompiles
without configuring again.

## Cantera Support

Cantera is an open-source suite of tools for problems involving chemical kinetics, thermodynamics, and transport processes. It ships a Fortran 90 interface, which FLINT links to compare its own results with Cantera's.

`--use-cantera` (`-DUSE_CANTERA=ON`) builds Cantera 3.2.0 from the submodule `lib/cantera` with SCons, at the first
build. It needs a C++ compiler, Python 3.10 or newer and SCons 4.0 or newer; the Boost headers are downloaded when the
machine has none. The next builds reuse it: the configure builds Cantera again only when its libraries or its Fortran
module are missing from `lib/cantera/build` (`scons clean` in `lib/cantera` removes them).

!!! tip "When the Cantera build fails"
    The build log shows the SCons step that failed, and `lib/cantera/config.log` the details of Cantera's configure
    checks.

## ORION and OSLO

`CMakeLists.txt` requires a range of releases of each, not a commit. The configure takes the newest release (a tag `vX.Y.Z`; pre-releases are skipped) that
satisfies the requirement and clones it in `lib/OSLO` and `lib/ORION`, out of the build folder: a new build folder
uses them again, without the network. A minor or patch release of ORION or OSLO reaches FLINT with no change to FLINT.

- **Sources in `lib/`.** A checkout at the selected release is built as it is, local changes included. When another
  release is selected, the configure moves the checkout to it, but only if it is as the configure left it (a detached
  HEAD without changes); on a branch or with changes it stops, and the error tells how to go on.
- **Offline.** As for FetchContent, `-DFETCHCONTENT_FULLY_DISCONNECTED=ON` forbids any clone or fetch, and
  `-DFETCHCONTENT_UPDATES_DISCONNECTED=ON` any in an existing checkout: the configure stops when one is needed.

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


## Verify the Build

After a build, the test programs are in `bin/test` and registered with CTest. From the root of the repository, run
the quick tier of the suite (about 30 s):

```bash
ctest --test-dir build -L quick
```

or every test with `ctest --test-dir build -j8`. Each test prints its checks and a verdict, and exits with code 1 on
failure; CTest summarizes them. A single test program can also be run by hand from its folder in `test/`, for
example the real-fluid (p, h) tables and their inverse (p, T) -> h:

```bash
cd test/thermo
../../bin/test/test-thermo-real
```

Without Cantera, the tests compare FLINT with Cantera references stored in the repository. See
[Testing](../development/testing.md) for the areas, the labels and the reference data.

## Using FLINT in Your Code

A project that uses FLINT fetches it with CMake and builds it in its own tree: see
[Using FLINT in a CMake Project](../user/cmake.md), and [Quick Start](quick-start.md) for a first program.

## Next Steps

- **[Quick Start Tutorial](quick-start.md)**: Learn FLINT basics with a hands-on example
- **[User Guide](../user/index.md)**: Embed FLINT in your code

---