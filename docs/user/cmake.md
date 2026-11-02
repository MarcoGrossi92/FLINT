# Using FLINT in a CMake Project

FLINT is built inside the tree of the code that uses it, fetched by CMake, not installed apart
([Quick Start](../getting-started/quick-start.md) builds a first program in this way). The project requires a range of
releases of FLINT with a copy of `cmake/SemVerDeps.cmake`, included by its path:

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
