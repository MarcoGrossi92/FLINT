# Changelog

All notable changes to FLINT are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html): see
[Semantic Versioning](https://marcogrossi92.github.io/FLINT/development/versioning/) in the documentation for what
a major, a minor and a patch release of FLINT may change.

## [Unreleased]

Changes on `main` since `v1.1.0`. Two of them break the projects that use FLINT (**Breaking**), so the next release is 2.0.0.

### Added

- `cmake/SemVerDeps.cmake` and its function `semver_dependency`: a project requires its in-house libraries by a
  range of releases (`VERSION 2.1`: 2.1.0 or a later 2.x), and CMake fetches them at configure time, one copy of each
  per build tree. The first project of the tree that requires a library selects its release; the others check it,
  and the configure stops, naming both, when a requirement is not satisfied. The releases selected are kept in
  `deps.lock`, in the top-level source folder, for the next configures, which then need no network;
  `-DSEMVER_UPDATE=ALL` (or `=<name>`) selects the newest ones again. FLINT does not commit its `deps.lock`.
- The static library target `FLINT` (alias `FLINT::FLINT`), next to the list of sources `FLINT_SRC`. It is
  compiled once per build tree, with the flags of the top-level project, and only for the targets that link it; it
  carries its dependencies (ORION, OSLO, BLAS/LAPACK, and Cantera and OpenMP when enabled) and its module folder.
- The options `USE_OPENMP` and `USE_MPI` of FLINT, both off, declared before ORION and OSLO are added. OSLO's own
  options, on by default, turned OpenMP and MPI on in the whole build tree. A project that adds FLINT sets them
  first, and FLINT's then do nothing.
- `install.sh build --update-deps`, which selects the newest releases of ORION and OSLO again. `install.sh` keeps
  TecIO and SUNDIALS in `build-deps/` (`-DORION_TECIO_DIR`, `-DOSLO_SUNDIALS_DIR`), so that a new build folder does
  not build them again.
- This changelog, `CONTRIBUTING.md`, templates for issues (question, bug report, feature request) and pull requests,
  and the page Semantic Versioning of the Developer Guide. The README, `CONTRIBUTING.md` and the About page say where
  to ask for help: the issue tracker.
- An API Reference in the User Guide: the procedures and module variables that a solver uses, module by module, with
  their arguments, units and return codes.
- An optional argument `status` of `CEA_solve`: 0 when the equilibrium iteration converged, 1 when it did not.

### Changed

- **Breaking.** ORION and OSLO are no longer submodules of FLINT (`lib/ORION`, `lib/OSLO`), added from the folders
  `ORION_PATH` and `OSLO_PATH`: FLINT requires OSLO 1.1.0 or a later 1.x and ORION 2.1.0 or a later 2.x, fetched at
  configure time (see Added). The variables `ORION_PATH` and `OSLO_PATH` are removed;
  `-DFETCHCONTENT_SOURCE_DIR_ORION=<folder>` (and `_OSLO`) builds a local checkout instead, which is what
  `install.sh build --include-orion=<folder>` and `--include-oslo=<folder>` now pass. A project that adds ORION or
  OSLO itself, with `add_subdirectory`, and then FLINT, gets two targets of the same name and stops at configure: it
  requires them with `semver_dependency`, as FLINT does, before it makes any of them available (Installation, Using
  FLINT in a CMake project).
- **Breaking.** ORION 2 in place of ORION 1 (the submodule was at 1.7.0). ORION 2.0 removed
  `tec_read_points_multivars`, and FLINT reads its tables with `tec_read_structured_multiblock`: the same files, with
  the same values. A build tree that holds FLINT and a project that requires ORION 1 stops at configure. ORION 2.1 is
  the first release with `ORION_TECIO_DIR`, which `install.sh` and the GitHub workflow pass.
- The Cantera submodule (`lib/cantera`) is at the release v3.2.0. It was at a development commit of Cantera 3.3.0a1
  (ef09fb4, 20 commits after v3.2.0) that removes the legacy C library on which Cantera's Fortran interface is built.
  FLINT does not rebuild a Cantera already built in `lib/cantera`: clean it (`scons clean` in `lib/cantera`) to
  build v3.2.0.
- The test `test-thermo-real` is `thermo/thermo-real`, with its tables in `test/thermo/INPUT`: the area `thermo` of
  the test suite holds the real-fluid test, and the area and label `realfluid` are removed.
- The GitHub workflow runs on every pull request and on demand, no longer on every push to `main`: the commits of a
  branch are tested through its pull request. It configures twice: without TecIO, to select the releases of ORION
  and OSLO, then with TecIO, kept in a cache for each release of ORION.

### Removed

- `test-runiv`, the test of the universal gas constant added in 1.1.0.

### Fixed

- `CEA_solve` returns mass fractions in `rhoi_out` also when the iteration does not converge: those of the input
  state, with `T_out = T_in`. It returned the partial densities given, `rhoi_in`, with nothing to tell them from the
  mass fractions of a converged call (for instance for a state whose energy lies beyond the thermo tables).
  `test-CEA` checks both cases.
- A `CEA_solve` call that does not converge writes one `[WARNING] FLINT CEA_solve` line with the input temperature
  and density, to standard output and standard error. It printed two bare integers, the iteration limit and the
  point index (`50 1`), to standard output only, and only when the iteration limit was reached.
- Cantera is built with its own copies of fmt and yaml-cpp (`system_fmt=n system_yamlcpp=n`). It took the copies it
  found on the machine, also those of an active conda environment, which `libcantera.a` does not hold: every
  executable linked with Cantera failed with undefined `fmt::` and `YAML::` symbols.
- `cantera_fortran` is linked before `cantera`, which it calls. GNU ld searches a static library only where it is on
  the command line, so in the former order the Cantera symbols of the Fortran interface were left undefined.
- The targets that use Cantera wait for its build (`cantera_world`), when the configure builds it: a parallel build
  could compile FLINT's sources before Cantera's Fortran modules existed.
- The README and the installation guide give the option of `install.sh` that selects the compilers as
  `--compilers=gnu` (or `intel`). They gave `--compiler`, which `install.sh` refuses, so the command of the Quick Start
  stopped at once. The README no longer advises installing Cantera apart: `--use-cantera` builds it from `lib/cantera`.
- Using the Library names the procedures as the code does: `read_realfluid_transport` (it gave
  `read_realgas_transport`), `chemistry_source(rhoi, T, wdot)` after `Assign_Mechanism` (it gave
  `Chemistry_Source(T, rhoi, wdot)`), `rhs_native(nz, time, Z, F)`, and `ph2pT()` before `pT2h`. The real-fluid
  `thermo.dat` of Native Format has 8 variables, with $(\partial\rho/\partial h)_p$ fourth (it listed 7, so every
  index from the fourth on was wrong), and its table `h_tab2D` is built by `ph2pT()`, not read from the file.

## [1.1.0] - 2026-10-09

### Added

- A CTest suite: every test has an area label (`thermo`, `realfluid`, `chemistry`, `equilibrium`, `batch`) and a
  tier label (`quick`, `long`), and runs in `test/<area>/` with its fixtures. `ctest -L quick` runs the subset that a
  GitHub workflow runs on every push to `main` and every pull request.
- The test drivers `test-contract`, `test-ranges`, `test-tables`, `test-falloff`, `test-orders`, `test-stoich`,
  `test-andersen` and `test-runiv`, each with its fixtures and, where it compares with Cantera, references written
  by Cantera.
- An optional last block of `chemistry-info.txt`, `Reaction orders`: explicit reaction orders, fractional or on
  species that are not reactants, which the general procedure applies to the Arrhenius-type reactions as Cantera
  does. Orders on falloff reactions are refused. A file without the block, which comes from an older table writer,
  gives a WARNING; without it the rates are those of 1.0.0, bit for bit.
- A contract of the compiled mechanism routines: for each mechanism name with a compiled routine, its species slots
  (name and elemental composition) and its counts of Arrhenius, Troe and Lindemann reactions, in the module
  `FLINT_Lib_Chemistry_contract`, generated by `utils/mechanism_contract.py` from
  `src/lib/Lib_ChemMech/mechanism_contract.json`. `Assign_Mechanism` checks the loaded species and tables against
  it, and stops with the expected and the loaded lists when they differ: an INPUT folder with the name of a compiled
  routine and another species order produced source terms on the wrong species, without a message.
  `mechanism_contract.py --check` verifies that the committed module and the contract match the sources.
- A strict mode of `Assign_Mechanism`, off by default: a mechanism name without a compiled routine stops the program
  instead of falling back to the general procedure. It is turned on by the module flag `FLINT_strict_mechanism` or
  by the environment variable `FLINT_STRICT_MECHANISM` (`1`, `true`, `yes` or `on`, in any letter case; another
  value is reported).

### Changed

- The universal gas constant `Runiv` (and `Rr` of the CEA data) is the exact SI value, 8314.46261815324 J/(kmol K),
  that Cantera uses and with which the tables of `database/` are made. It was 8314.51: every pressure computed from
  a density was 5.7e-6 (relative) higher.
- `read_chemistry` reads the mechanism name as the whole first line of `chemistry-info.txt`. A name with a blank, a
  comma or a slash was cut at the separator ('Aramco 2.0' was read as 'Aramco').
- The loaders refuse a table with fewer zones than its reactions, and a table whose temperature step is not 1 K (with
  an `[ERROR]` line and a nonzero return value); the rate, transport and diffusion tables of a folder must also cover
  the range of its thermo tables. Every refusal is printed on the error unit as well as on standard output. Tables
  on a 1 K step with one zone per reaction, as the table writers produce them, are read as before.
- The WARNING of `Assign_Mechanism` for a name without a compiled routine is printed on the error unit as well.
- The species of `database/WD` are in the slot order of the routine WD (CH4 O2 CO2 H2O CO); they were in alphabetical
  order (see Fixed).
- The INPUT folders of the tests are in `database/`.
- ORION 1.7.0 and a newer OSLO, as submodules.

### Fixed

- The hand-written mechanism routines (WD, JLR, singh, coronetti, global-H2) read the rate at T kelvin from tables
  that do not start at 1 K; they read the rate at T + T0 - 1 for a table that starts at T0. `comp_ch_tabT` keeps its
  interface. The tables of `database/` start at 1 K, or are read by generated routines, so their results do not
  change.
- The general procedure raises the concentrations of Arrhenius reactions to their real stoichiometric coefficients;
  it rounded them to integers (H2 + 0.5 O2 <=> H2O ran with [O2]^1, 0.25 O2 with [O2]^0). The Troe and Lindemann
  reactions already used the real coefficients.
- With `database/WD` in alphabetical order the routine WD took CO as the oxidizer and wrote the source of O2 into the
  column of CO.
- WD-Andersen computes step 3 with the Andersen law, kf [CO2] [H2O]^0.5 [O2]^-0.25 (it used kf [CO2]^1.25), so steps
  2 and 3 reach the equilibrium of CO + 0.5 O2 <=> CO2. A step with a species of nonzero order at zero concentration
  has a zero rate, Cantera's convention, and no negative power of zero is evaluated: it stopped the DEBUG builds.
- Falloff reactions give finite rates where the tables vanish (with no third body, log10(0) gave a NaN broadening
  factor), and a reaction with k_c <= 0 is irreversible.
- `rhs_native` and `jac_native` return their bail-out value also for a temperature outside the rate tables (they
  read below the first row) and for a NaN temperature, now tested on its bit pattern: with `-ffast-math`, `isnan`
  was folded to false and a NaN indexed the tables at `-huge`.
- `jac_native` zeroes the source and Jacobian blocks before it calls the mechanism routine, and the routines define
  their whole INTENT(OUT) arguments: with species loaded after the slots of a routine, undefined values entered the
  energy row.
- The analytical Jacobian of ONERA-7 accepts species loaded after its seven slots (a shape mismatch, which
  `-fcheck=all` stops).
- The test drivers are linked as Fortran programs with TecIO on: with ifx the C++ link dropped their main program.
- `test-equilCXX` stops at a state without equilibrium and writes a reference only when it is complete.

## [1.0.0] - 2026-09-02

First release: ideal-gas and real-fluid (pressure-enthalpy tables) thermodynamics, transport properties (viscosity,
thermal conductivity, diffusion), finite-rate chemical kinetics with compiled routines for the mechanisms of
`database/` and a general procedure for the others, equilibrium chemistry, and 0-D reactors, with an optional
Cantera backend.
