# Testing Infrastructure

FLINT includes a suite of standalone Fortran programs used for:

- Numerical verification
- Cross-validation against Cantera
- Performance benchmarking
- Regression testing during development

They are compiled with the library (into `bin/test/`) and registered as CTest tests.

---

## Testing Philosophy

FLINT testing follows three principles:

1. **Numerical consistency**  
   Native routines and Cantera-interface ones must reproduce reference solutions within defined tolerances.

2. **Backend consistency**  
   Dedicated chemistry kernels, general routines, and optional Cantera interfaces must produce equivalent physical results.

3. **Regression stability**  
   Changes to the codebase must not alter validated results unexpectedly.

---

## Organization

The tests are grouped in five areas. An area has the same name in the sources, in the data and in the tests:

```
src/test/<area>/test-<name>.f90    the driver, built as bin/test/test-<name>
test/<area>/                       its fixtures and references (committed), its outputs (ignored by git)
<area>/<name>                      the CTest test, run in test/<area>/
```

| Area | What it covers | Tests |
|---|---|---|
| `thermo` | ideal-gas thermodynamics | `runiv` |
| `realfluid` | real-fluid (p, h) tables | `thermo-real` |
| `chemistry` | rate tables, mechanism contract, rate laws against Cantera references | `tables`, `ranges`, `rhs-range`, `contract`, `inert`, `falloff`, `orders`, `orders-warn`, `stoich`, `andersen` |
| `equilibrium` | the CEA solver against Cantera values | `CEA` |
| `batch` | constant-volume reactor against Cantera references, one test per case | `WD`, `Troyes`, `Ecker`, `Cross`, `Smooke`, `CORIA`, `TSR-CDF-13`, `Pelucchi`, `ZK`, `TSR-GP-24`, `TSR-Rich-31`, `Gerlinger` |

```
src/test/
  CMakeLists.txt            the drivers and the registration of every test with its labels
  thermo/                   test-runiv, test-thermo-ideal
  realfluid/                test-thermo-real
  chemistry/                test-tables, test-ranges, test-rhs-range, test-contract, test-inert, test-falloff,
                            test-orders, test-orders-warn, test-stoich, test-andersen, test-wdot
  equilibrium/              test-CEA, test-equilCXX (Cantera C++, USE_CANTERA)
  batch/                    test-batchF, test-batchCXX (Cantera C++, USE_CANTERA)
test/
  realfluid/INPUT/          real-fluid tables
  chemistry/tables/ ranges/ orders/ stoich/ andersen/ inert/
                            one fixture folder per topic, each with the make_*.py that wrote it
  equilibrium/              reference/<mechanism>.dat (the Cantera sweeps), eq-verification.py
  batch/                    cases.txt, reference/<case>.dat, element-standard-entropies.yaml,
                            batch-verification.py, batch-performance.py
```

test-thermo-ideal and test-wdot are Cantera benchmarks without a verdict: they are built, they are not tests.

### Labels

Every test has two labels:

* its **area**: `thermo`, `realfluid`, `chemistry`, `equilibrium`, `batch`
* its **tier**: `quick` (the subset the GitHub workflow runs, about 30 s) or `long` (the rest)

Every test is `quick` except the batch cases Ecker, Cross, Smooke, CORIA, TSR-CDF-13, ZK and TSR-Rich-31
(10-20 s each, most of it reading their text tables). The quick batch cases cover a global mechanism with
orders (WD), tables in TecIO's `.szplt` (Troyes, TSR-GP-24), Troe (TSR-GP-24), Troe and Lindemann (Pelucchi)
and H2/air (Gerlinger). The tier of a batch case is the `BATCH_QUICK` list of `src/test/CMakeLists.txt`.

---

## Running the Test Suite

From the build directory, after the build:

```bash
ctest -j8                     # every test
ctest -L quick                # what the GitHub workflow runs
ctest -L chemistry            # one area
ctest -L batch -LE long       # an area, quick tier only
ctest -R batch/ZK -V          # one test, with its output
```

Without `USE_TECIO`, the tests that read TecIO `.szplt` tables (batch Troyes, TSR-GP-24, TSR-Rich-31 and
equilibrium/CEA) are disabled.

A driver can also be run by hand from its area folder, for example:

```bash
cd test/chemistry && ../../bin/test/test-orders
cd test/batch && ../../bin/test/test-batchF check ZK
```

Each test reports:

* Key computed quantities
* Reference values (if applicable)
* Success/failure verdict (exit code 1 on failure)

---

## Continuous Integration

`.github/workflows/tests.yml` builds FLINT on Ubuntu (gfortran, `USE_TECIO=ON`, no Cantera, no SUNDIALS) and runs:

* `ctest -L quick` on every push to `main` and every pull request (changes to the documentation only are skipped)
* every test once a week (Monday 03:00 UTC)
* on demand (Actions > Tests > Run workflow) the tests of a label regex: `quick`, `batch`, `.` for every test

TecIO, built by ORION at configure time, is cached for each ORION commit.

---

## Areas

### thermo

```
test-runiv      the universal gas constant is the exact SI value (8314.46261815324 J/(kmol K), the
                value of Cantera 3.0.1) in FLINT_Lib_Thermodynamic and in the CEA data, Ri_tab derives
                from it, the pressure p = rho R_mix T of a Cantera state of database/WD agrees with
                Cantera to 1e-9, and the compiled Frolov routine reproduces its (p/p_atm)^-1.15 law at
                the Cantera pressure to 1e-9 (5.7e-6 / 6.6e-6 off with the former 8314.51)
```

`test-thermo-ideal` (Cantera): specific heat of FLINT and of Cantera over 1..2000 K, times and relative difference.

### realfluid

```
test-thermo-real  the (p, h) tables of test/realfluid/INPUT: interpolation, and the inverse (p, T) -> h
                  of ph2pT, whose round trip T -> h -> T must stay within 0.5 K (0.16 K now)
```

### chemistry

The drivers need no Cantera; the fixtures are in `test/chemistry/<topic>/`, each written by the `make_*.py`
next to it (run from the repository root):

```
test-tables     rate tables are indexed by temperature whatever their first row:
                f_kf/f_kb return row T on database/WD and on the 100..400 K copy, and every
                hand-written routine gives bit-identical omegadot with tables starting at 50, 100,
                300 and 799 K (positive control: the assumed-shape accessor of FLINT <= 2223136);
                the public comp_ch_tabT equals f_kf/f_kb on the 1 K, the 100..400 K and the
                synthetic 100 K tables; the analytical Jacobians define their whole block with
                species after the slots
test-contract   the mechanism contract check of Assign_Mechanism (and, with python3 on the PATH, that
                src/lib/Lib_Chemistry_contract.f90 is the output of utils/mechanism_contract.py):
                routine order accepted, a hook before the tables checked at the first call,
                swapped slots / wrong reaction count refused, appended inert species accepted,
                calibrated species with the same composition accepted, name read whole; in child
                processes: the strict fallback policy (FLINT_STRICT_MECHANISM) and the refusal channels
test-falloff    the Troe/Lindemann rates where the tables vanish (Pr = 0, k_inf = 0 give a
                zero rate, no NaN) and the k_c <= 0 convention for irreversible falloff
test-ranges     the temperature-grid contract of the tables: falloff tables on the Arrhenius grid, transport
                and diffusion tables on the thermo grid or starting below it (accepted, same values as on
                the thermo grid) and starting above it (refused), rate tables that do not cover the thermo
                grid (refusals), a rate table wider than the thermo grid (accepted, same rows as on the
                thermo grid), rate tables with a later zone on another grid or an interior row missing
                (refusals), a falloff-Troe table with a NaN F_cent (refused, also under -ffast-math),
                thermo, transport and diffusion tables with an interior row missing (refusals)
test-inert      species appended after the slots of a compiled routine are inert on every path (direct
                call, rhs_native, analytical Jacobian, jac_native) with sentinel-filled outputs
test-rhs-range  rhs_native/jac_native bail out (F = -1, zero Jacobian) outside the RATE tables as
                they do outside the thermo tables (test/chemistry/tables/WD-100K loaded on one grid, the
                rate range then narrowed in memory: a defence for tables set by another path)
test-orders     the general procedure with the 'Reaction orders' block reproduces Cantera's
                rates for JLR-frassoldati (yaml orders) and, without it, Cantera's law of the stoichiometric
                coefficients (the yaml orders removed), not the integer-rounded one of older versions;
                the general procedure warns about a file without the block, not about a file with it;
                a block with a negative row count is refused
                (fixture test/chemistry/orders/JLR-frassoldati: tables of a table writer, references embedded)
test-orders-warn the WARNING of the general procedure for a chemistry-info.txt without the 'Reaction
                orders' block: none for a block with zero rows (same omegadot, bit for bit, as without the
                block), one when the general procedure is selected before or after the tables are loaded,
                none for a hooked name; once per load, on standard output and on the error unit (child
                process). The files without the block are copies made at run time (orders/noblock from
                test/chemistry/orders/JLR-frassoldati/chemistry-info-noblock.txt, orders/wd-noblock from database/WD):
                every chemistry-info.txt of database/ and of the fixtures ends with the block (0 rows when
                the mechanism has no orders), as a table writer writes it
test-stoich     the general procedure reproduces Cantera's net production rates (to 1e-10 of the gross
                rates) for fractional stoichiometric coefficients without orders: Arrhenius reactants and
                products, three-body, Troe and Lindemann, 24-30 states each, and its net rates vanish at
                Cantera's equilibrium composition; integer control 2 H2 + O2 <=> 2 H2O; general has no
                analytical Jacobian (fixtures test/chemistry/stoich/<name> made by test/chemistry/stoich/make_stoich.py from
                constructed yaml mechanisms, tables of a table writer from the yaml thermo, references
                embedded)
test-andersen   the WD-Andersen routine (step 3 with the Andersen orders [CO2] [H2O]^0.5 [O2]^-0.25)
                reproduces Cantera's net production rates on the tables of a table writer (fixture
                test/chemistry/andersen/WD-Andersen, references embedded), the zero rate at O2 = 0 and the
                zero-concentration convention of Coronetti (H2 = 0: finite, no divide-by-zero)
```

`test-wdot` (Cantera): FLINT's and Cantera's net production rates over 100..2000 K for the mechanisms,
written to `test/chemistry/wdot/<mechanism>/`.

### equilibrium

```
test-CEA        equilibrium at constant internal energy and volume of the CEA solver for the species of
                WD, ZK, TSR-GP-24 (1000 O2/CH4 mixture ratios from 0.01 to 100, 1000 K, 3.25 kg/m3) and
                Ecker (1000 pressures from 1e-5 to 100 bar, 3000 K): every equilibrium temperature of the
                sweep within 2e-3 of the Cantera reference test/equilibrium/reference/<mechanism>.dat
                (8.4e-4 at most now), and at one state per mechanism the temperature and a species mass
                fraction within 1 % of the Cantera values in the driver
```

The sweeps are written to `test/equilibrium/<mechanism>/FLINT-CEA.txt` and plotted against the references by
`eq-verification.py` (run from `test/equilibrium`). The references are written once by `test-equilCXX`
(Cantera C++, `equilibrate("UV")` from the same states) and committed; regenerate them, after a change to the
sweeps of test-CEA (also in test-equilCXX.cpp), to a mechanism or to Cantera, with
`cmake --build . --target equilibrium-reference` in a build configured with `USE_CANTERA`.

### batch

- Constant-volume batch reactor, the cases of `test/batch/cases.txt` (mechanism, Cantera yaml, end time,
  initial pressure, temperature and mass fractions; one line per case, read by both drivers and by CMake)
- Compares temperature evolution across:
  - General chemistry routines
  - Dedicated chemistry kernels
  - Cantera backend (if enabled)
- `test-batchF check <case>` (the CTest test `batch/<case>`) needs no Cantera: it compares the dedicated
  routine and the general procedure with the Cantera reference `test/batch/reference/<case>.dat`, stored in
  the repository, and the general procedure with the dedicated routine

Executables, run in `test/batch/`:
```
test-batchF      [verification | performance | check <case>]   (no argument: asks for the mode)
test-batchCXX    [verification | performance | --reference [case...]]   (Cantera C++, USE_CANTERA)
```

`verification` writes `test/batch/<case>/batch-<backend>.dat` (plotted by `batch-verification.py`),
`performance` the times in `test/batch/comp-batch-<backend>.dat` (`batch-performance.py`).
With Cantera, both drivers load the Cantera phase in `test/batch`, where Cantera reads
`element-standard-entropies.yaml` (the standard entropies of the elements; Cantera searches the working
directory first, and leaves the entropies unknown without the file).

Acceptance of `test-batchF check` (FLINT at RT = AT = 1e-7, reference at rtol 1e-10, atol 1e-15; the
tolerances are at the top of `test-batchF.f90`):

| Check | Tolerance |
|---|---|
| final temperature, relative | 5e-4 |
| time of half the temperature rise, relative | 2e-3 |
| mean \|T - T_ref\| over the run / temperature rise | 1e-3 |
| max \|T_general - T_explicit\| / temperature rise | 1e-5 |

#### Batch reactor references

The files `test/batch/reference/<case>.dat` are written once with Cantera and committed; CTest runs Fortran only.
Regenerate them after a change to a line of `test/batch/cases.txt`, to a mechanism of `database/` or to Cantera,
with a build configured with `USE_CANTERA`:

```bash
cmake --build . --target batch-reference                            # every case
cd ../test/batch && ../../bin/test/test-batchCXX --reference ZK Gerlinger   # some cases
```

The references are integrated at rtol 1e-10, atol 1e-15: with an absolute tolerance of 1e-7 on the mass fractions,
CVODE does not resolve the radicals that start at 0 (Gerlinger, 1200 K, does not ignite in 2e-4 s).

---

## Building with Cantera and regenerating the figures

The figures of the [verification page](../examples/verification.md) need the Cantera branches: Cantera's own
reactor and equilibrium (C++), and FLINT's integrator with the Cantera source terms (Fortran interface).
A second build with Cantera (it needs SUNDIALS, built by OSLO) keeps its drivers apart from the default build:

```bash
cmake -S . -B build-cantera -DUSE_CANTERA=ON -DUSE_SUNDIALS=ON -DUSE_TECIO=ON -DUSE_MPI=OFF -DUSE_OPENMP=OFF \
      -DFLINT_TEST_BINDIR=$PWD/bin/test-cantera
cmake --build build-cantera -j8

cd test/batch                                   # batch reactor: 4 datasets per case
../../bin/test-cantera/test-batchCXX verification     # Cantera            -> <case>/batch-CXX.dat
../../bin/test-cantera/test-batchF verification       # FLINT Cantera, FLINT Explicit, FLINT General
python3 batch-verification.py                   # docs/examples/images/<case>.svg

cd ../equilibrium                               # equilibrium: FLINT against the references
../../bin/test-cantera/test-CEA
python3 eq-verification.py                      # docs/examples/images/<mechanism>-eq.svg
```

---

## Adding a New Mechanism to the Test Suite

1. Add the mechanism data under `database/<Mechanism>/` (YAML file plus supporting `.dat`/`.txt` files).
2. Generate the dedicated explicit routine with `utils/YTF.py <Mechanism>` and place the resulting source in `src/lib/Lib_ChemMech/`.
3. Add a line for the mechanism to `test/batch/cases.txt` (CMake then adds the CTest test `batch/<Mechanism>`,
   in the `long` tier unless it is added to `BATCH_QUICK` in `src/test/CMakeLists.txt`).
4. Write its Cantera reference with `test-batchCXX --reference <Mechanism>` (from `test/batch/`) and commit
   `test/batch/reference/<Mechanism>.dat`; check it with `ctest -R batch/<Mechanism> -V`.

---

## Adding a New Test

1. Write the driver `src/test/<area>/test-<name>.f90` (every driver of an area folder is built). It runs in
   `test/<area>/`: the database is `../../database/`, its fixtures are in `test/<area>/`.
2. Define reference values or comparison logic, with clear tolerances; print ` Verdict -> pass|fail` and exit
   with code 1 on failure (`stop 1`).
3. Register it in `src/test/CMakeLists.txt` with its area and tier:

   ```cmake
   flint_add_test(<area> <name> quick test-<name>)
   ```

4. Commit its fixtures, with the script that writes them; add the files it writes at run time to `test/.gitignore`.
   Two tests that write the same files need a common `RESOURCE_LOCK`.

Tests should:

* Be deterministic
* Avoid unnecessary I/O
* Use clear tolerances
* Focus on a single capability

---

## Regression Strategy

Reference ("blessed") values are embedded in the test drivers, or stored next to their fixtures
(`reference.txt` of test/chemistry/orders, stoich, andersen; `test/batch/reference/` of the batch reactor).
A test fails (exit code 1) if:

* The solution is not finite
* Relative error exceeds defined tolerance
* Unexpected numerical behavior is detected

This ensures that modifications to:

* Thermodynamic routines
* Chemistry kernels
* Solver infrastructure

do not silently alter validated behavior.

<!-- # Testing & Verification

FLINT capabilities are verified through a suite of small Fortran programs compiled and executed as part of the test infrastructure.

The test suite focuses on numerical verification of thermodynamic properties, chemical source terms, reactor integration, and chemical equilibrium calculations.

Where available, Cantera is used as a reference implementation to assess numerical consistency.

These tests are primarily intended for developers and continuous integration, and are not meant to serve as end-user examples.


---

## Thermodynamic Properties

**Purpose**

Benchmark and verify the computation of the specific heat at constant pressure (`cp`) for ideal gas mixtures.

**Description**

This test compares:

* FLINT native thermodynamic routines
* Cantera thermodynamic routines

The benchmark evaluates both numerical agreement and computational performance over a wide temperature range.

**Simulation Modes**

The test must be executed from the `./test` directory, where the chemical mechanism data are located.
From a shell, run:
```
./../bin/test/test-thermo
```

**Test Procedure**

1. Load ideal-gas thermodynamic data.
2. Load the same mechanism into Cantera from.
3. Define a fixed gas mixture.
4. Loop over temperatures.
5. Compute thermodynamic properties.
6. Measure:
   * Wall-clock CPU time
   * Relative error between FLINT and Cantera `cp`

**Output**

* Execution time for FLINT
* Execution time for Cantera
* Relative error (%)

**Verification Criteria**

* Relative error typically below machine precision
* FLINT expected to outperform Cantera in raw throughput

---

## Chemical Source Terms

**Purpose**

Validate and benchmark species production rates (`\dot{\omega}`) for mechanisms including third-body reactions.

**Description**

This test compares:

* FLINT explicit chemistry kernel
* Cantera net production rates

Multiple chemical mechanisms are evaluated.

**Tested Mechanisms**

* Westbrook & Dryer
* Troyes
* Ecker
* Cross
* Pelucchi
* Smooke
* CORIA-CNRS
* TSR-CDF-13
* TSR-GP-24
* TSR-Rich-31

Each mechanism is tested independently using identical initial compositions.

**Simulation Modes**

The test must be executed from the `./test` directory, where the chemical mechanism data are located.
From a shell, run:
```
./../bin/test/test-wdot
```

**Test Procedure**

1. Load thermodynamics and chemistry data.
2. Loop over temperature range
3. Compute:

   * `Chemistry_Source` (FLINT)
   * `getNetProductionRates` (Cantera)
4. Store results for post-processing.

**Output**

For each mechanism:

* `OUTPUT/wdot-explicit.dat` — FLINT results
* `OUTPUT/wdot-cantera.dat` — Cantera results (if enabled)

Each file contains:
```
T  wdot_1  wdot_2  ...  wdot_ns
```

**Verification Criteria**

* Pointwise agreement between FLINT and Cantera

---

## Batch Reactor Integration

**Purpose**

Verify and benchmark time integration of reacting systems using FLINT’s ODE solvers.

**Description**

This test integrates a constant-volume batch reactor and compares:

* FLINT native RHS
* FLINT–Cantera RHS
* FLINT general (non-coded) chemistry RHS

Both accuracy and performance are evaluated.

**Tested Mechanisms**

* Westbrook & Dryer
* Troyes
* Ecker
* Cross
* Pelucchi
* Smooke
* Zhukov & Kong
* CORIA-CNRS
* TSR-CDF-13
* TSR-GP-24
* TSR-Rich-31

**Simulation Modes**

The test must be executed from the `./test` directory, where the chemical mechanism data are located.
From a shell, run:
```
./../bin/test/test-batchF
```

and select the desired simulation mode:

1. **Verification mode**

   * Many time steps
   * Accuracy-focused

2. **Performance mode**

   * Single time step
   * Timing-focused

**Output**

For each mechanism, the following files are written in the `./test/<mechanism>/OUTPUT` directory:

* `batch-general.dat`
* `batch-explicit.dat`
* `batch-cantera.dat` (if enabled)

Each file contains:
```
time  temperature
```

Performance summary files are written in `./test`:

* `comp-batch-general.dat`
* `comp-batch-explicit.dat`
* `comp-batch-canteraFor.dat`

**Verification Criteria**

* Consistent temperature evolution
* Agreement with Cantera within solver tolerances
* Expected speed-up from coded mechanisms

For comparison with Cantera, the C++ test program must be executed to generate the reference results:
```
./../bin/test/test-batchCXX
```

## Chemical Equilibrium

**Purpose**

Validate the FLINT equilibrium solver against Cantera results.

**Description**

This test computes chemical equilibrium at constant volume and compares:

* Equilibrium temperature
* Selected species mass fractions

**Tested Mechanisms**

* Westbrook & Dryer
* Zhukov & Kong
* TSR-GP-24
* Ecker

Note that mechanisms are loaded just to have a set of species.

**Simulation Modes**

The test must be executed from the `./test` directory, where the chemical mechanism data are located.
From a shell, run:
```
./../bin/test/test-CEA
```

**Test Procedure**

1. Load thermodynamic data.
2. Define initial composition.
3. Solve equilibrium.
4. Compare results with precomputed Cantera reference values.

**Verification Criteria**

* Equilibrium temperature
* Key species mass fraction (e.g. CO, OH)

A test is marked success if:

* Solution is finite
* Relative error < 1% compared to Cantera

---

## Summary

The FLINT test suite provides:

* Numerical verification against Cantera
* Performance benchmarks
* Validation across multiple chemical mechanisms

Together, these tests ensure correctness, robustness, and high performance of the FLINT library. -->
