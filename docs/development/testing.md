# Testing Infrastructure

FLINT includes a suite of standalone Fortran programs used for:

- Numerical verification
- Cross-validation against Cantera
- Performance benchmarking
- Regression testing during development

These programs are compiled as part of the standard build process and are located in the `bin/test` directory.

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

## Categories of Tests

### 1. Thermodynamic Verification

- Compares specific heat and thermodynamic properties
- Native implementation vs Cantera
- Checks relative error and execution time
- Split into ideal-gas and real-fluid variants

Executables:
```

test-thermo-ideal
test-thermo-real

```

---

### 2. Chemical Source Terms

- Validates species production rates
- Explicit kernel vs Cantera net production rates

Executable:
```

test-wdot

```

---

### 3. Reactor Integration

- Constant-volume batch reactor
- Compares temperature evolution across:
  - General chemistry routines
  - Dedicated chemistry kernels
  - Cantera backend (if enabled)

Executable:
```

test-batchF

```

---

### 4. Chemical Equilibrium

- Validates CEA-based equilibrium solver
- Compares equilibrium temperature and selected species

Executable:
```

test-CEA

```

---

### 5. Contract and Table Unit Tests

Three drivers that need no Cantera and no fixture beyond `database/WD` and `test/tables/WD-100K`
(made by `test/tables/make_WD-100K.py`); each prints ` Verdict -> pass|fail` and exits 1 on failure:

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
                routine order accepted,
                swapped slots / wrong reaction count refused, appended inert species accepted,
                calibrated species with the same composition accepted, name read whole; in child
                processes: the strict fallback policy (FLINT_STRICT_MECHANISM) and the refusal channels
test-falloff    the Troe/Lindemann rates where the tables vanish (Pr = 0, k_inf = 0 give a
                zero rate, no NaN) and the k_c <= 0 convention for irreversible falloff
test-ranges     the temperature-grid contract of the tables: falloff tables on the Arrhenius grid, transport
                tables on the thermo grid, rate tables that do not cover the thermo grid (refusals), a
                rate table wider than the thermo grid (accepted, same rows as on the thermo grid)
test-inert      species appended after the slots of a compiled routine are inert on every path (direct
                call, rhs_native, analytical Jacobian, jac_native) with sentinel-filled outputs
test-rhs-range  rhs_native/jac_native bail out (F = -1, zero Jacobian) outside the RATE tables as
                they do outside the thermo tables (test/tables/WD-100K loaded on one grid, the
                rate range then narrowed in memory: a defence for tables set by another path)
test-orders     the general procedure with the optional 'Reaction orders' block reproduces Cantera's
                rates for JLR-frassoldati (yaml orders) and, without it, Cantera's law of the stoichiometric
                coefficients (the yaml orders removed), not the integer-rounded one of older versions
                (fixture test/orders/JLR-frassoldati: tables of a table writer, references embedded)
test-stoich     the general procedure reproduces Cantera's net production rates (to 1e-10 of the gross
                rates) for fractional stoichiometric coefficients without orders: Arrhenius reactants and
                products, three-body, Troe and Lindemann, 24-30 states each, and its net rates vanish at
                Cantera's equilibrium composition; integer control 2 H2 + O2 <=> 2 H2O; general has no
                analytical Jacobian (fixtures test/stoich/<name> made by test/stoich/make_stoich.py from
                constructed yaml mechanisms, tables of a table writer from the yaml thermo, references
                embedded)
test-andersen   the WD-Andersen routine (step 3 with the Andersen orders [CO2] [H2O]^0.5 [O2]^-0.25)
                reproduces Cantera's net production rates on the tables of a table writer (fixture
                test/andersen/WD-Andersen, references embedded), the zero rate at O2 = 0 and the
                zero-concentration convention of Coronetti (H2 = 0: finite, no divide-by-zero)
test-runiv      the universal gas constant is the exact SI value (8314.46261815324 J/(kmol K), the
                value of Cantera 3.0.1) in FLINT_Lib_Thermodynamic and in the CEA data, Ri_tab derives
                from it, the pressure p = rho R_mix T of a Cantera state of database/WD agrees with
                Cantera to 1e-9, and the compiled Frolov routine reproduces its (p/p_atm)^-1.15 law at
                the Cantera pressure to 1e-9 (5.7e-6 / 6.6e-6 off with the former 8314.51)
```

## Running the Test Suite

From the `test` directory:

```bash
./../bin/test/<test-name>
```

Each test reports:

* Key computed quantities
* Reference values (if applicable)
* Success/failure verdict

---

## Test Data Layout and Management

Test data is split across two directories:

* `database/` — raw mechanism data (Cantera YAML, thermodynamic and chemistry data files) used to build the library and to generate the mechanism-specific explicit routines via `utils/YTF.py`.
* `test/` — per-mechanism inputs and reference ("blessed") outputs used by the test executables at runtime.

Both directories are organized per mechanism, using the same mechanism name as the subfolder (e.g. `WD/`, `ZK/`, `Cross/`):

```
database/<Mechanism>/   # mechanism definition and raw data
test/<Mechanism>/       # INPUT/ + reference outputs for that mechanism
```

### Adding a New Mechanism to the Test Suite

1. Add the mechanism data under `database/<Mechanism>/` (YAML file plus supporting `.dat`/`.txt` files).
2. Generate the dedicated explicit routine with `utils/YTF.py <Mechanism>` and place the resulting source in `src/lib/Lib_ChemMech/`.
3. Create a matching `test/<Mechanism>/` folder with an `INPUT/` subfolder for any runtime input files.
4. Run the relevant test executables against the new mechanism and save the resulting outputs (e.g. `batch-explicit.dat`, `wdot-explicit.dat`, `eq-ref.txt`) as the reference values for future regression checks.
5. Keep `test/<Mechanism>/` and `database/<Mechanism>/` in sync — removing or renaming a mechanism should update both locations.

---

## Regression Strategy

Reference ("blessed") values are embedded in the test drivers.
A test fails if:

* The solution is not finite
* Relative error exceeds defined tolerance
* Unexpected numerical behavior is detected

This ensures that modifications to:

* Thermodynamic routines
* Chemistry kernels
* Solver infrastructure

do not silently alter validated behavior.

---

## Adding a New Test

To add a new regression test:

1. Create a standalone Fortran driver.
2. Load the required mechanism and data.
3. Define reference values or comparison logic.
4. Add success/failure criteria.
5. Register the executable in the CMake configuration.

Tests should:

* Be deterministic
* Avoid unnecessary I/O
* Use clear tolerances
* Focus on a single capability

---

## Continuous Integration (Optional)

When integrated into CI workflows, test executables can be run automatically after each build to ensure numerical stability across commits.




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
