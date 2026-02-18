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

Executable:
```

test-thermo

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
