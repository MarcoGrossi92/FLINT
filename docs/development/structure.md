# Code Structure

FLINT is organized as a modular Fortran library with a clear separation between:

* Core thermodynamic and chemistry kernels
* Generated mechanism-specific routines
* Equilibrium (CEA) solver
* Optional Cantera interface
* Test and validation programs
* Utilities and mechanism generation tools

The overall structure is shown below.

## Top-Level Layout

```
CMakeLists.txt
src/
lib/
database/
test/
bin/
utils/
docs/
cmake/
```

### Key Directories

| Directory               | Purpose                                             |
| ----------------------- | --------------------------------------------------- |
| `src/`                  | FLINT source code                                   |
| `src/lib/`              | Core library modules                                |
| `src/lib/Lib_ChemMech/` | Mechanism-specific explicit routines                |
| `src/test/`             | Test drivers, one folder per area                   |
| `lib/`                  | External submodules (OSLO, ORION, optional Cantera) |
| `database/`             | Raw chemical mechanism data (per-mechanism)          |
| `test/`                 | Test fixtures, references and outputs, per area     |
| `bin/`                  | Compiled test executables (`bin/test/`)             |
| `utils/`                | Mechanism generation tools                          |
| `cmake/`                | Build configuration modules                         |
| `docs/`                 | Documentation                                       |

---

## Core Library (`src/lib/`)

This directory contains the main FLINT computational kernels.

### Thermodynamics & Transport

```
Lib_ThermoTransport.f90
Load_ThermoTransport.f90
```

Responsible for:

* NASA polynomial evaluation
* Thermodynamic properties (`cp`, `h`, `s`, etc.)
* Transport property support
* Ideal-gas mixture handling

### Chemistry Kernel

```
Lib_Chemistry_data.f90
Lib_Chemistry_rhs.f90
Lib_Chemistry_wdot.f90
Lib_Chemistry_falloff.f90
Load_Chemistry.f90
```

Provides:

* Reaction data structures
* Source term computation (`wdot`)
* RHS evaluation for ODE integration
* Support for Arrhenius, Lindemann, and Troe formulations
* Mechanism loading from input files

These routines are mechanism-agnostic and operate on general chemistry data.

### Mechanism-Specific Explicit Routines

```
src/lib/Lib_ChemMech/
```

Contains dedicated Fortran source files such as:

```
WD.f90
ZK.f90
TSR-GP-24.f90
ecker.f90
...
```

These files implement:

* Hard-coded reaction kernels
* Optimized source term evaluation
* Mechanism-specific RHS routines

They are generated using the mechanism generation tool (see `utils/YTF.py`).

These routines provide:

* Maximum performance
* Production-level chemistry evaluation

### Chemical Equilibrium (CEA Solver)

```
Lib_CEA_data.f90
Lib_CEA_setup.f90
Lib_CEA_solver.f90
```

Implements:

* NASA CEA-based equilibrium solver
* Constant-volume (UV) equilibrium
* Species mass fraction update
* Equilibrium temperature calculation

This solver operates independently of the finite-rate chemistry kernel.

### Optional Cantera Interface

```
Load_Cantera.f90
```

Provides:

* Interface to Cantera routines
* Reference solution comparison
* Cross-validation capability

Cantera is optional and not required for production use.

## Test Programs (`src/test/`)

Test programs are separated from the core library, one folder per area:

```
src/test/CMakeLists.txt     builds the drivers into bin/test/ and registers the CTest tests with their labels
src/test/thermo/            test-runiv, test-thermo-ideal
src/test/realfluid/         test-thermo-real
src/test/chemistry/         test-tables, test-ranges, test-rhs-range, test-contract, test-inert, test-falloff,
                            test-orders, test-orders-warn, test-stoich, test-andersen, test-wdot
src/test/equilibrium/       test-CEA, test-equilCXX.cpp (Cantera C++: the equilibrium references)
src/test/batch/             test-batchF, test-batchCXX.cpp (Cantera C++: the batch references)
```

These validate:

* Thermodynamic properties (ideal-gas and real-fluid)
* Rate tables, the mechanism contract and the rate laws
* Batch reactor integration
* Equilibrium solver

Compiled executables from both Fortran and C++ tests are placed in `bin/test/`.

## Mechanism Database (`database/`)

The `database/` directory stores the raw data for each chemical mechanism supported by FLINT. Each mechanism has its own subfolder, named after the mechanism (e.g. `WD/`, `ZK/`, `Cross/`, `Pelucchi/`):

```
<Mechanism>/
    <Mechanism>.yaml
    chemistry-Arrhenius.dat
    chemistry-Troe.dat
    chemistry-info.txt
    composition.txt
    phase.txt
    thermo.dat
```

* The `<Mechanism>.yaml` file is the Cantera-format mechanism definition and is the input consumed by `utils/YTF.py` to generate the dedicated explicit routine in `src/lib/Lib_ChemMech/`.
* The remaining files (`chemistry-*.dat`, `composition.txt`, `phase.txt`, `thermo.dat`) are the plain-text data FLINT reads at runtime to build the general (non-explicit) chemistry and thermodynamic model for that mechanism.

Adding support for a new mechanism starts by adding a new subfolder here.

## Test Cases and Validation Data (`test/`)

The tests are grouped in five areas: `thermo`, `realfluid`, `chemistry`, `equilibrium`, `batch`. The drivers
are in `src/test/<area>/`, their data in `test/<area>/`, where CTest runs them:

```
test/
  realfluid/INPUT/          real-fluid tables
  chemistry/<topic>/        fixtures of the chemistry drivers (tables, ranges, orders, stoich, andersen, inert),
                            each written by the make_*.py next to it
  equilibrium/              reference/ (the Cantera equilibrium sweeps), eq-verification.py
  batch/                    cases.txt (the batch reactor cases), reference/ (the Cantera references),
                            element-standard-entropies.yaml (read by Cantera), batch-verification.py,
                            batch-performance.py
```

The outputs of the drivers are written next to their data and ignored by git: `batch/<case>/batch-*.dat`,
`batch/comp-batch-*.dat`, `equilibrium/<mechanism>/FLINT-CEA.txt`, `chemistry/wdot/<mechanism>/`.
See the [Testing](testing.md) page for the tests, their labels and the GitHub workflow.

---

## Mechanism Generation (`utils/`)

```
utils/YTF.py
```

This tool:

* Parses mechanism definitions
* Generates optimized Fortran source files
* Writes new modules into `Lib_ChemMech/`

The generated files expand FLINT’s set of dedicated explicit routines.

Mechanism generation is part of the development workflow and is documented in:

```
development/chemistry_generation.md
```

## External Dependencies (`lib/`)

```
lib/OSLO
lib/ORION
lib/cantera
```

These are managed as Git submodules.

* **OSLO / ORION**: Numerical infrastructure
* **Cantera**: Optional reference implementation

## Architectural Overview

FLINT follows a layered architecture:

```
Applications / Test Programs
          ↓
Mechanism-Specific Routines (Generated)
          ↓
General Chemistry Kernel
          ↓
Thermodynamic & Transport Layer
```

The equilibrium solver (CEA) operates as a parallel module using thermodynamic data.

## Design Principles

* Separation of data loading and computation
* Mechanism-agnostic core
* Optional reference backend (Cantera)
* Generated high-performance chemistry kernels
* Strict verification against reference implementations

---
