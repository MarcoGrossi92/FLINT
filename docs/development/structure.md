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
| `src/test/`             | Test program sources                                |
| `test/`                 | Test cases and validation data                      |
| `bin/`                  | Compiled test executables                           |
| `lib/`                  | External submodules (OSLO, ORION, optional Cantera) |
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

Test programs are separated from the core library.

```
src/test/Fortran/
src/test/CXX/
```

### Fortran Tests

* `test-thermo.f90`
* `test-wdot.f90`
* `test-batchF.f90`
* `test-CEA.f90`

These validate:

* Thermodynamic properties
* Chemical source terms
* Batch reactor integration
* Equilibrium solver

### C++ Test

```
test-batchCXX.cpp
```

Used to generate Cantera reference batch-reactor solutions.

## Test Cases and Validation Data (`test/`)

The `test/` directory contains:

* Mechanism input files
* YAML files
* Thermodynamic data
* Output files
* Performance comparison data

Each mechanism has the structure:

```
<Mechanism>/
    INPUT/
    OUTPUT/
```

This allows systematic validation across multiple chemical mechanisms.

Python scripts:

```
verification.py
performance.py
```

are used for post-processing and benchmarking.

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
