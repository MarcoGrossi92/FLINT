# Using the Library

FLINT is designed to be embedded inside a Fortran-based reacting-flow solver.
This section describes the typical workflow for:

* Thermodynamic property evaluation
* Finite-rate chemistry
* Chemical equilibrium

For more practical examples on how to use FLINT refers to the [Examples](../examples/definition.md).

---

## Thermodynamic Properties

**Required Modules**

```fortran
use FLINT_Lib_Thermodynamic
use FLINT_Load_ThermoTransport
```

**Initialization**

```fortran
integer :: err
err = read_idealgas_thermo("folder")
err = read_idealgas_transport("folder") ! optional
```

This loads:

* Species names
* Molecular weights
* Thermodynamic properties tables
* Transport properties tables (optional)

**Property Evaluation**

Once loaded, mixture properties can be evaluated using the provided routines.

Typical operations include:

* `cp` (specific heat at constant pressure)
* `h` (enthalpy)
* `s` (entropy)
* Gas constant and mixture molecular weight

Thermodynamics is required by both finite-rate chemistry and equilibrium calculations.

### Source Files

- `src/lib/Load_ThermoTransport.f90`
- `src/lib/Lib_ThermoTransport.f90`

---

## Finite-Rate Chemistry

Finite-rate chemistry requires both thermodynamic and reaction data.

**Required Modules**

```fortran
use FLINT_Load_Chemistry
use FLINT_Lib_Chemistry_rhs
use FLINT_Lib_Chemistry_wdot
```

**Loading Chemistry Data**

```fortran
err = read_chemistry( "folder", mech_name )
```

This reads:

* Mechanism name (`mech_name`)
* Reaction definitions
* Arrhenius tables
* Falloff tables (Lindemann, Troe)
* Third-body efficiencies

**Computing Species Source Terms**

To compute net production rates:

```fortran
call Chemistry_Source(T, rhoi, wdot)
```

where:

* `T` = temperature
* `rhoi` = species partial densities
* `wdot` = species production rates

This subroutine is built upon an abstract interface and a pointer concretization. At runtime, the pointer procedure is assigned to one of the subroutine realizations. These are divided into:

- a general routine:  

  -- mechanism loaded at runtime
  -- flexible
  -- slightly slower

- several explicit routines:

        - Hard-coded Fortran kernels
        - Maximum performance
        - Recommended for production

Explicit routines are generated using the mechanism generation tool. See [Chemistry Database](database.md) and [Chemistry Generation](chemistry_routines.md) for more details.

**RHS Evaluation for ODE Integration**

For batch reactor simulations:

```fortran
call Chemistry_RHS(...)
```

This routine computes the full time derivative vector composed of species and temperature evolution.

### Source Files

- `src/lib/Lib_ChemMech/`
- `src/lib/Load_Chemistry.f90`
- `src/lib/Lib_Chemistry_data.f90`
- `src/lib/Lib_Chemistry_wdot.f90`
- `src/lib/Lib_Chemistry_rhs.f90`

---

## Chemical Equilibrium (CEA Solver)

FLINT provides a NASA-CEA-based equilibrium solver.

**Required Modules**

```fortran
use FLINT_CEA_setup
use FLINT_CEA_solver
```

Thermodynamic data must be loaded before using the equilibrium solver.

**Initialization**

```fortran
call CEA_initialize_global()
```

**Solving Equilibrium**

```fortran
call CEA_solve(T_initial, rhoi, T_eq, y_eq)
```

Inputs:

* `T_initial` — initial temperature
* `rhoi` — initial partial densities

Outputs:

* `T_eq` — equilibrium temperature
* `y_eq` — equilibrium mass fractions

The solver performs constant-volume (UV) equilibrium.

### Source Files

- `src/lib/Lib_CEA_data.f90`
- `src/lib/Lib_CEA_setup.f90`
- `src/lib/Lib_CEA_solver.f90`

---

## Optional Cantera Interface

FLINT can interface with Cantera for validation purposes.

Cantera is:

* Optional
* Not required for production runs
* Used for verification and benchmarking

When enabled, Cantera routines can be used to compute:

* Thermodynamic properties
* Net production rates
* Reactor integration results

The FLINT native implementation is the intended production backend.

See the [Cantera documentation](https://cantera.org/index.html) about Fortran interface for more details on this topic. 
