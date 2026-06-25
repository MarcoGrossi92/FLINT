# Using the Library

FLINT is designed to be embedded inside a Fortran-based reacting-flow solver.
This section describes the typical workflow for:

* Thermodynamic and transport property evaluation (ideal gas and real fluid)
* Finite-rate chemistry
* Chemical equilibrium

For more practical examples on how to use FLINT refers to the [Examples](../examples/definition.md).

---

## Thermodynamic Properties

Both the ideal gas and real fluid models share the same module interface. The model is selected implicitly by calling the appropriate load routine; all downstream property functions operate on whichever tables are in memory.

**Required Modules**

```fortran
use FLINT_Lib_Thermodynamic
use FLINT_Load_ThermoTransport
```

### Ideal Gas

**Loading**

```fortran
integer :: err
err = read_idealgas_thermo("folder")
err = read_idealgas_transport("folder") ! optional
```

This reads species names, molecular weights, temperature-indexed thermodynamic tables ($c_p$, $h$, $s$), and — optionally — transport tables ($\mu$, $k$). Thermodynamics is also required by both finite-rate chemistry and equilibrium calculations.

### Real Fluid

For single-component fluids at high pressure or near the critical point, FLINT provides a real fluid model backed by uniform 2D $(p, h)$ lookup tables.

**Loading**

```fortran
integer :: err
err = read_realfluid_thermo("folder")
err = read_realgas_transport("folder") ! optional
```

`read_realfluid_thermo` reads the species name from `phase.txt` and populates a $(p, h)$ grid with density, temperature, entropy, speed of sound, $c_p$, and two partial derivatives. `read_realgas_transport` populates the matching viscosity and thermal conductivity grid.

**Error codes**

Both routines signal errors through their integer return value:

| Code | `read_realfluid_thermo` | `read_realgas_transport` |
|------|------------------------|--------------------------|
| 0 | No error | No error |
| 1 | `phase.txt` not found | Transport file not found |
| 2 | Error reading `phase.txt` | Error reading transport data |
| 3 | Thermo file not found | More than one block |
| 4 | Error reading thermo file | Mesh size mismatch with thermo |
| 5 | More than one block in thermo | — |

**Property evaluation**

Properties at a given state $(p, h)$ are retrieved by bilinear interpolation using `ph2vars`:

```fortran
real(8) :: rho, T_fluid, spd

rho      = ph2vars(p, h, rho_tab)
T_fluid  = ph2vars(p, h, T_tab)
spd      = ph2vars(p, h, sound_tab)
```

Available tables exposed by `FLINT_Lib_Thermodynamic`:

| Table | Property | Unit |
|-------|----------|------|
| `rho_tab` | Density | kg/m³ |
| `T_tab` | Temperature | K |
| `dT_tab` | $(\partial\rho/\partial T)_p$ | kg/(m³·K) |
| `hT_tab` | $c_p = (\partial h/\partial T)_p$ | J/(kg·K) |
| `s_tab2D` | Specific entropy | J/(kg·K) |
| `rp_tab` | $(\partial\rho/\partial p)_h$ | kg/(m³·Pa) |
| `sound_tab` | Speed of sound | m/s |
| `mi_tab2D` | Dynamic viscosity | Pa·s |
| `k_tab2D` | Thermal conductivity | W/(m·K) |

To initialise from primitive variables $(p, T)$ rather than $(p, h)$, convert with `pT2h` before querying the tables:

```fortran
real(8) :: h_init
h_init = pT2h(p, T)
```

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

  - hard-coded Fortran kernels  
  - maximum performance  
  - recommended for production  

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
