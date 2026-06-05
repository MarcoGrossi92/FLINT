# Examples Definition

This section provides procedural description and pseudo-code for several applications of FLINT. Only native routines are presented. 

## Batch Reactor Integration

This example demonstrates the use of FLINT’s finite-rate chemical solver.

The procedure relies on the following steps:

- Load thermodynamic and chemistry data.
- Definition of initial state.
- Inizialization of the solver.
- Run the solver.

**Pseudo-code**

```fortran
  use oslo
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  use FLINT_Lib_Chemistry_data
  use FLINT_Lib_Chemistry_rhs
  use FLINT_Load_chemistry

  ! Load data
  err = read_idealgas_thermo('mechanism-folder')
  err = read_chemistry( folder='mechanism-folder', mech_name=mech_name )

  ! Initialize tollerances and state
  ! ...

  ! Initialize ODE solver (OSLO)
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  ! Initialize chemical mechanism specific routine
  call Assign_Mechanism(mech_name)

  ! Run the solver (OSLO)
  call run_odesolver(neq,timein,timeout,Y,rhs_native,err)

  ! Output stored in Y (N species + temperature)

```

## Equilibrium Calculations

This example demonstrates the use of FLINT’s chemical equilibrium solver based on the NASA CEA algorithm for constant-volume systems.

The procedure relies on the following steps:

- Load thermodynamic data.
- Definition of initial state.
- Initialization of the solver.
- Run the solver.

**Pseudo-code**

```fortran
use FLINT_Lib_Thermodynamic
use FLINT_Load_ThermoTransport
use FLINT_CEA_setup
use FLINT_CEA_solver

! Load data
err = read_idealgas_thermo('mechanism-folder')

! Initialize solver
call CEA_initialize_global()

! Initialize state values of p, T, rhoi
! ...

! Run the solver
call CEA_solve(p, T, rhoi, teq, y_eq)

! Output stored in N species (y_eq) plus temperature (teq)
```

---