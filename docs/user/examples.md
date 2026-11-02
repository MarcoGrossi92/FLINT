# Examples

This page outlines, in pseudo-code, two applications of FLINT with its native routines. [Quick Start](../getting-started/quick-start.md)
gives a complete program that you can build and run; the test drivers `src/test/batch/test-batchF.f90` and
`src/test/equilibrium/test-CEA.f90` are complete programs for the two applications below.

## Batch Reactor Integration

This example demonstrates the use of FLINT’s finite-rate chemical solver.

The procedure relies on the following steps:

- Load thermodynamic and chemistry data.
- Definition of initial state.
- Initialization of the solver.
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

  ! Initialize tolerances and state
  ! ...

  ! Initialize ODE solver (OSLO)
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  ! Initialize chemical mechanism specific routine, and its analytical Jacobian if it has one
  call Assign_Mechanism(mech_name)
  call set_analytical_jacobian(.true.)

  ! Run the solver (OSLO)
  call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)

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

! Initialize solver (CEA_solve also does it, at its first call)
call CEA_initialize_global()

! Initialize state values of T, rhoi
! ...

! Run the solver: equilibrium at the internal energy and volume of (T, rhoi)
call CEA_solve(T, rhoi, teq, y_eq, status)

! Output: equilibrium temperature (teq) and mass fractions (y_eq); status = 0 if converged
```

---