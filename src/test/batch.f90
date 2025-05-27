program batch_reactor
  use cantera
  implicit none

  type(phase_t) :: gas
  type(reactor_t) :: r
  type(sim1d_t) :: sim
  double precision :: t, dt, tfinal
  integer :: i, nsp
  double precision, allocatable :: x(:)
  character(len=20) :: name

  ! Load the gas phase (update the path/phase name if needed)
  gas = importPhase('h2o2.yaml', 'ohmech')

  ! Set initial state
  call setState_TPX(gas, 1000.0d0, 101325.0d0, 'H2:2, O2:1, AR:4')

  ! Create reactor with the gas
  r = reactor(gas)

  ! Create an integrator (1D time integration)
  sim = sim1d()
  call addReactor(sim, r)

  ! Set simulation time parameters
  t = 0.0d0
  tfinal = 1.0d-3  ! seconds
  dt = 1.0d-5      ! time step

  ! Allocate species mole fraction array
  nsp = nSpecies(gas)
  allocate(x(nsp))

  write(*,'(a)') '  time [s]     T [K]     X(H2)     X(O2)     X(H2O)'

  do while (t < tfinal)
    call advance(sim, t + dt)
    t = time(sim)

    call moleFractions(gas, x)

    ! Print time, temperature, and some species (if available)
    call getSpeciesIndex(gas, 'H2', i); write(*,'(f10.6,1x,f8.2,1x,f8.4)', advance='no') t, temperature(gas), x(i)
    call getSpeciesIndex(gas, 'O2', i); write(*,'(1x,f8.4)', advance='no') x(i)
    call getSpeciesIndex(gas, 'H2O', i); write(*,'(1x,f8.4)') x(i)
  end do

  deallocate(x)

end program batch_reactor
