program test2
  use oslo
  use U_cantera_load
  use U_Lib_Thermodynamic
  use U_IO_Table
  use U_Lib_Chemistry, only: rhs_cantera, rhs_native, gas
  use U_Lib_ChemMech, only: Assign_Mechanism
  use U_IO_chemistry
  implicit none
  real(8)                    :: R, Tout, pout, pin, Tin
  real(8), allocatable       :: sp_Y(:)
  real(8), allocatable       :: Y(:)
  real(8), allocatable       :: RT(:), AT(:)
  integer                    :: err, neq, nstep, n
  real(8)                    :: timein, timeout, dt=0d0, tlim=0d0, time1, time2
  character(32)              :: solver, mech_name
  integer                    :: iopt(3)

  call load_phase(gas, 'coria.yaml')

  ! WD
  ! tlim = 0.01
  ! coria
  tlim = 0.5

  pin = 1.0d+5
  Tin = 1000

  nstep = 1000

  call Read_IdealGas_Properties()
  err = read_chemistry_file( mech_name=mech_name )
  if (err==0) then
    call Assign_Mechanism(mech_name)
  else
    write(*,*) 'Error during chemistry file loading'
    stop
  endif

  neq = nsc + 1
  allocate(Y(neq))
  allocate(sp_Y(nsc))

  sp_Y = 1d-20
  ! WD
  !sp_Y(1) = 0.2
  !sp_Y(2) = 0.8
  ! coria
  sp_Y(9) = 0.2
  sp_Y(4) = 0.8

  iopt = 0
  allocate(RT(neq),AT(neq))
  iopt(1)=100000
  RT(1:nsc)=1d-7
  RT(neq)=1d-7
  AT(1:nsc)=1d-7
  AT(neq)=1d-7
  solver = 'dvodef90'
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  dt = tlim/nstep

  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,err)
    R = Rtot(Y(1:neq-1))
    Tout = y(neq) !E2T(y(neq),Y(1:nsc)) !y(neq) !E2T(y(neq),Y(1:nsc)) !y(neq) !/ (sum(Y(1:neq-1))*R)
    pout = Tout*sum(Y(1:neq-1))*R
    !write(*,*) timeout, pout, Tout
  enddo
  call cpu_time(time2)

  print*, time2-time1
  write(*,*) timeout, pout, Tout

  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_cantera,err)
    R = Rtot(Y(1:neq-1))
    Tout = y(neq) !E2T(y(neq),Y(1:nsc)) !y(neq) !E2T(y(neq),Y(1:nsc)) !y(neq) !/ (sum(Y(1:neq-1))*R)
    pout = Tout*sum(Y(1:neq-1))*R
    !write(*,*) timeout, pout, Tout
  enddo
  call cpu_time(time2)

  print*, time2-time1
  write(*,*) timeout, pout, Tout

contains

  subroutine initialize()
    implicit none
    real(8) :: rho

    R = Rtot(sp_Y)
    rho = pin/(R*Tin)

    Y(1:nsc) = rho*sp_Y
    Y(neq) = Tin
    
    timein  = 0.D0
    timeout = 0.D0
  end subroutine

end program test2