program batch
  use oslo
  use U_cantera_load
  use U_Lib_Thermodynamic
  use U_IO_Table
  use U_Lib_Chemistry, only: rhs_cantera, rhs_native, gas, Assign_Mechanism
  use U_IO_chemistry
  implicit none
  type(file_ini)             :: fini
  real(8)                    :: press, Temp, of, R, Tout, pout, pin, Tin
  real(8), allocatable       :: sp_Y(:)
  real(8), allocatable       :: Y(:)
  real(8), allocatable       :: RT(:), AT(:)
  character(20), allocatable :: sp_name(:)
  integer                    :: err, s, Ns, neq, nstep, n, isc
  real(8)                    :: rtol, atol, timein, timeout, dt=0d0, tlim=0d0, time1, time2
  character(256)             :: wholestring
  character(20)              :: args(20), solver(20)
  integer                    :: uni, iopt(3)

  ! Load input.ini
  call fini%load(filename='input.ini')

  ! Read input.ini
  call obj_tab%ini_setup(fini=fini)
  call obj_chem%ini_setup(fini=fini)

  call fini%get(section_name='MOSE-Batch', option_name='p', val=pin, error=err)
  call fini%get(section_name='MOSE-Batch', option_name='T', val=Tin, error=err)
  call fini%get(section_name='MOSE-Batch', option_name='of', val=of, error=err)
  press = press*1d5

  call fini%get(section_name='MOSE-Batch', option_name='dt', val=dt, error=err)
  call fini%get(section_name='MOSE-Batch', option_name='to', val=tlim, error=err)
  call fini%get(section_name='MOSE-Batch', option_name='ns', val=nstep, error=err)

  neq = nsc + 1
  allocate(Y(neq))

  allocate(sp_name(nsc))
  allocate(sp_Y(nsc))

  open(newunit=uni,file='INPUT/phase.txt',status='old')
  read(uni,*)
  do s = 1, nsc
    read(uni,'(A)') wholestring
    call parse(wholestring,' ',args)
    sp_name(s) = args(1)
  end do
  close(uni)

  sp_Y = 1d-20
  do s = 1, nsc
    call fini%get(section_name='MOSE-Batch', option_name='y'//trim(sp_name(s)), val=sp_Y(s), error=err)
  enddo

  iopt = 0
  
  Ns = 0
  call fini%get(section_name='MOSE-ODE', option_name='solver', val=wholestring, error=err)
  call parse(wholestring,' ',args)
  do s = 1, size(args)
    if (trim(args(s))=='') exit
    Ns = Ns + 1
    solver(Ns) = trim(args(s))
  enddo

  allocate(RT(neq),AT(neq))

  call fini%get(section_name='MOSE-ODE', option_name='max-steps-ode', val=iopt(1), error=err)
  if (err/=0) iopt(1)=100000

  call fini%get(section_name='MOSE-ODE', option_name='relative-tol-species', val=rtol, error=err)
  if (err/=0) then
    RT(1:nsc)=1d-5
  else
    RT(1:nsc)=rtol
  end if
  
  call fini%get(section_name='MOSE-ODE', option_name='relative-tol-pressure', val=rtol, error=err)
  if (err/=0) then
    RT(neq)=1d-5
  else
    RT(neq)=rtol
  end if

  call fini%get(section_name='MOSE-ODE', option_name='absolute-tol-species', val=atol, error=err)
  if (err/=0) then
    AT(1:nsc)=1d-5
  else
    AT(1:nsc)=atol
  end if

  call fini%get(section_name='MOSE-ODE', option_name='absolute-tol-pressure', val=atol, error=err)
  if (err/=0) then
    AT(neq)=1d-5
  else
    AT(neq)=atol
  end if

  if (tlim/=0 .and. dt/=0) then
    nstep = tlim/dt
  elseif (tlim/=0 .and. nstep/=0) then
    dt = tlim/nstep
  elseif (dt/=0 .and. nstep/=0) then
    tlim = nstep*dt
  endif

!   ! ODE solvers

  call initialize
  ! write(*,*) ' rho = ', sum(Y(1:neq-1))
  ! write(*,*) ' T   = ',Tin
  ! write(*,*) ' p   = ',pin
  ! write(*,*)

  open(newunit=uni,file='results.txt')

  do s = 1, Ns
    write(*,*) ' Running ',trim(solver(s))

    call setup_odesolver(N=neq,solver=solver(s),RT=RT,AT=AT,iopt=iopt)
    call initialize

    call cpu_time(time1)
    do n = 1, nstep
      timein = timeout; timeout = timeout+dt
      call run_odesolver(neq,timein,timeout,Y,rhs,err)
      call co_Rtot(Y(1:neq-1),R)
      Tout = y(neq) / (sum(Y(1:neq-1))*R)
      pout = Tout*sum(Y(1:neq-1))*R
      write(*,*) timeout, pout, Tout
    enddo
    call cpu_time(time2)

    call co_Rtot(Y(1:neq-1),R)
    Tout = y(neq) / (sum(Y(1:neq-1))*R)
    pout = Tout*sum(Y(1:neq-1))*R

    write(uni,*) ' ODE solver: ', solver(s)
    write(uni,*) ' Execution time: ', (time2-time1)*1d3, 'msec'
    do isc = 1, nsc
      write(uni,'(A10,E12.4)') 'y'//trim(sp_name(isc)), Y(isc)/sum(Y(1:neq-1))
    end do
    write(uni,*)
    write(uni,*) 'p', pout
    write(uni,*) 'T', Tout

  enddo

  close(uni)
  ! write(*,*)
  ! write(*,*) ' Done. Check output in result.txt'

contains

  subroutine initialize()
    implicit none
    real(8) :: rhoi(neq-1), rho

    rhoi = sp_Y
    call co_Rtot(rhoi,R)
    rho = pin/(R*Tin)

    Y(1:nsc) = sp_Y*rho
    Y(neq) = pin!Tin
    
    timein  = 0.D0
    timeout = 0.D0
  end subroutine

end program batch