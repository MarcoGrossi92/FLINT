program test
  use oslo
  use U_cantera_load
  use U_Lib_Thermodynamic
  use U_IO_Table
  use U_Lib_Chemistry_data
  use U_Lib_Chemistry_wdot
  use U_Lib_Chemistry_rhs, only: rhs_cantera, rhs_native, gas
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

  solver = 'dvodef90'

  open(unit=10, file='comp-batch-general.dat', status='replace', form='formatted')
  open(unit=20, file='comp-batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(unit=30, file='comp-batch-cantera.dat', status='replace', form='formatted')
# endif

  !-------------------------------------------------------------------------------------------------
  ! WD
  !-------------------------------------------------------------------------------------------------

  call Read_IdealGas_Properties('WD')
  err = read_chemistry_file( folder='WD', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, 'WD/WD.yaml')
# endif

  open(200, file='WD-batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='WD-batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 0.01
  pin = 1.0d+5
  Tin = 1000
  nstep = 1000
  dt = tlim/nstep

  neq = nsc + 1
  allocate(Y(neq))
  allocate(sp_Y(nsc))

  sp_Y = 1d-20
  sp_Y(1) = 0.2
  sp_Y(2) = 0.8

  iopt = 0
  allocate(RT(neq),AT(neq))
  iopt(1)=100000
  RT(1:nsc)=1d-7
  RT(neq)=1d-7
  AT(1:nsc)=1d-7
  AT(neq)=1d-7
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  !! Cantera
# if defined (CANTERA)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_cantera,err)
    Tout = y(neq)
    write(300,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'WD Cantera time =', time2-time1
  write(30,*) nsc, time2-time1
# endif

  !! Native with coded mechanism
  call Assign_Mechanism(mech_name)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,err)
    Tout = y(neq)
    write(200,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'WD explicit time =', time2-time1
  write(20,*) nsc, time2-time1

  close(200)
# if defined (CANTERA)
  close(300)
# endif

  deallocate(Y); deallocate(sp_Y)
  deallocate(AT); deallocate(RT)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(mi_tab); deallocate(k_tab)
  deallocate(kf_tab); deallocate(kb_tab); deallocate(ni1_tab); deallocate(ni2_tab); deallocate(epsch_tab)

  !-------------------------------------------------------------------------------------------------
  ! Troyes
  !-------------------------------------------------------------------------------------------------

  call Read_IdealGas_Properties('Troyes')
  err = read_chemistry_file( folder='Troyes', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, 'Troyes/troyes.yaml')
# endif

  open(100, file='Troyes-batch-general.dat', status='replace', form='formatted')
  open(200, file='Troyes-batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='Troyes-batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 0.02
  pin = 1.0d+5
  Tin = 1000
  nstep = 1000
  dt = tlim/nstep

  neq = nsc + 1
  allocate(Y(neq))
  allocate(sp_Y(nsc))

  sp_Y = 1d-20
  sp_Y(6) = 0.00534
  sp_Y(10) = 0.18796
  sp_Y(12) = 0.80670

  iopt = 0
  allocate(RT(neq),AT(neq))
  iopt(1)=100000
  RT(1:nsc)=1d-7
  RT(neq)=1d-7
  AT(1:nsc)=1d-7
  AT(neq)=1d-7
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  !! Cantera
# if defined (CANTERA)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_cantera,err)
    Tout = y(neq)
    write(300,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'Troyes Cantera time =', time2-time1
  write(30,*) nsc, time2-time1
# endif

  !! Native with coded mechanism
  call Assign_Mechanism(mech_name)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,err)
    Tout = y(neq)
    write(100,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'Troyes explicit time =', time2-time1
  write(20,*) nsc, time2-time1

  !! Native without coded mechanism
  call Assign_Mechanism('nemo')
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,err)
    Tout = y(neq)
    write(200,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'Troyes general time =', time2-time1
  write(10,*) nsc, time2-time1

  close(100); close(200)
# if defined (CANTERA)
  close(300)
# endif

  deallocate(Y); deallocate(sp_Y)
  deallocate(AT); deallocate(RT)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(mi_tab); deallocate(k_tab)
  deallocate(kf_tab); deallocate(kb_tab); deallocate(ni1_tab); deallocate(ni2_tab); deallocate(epsch_tab)

  !-------------------------------------------------------------------------------------------------
  ! Ecker
  !-------------------------------------------------------------------------------------------------

  call Read_IdealGas_Properties('Ecker')
  err = read_chemistry_file( folder='Ecker', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, 'Ecker/ecker.yaml')
# endif

  open(100, file='Ecker-batch-general.dat', status='replace', form='formatted')
  open(200, file='Ecker-batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='Ecker-batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 0.02
  pin = 1.0d+5
  Tin = 1000
  nstep = 1000
  dt = tlim/nstep

  neq = nsc + 1
  allocate(Y(neq))
  allocate(sp_Y(nsc))

  sp_Y = 1d-20
  sp_Y(12) = 0.18798856d0
  sp_Y(1) = 0.00534534d0
  sp_Y(14) = 0.8066661d0

  iopt = 0
  allocate(RT(neq),AT(neq))
  iopt(1)=100000
  RT(1:nsc)=1d-7
  RT(neq)=1d-7
  AT(1:nsc)=1d-7
  AT(neq)=1d-7
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  !! Cantera
# if defined (CANTERA)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_cantera,err)
    Tout = y(neq)
    write(300,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'Ecker Cantera time =', time2-time1
  write(30,*) nsc, time2-time1
# endif

  !! Native with coded mechanism
  call Assign_Mechanism(mech_name)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,err)
    Tout = y(neq)
    write(100,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'Ecker explicit time =', time2-time1
  write(20,*) nsc, time2-time1

  !! Native without coded mechanism
  call Assign_Mechanism('nemo')
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,err)
    Tout = y(neq)
    write(200,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'Ecker general time =', time2-time1
  write(10,*) nsc, time2-time1

  close(100); close(200)
# if defined (CANTERA)
  close(300)
# endif

  deallocate(Y); deallocate(sp_Y)
  deallocate(AT); deallocate(RT)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(mi_tab); deallocate(k_tab)
  deallocate(kf_tab); deallocate(kb_tab); deallocate(ni1_tab); deallocate(ni2_tab); deallocate(epsch_tab)

  !-------------------------------------------------------------------------------------------------
  ! Cross
  !-------------------------------------------------------------------------------------------------

  call Read_IdealGas_Properties('Cross')
  err = read_chemistry_file( folder='Cross', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, 'Cross/cross.yaml')
# endif

  open(100, file='Cross-batch-general.dat', status='replace', form='formatted')
  open(200, file='Cross-batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='Cross-batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 0.02
  pin = 1.0d+5
  Tin = 1005
  nstep = 1000
  dt = tlim/nstep

  neq = nsc + 1
  allocate(Y(neq))
  allocate(sp_Y(nsc))

  sp_Y = 1d-20
  sp_Y(9) = 0.18798856d0
  sp_Y(12) = 0.00534534d0
  sp_Y(17) = 0.8066661d0

  iopt = 0
  allocate(RT(neq),AT(neq))
  iopt(1)=100000
  RT(1:nsc)=1d-7
  RT(neq)=1d-7
  AT(1:nsc)=1d-7
  AT(neq)=1d-7
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  !! Cantera
# if defined (CANTERA)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_cantera,err)
    Tout = y(neq)
    write(300,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'Cross Cantera time =', time2-time1
  write(30,*) nsc, time2-time1
# endif

  !! Native with coded mechanism
  call Assign_Mechanism(mech_name)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,err)
    Tout = y(neq)
    write(100,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'Cross explicit time =', time2-time1
  write(20,*) nsc, time2-time1

  !! Native without coded mechanism
  call Assign_Mechanism('nemo')
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,err)
    Tout = y(neq)
    write(200,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'Cross general time =', time2-time1
  write(10,*) nsc, time2-time1

  close(100); close(200)
# if defined (CANTERA)
  close(300)
# endif

  deallocate(Y); deallocate(sp_Y)
  deallocate(AT); deallocate(RT)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(mi_tab); deallocate(k_tab)
  deallocate(kf_tab); deallocate(kb_tab); deallocate(ni1_tab); deallocate(ni2_tab); deallocate(epsch_tab)

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

end program test