program test
  use oslo
  use FLINT_cantera_load
  use FLINT_Lib_Thermodynamic
  use FLINT_IO_ThermoTransport
  use FLINT_Lib_Chemistry_data
  use FLINT_Lib_Chemistry_wdot
  use FLINT_Lib_Chemistry_rhs
  use FLINT_IO_chemistry
  implicit none
  real(8)                    :: R, Tout, pin, Tin, rho
  real(8), allocatable       :: sp_Y(:)
  real(8), allocatable       :: Y(:)
  real(8), allocatable       :: RT(:), AT(:)
  integer                    :: err, neq, nstep, n
  real(8)                    :: timein, timeout, dt=0d0, tlim=0d0, time1, time2
  character(32)              :: solver, mech_name
  integer                    :: iopt(3)

  solver = 'cvode'
  nstep = 1000
  iopt = 0
  iopt(1) = 1000000

  open(unit=10, file='comp-batch-general.dat', status='replace', form='formatted')
  open(unit=20, file='comp-batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(unit=30, file='comp-batch-canteraFor.dat', status='replace', form='formatted')
# endif

  !-------------------------------------------------------------------------------------------------
  ! WD
  !-------------------------------------------------------------------------------------------------

  err = read_idealgas_properties('WD/INPUT')
  err = read_chemistry_file( folder='WD/INPUT', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, 'WD/INPUT/WD.yaml')
# endif

  open(200, file='WD/OUTPUT/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='WD/OUTPUT/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 8.d-3
  pin = 1.0d+5
  Tin = 1000
  dt = tlim/nstep

  neq = nsc + 1
  allocate(Y(neq))
  allocate(sp_Y(nsc))

  sp_Y = 1d-20
  sp_Y(1) = 0.2
  sp_Y(2) = 0.8

  allocate(RT(neq),AT(neq))
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
  write(30,*) 'WD', time2-time1
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
  write(20,*) 'WD', time2-time1

  close(200)
# if defined (CANTERA)
  close(300)
# endif

  deallocate(Y); deallocate(sp_Y)
  deallocate(AT); deallocate(RT)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! Troyes
  !-------------------------------------------------------------------------------------------------

  err = read_idealgas_properties('Troyes/INPUT')
  err = read_chemistry_file( folder='Troyes/INPUT', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, 'Troyes/INPUT/troyes.yaml')
# endif

  open(100, file='Troyes/OUTPUT/batch-general.dat', status='replace', form='formatted')
  open(200, file='Troyes/OUTPUT/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='Troyes/OUTPUT/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 5d-3
  pin = 1.0d+5
  Tin = 1000
  dt = tlim/nstep

  neq = nsc + 1
  allocate(Y(neq))
  allocate(sp_Y(nsc))

  sp_Y = 1d-20
  sp_Y(6) = 0.00534
  sp_Y(10) = 0.18796
  sp_Y(12) = 0.80670

  allocate(RT(neq),AT(neq))
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
  write(30,*) 'Troyes', time2-time1
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

  write(*,*) 'Troyes explicit time =', time2-time1
  write(20,*) 'Troyes', time2-time1

  !! Native without coded mechanism
  call Assign_Mechanism('nemo')
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,err)
    Tout = y(neq)
    write(100,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'Troyes general time =', time2-time1
  write(10,*) 'Troyes', time2-time1

  close(100); close(200)
# if defined (CANTERA)
  close(300)
# endif

  deallocate(Y); deallocate(sp_Y)
  deallocate(AT); deallocate(RT)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! Ecker
  !-------------------------------------------------------------------------------------------------

  err = read_idealgas_properties('Ecker/INPUT')
  err = read_chemistry_file( folder='Ecker/INPUT', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, 'Ecker/INPUT/ecker.yaml')
# endif

  open(100, file='Ecker/OUTPUT/batch-general.dat', status='replace', form='formatted')
  open(200, file='Ecker/OUTPUT/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='Ecker/OUTPUT/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 2d-2
  pin = 1.0d+5
  Tin = 1000d0
  dt = tlim/nstep

  neq = nsc + 1
  allocate(Y(neq))
  allocate(sp_Y(nsc))

  sp_Y = 1d-20
  sp_Y(12) = 0.18798856d0
  sp_Y(1) = 0.00534534d0
  sp_Y(14) = 0.8066661d0

  allocate(RT(neq),AT(neq))
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
  write(30,*) 'Ecker', time2-time1
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

  write(*,*) 'Ecker explicit time =', time2-time1
  write(20,*) 'Ecker', time2-time1

  !! Native without coded mechanism
  call Assign_Mechanism('nemo')
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,err)
    Tout = y(neq)
    write(100,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'Ecker general time =', time2-time1
  write(10,*) 'Ecker', time2-time1

  close(100); close(200)
# if defined (CANTERA)
  close(300)
# endif

  deallocate(Y); deallocate(sp_Y)
  deallocate(AT); deallocate(RT)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! Cross
  !-------------------------------------------------------------------------------------------------

  err = read_idealgas_properties('Cross/INPUT')
  err = read_chemistry_file( folder='Cross/INPUT', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, 'Cross/INPUT/cross.yaml')
# endif

  open(100, file='Cross/OUTPUT/batch-general.dat', status='replace', form='formatted')
  open(200, file='Cross/OUTPUT/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='Cross/OUTPUT/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 1d-2
  pin = 1.0d+5
  Tin = 1010
  dt = tlim/nstep

  neq = nsc + 1
  allocate(Y(neq))
  allocate(sp_Y(nsc))

  sp_Y = 1d-20
  sp_Y(9) = 0.18798856d0
  sp_Y(12) = 0.00534534d0
  sp_Y(17) = 0.8066661d0

  allocate(RT(neq),AT(neq))
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
  write(30,*) 'Cross', time2-time1
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

  write(*,*) 'Cross explicit time =', time2-time1
  write(20,*) 'Cross', time2-time1

  !! Native without coded mechanism
  call Assign_Mechanism('nemo')
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,err)
    Tout = y(neq)
    write(100,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'Cross general time =', time2-time1
  write(10,*) 'Cross', time2-time1

  close(100); close(200)
# if defined (CANTERA)
  close(300)
# endif

  deallocate(Y); deallocate(sp_Y)
  deallocate(AT); deallocate(RT)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! Smooke
  !-------------------------------------------------------------------------------------------------

  err = read_idealgas_properties('Smooke/INPUT')
  err = read_chemistry_file( folder='Smooke/INPUT', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, 'Smooke/INPUT/smooke.yaml')
# endif

  open(100, file='Smooke/OUTPUT/batch-general.dat', status='replace', form='formatted')
  open(200, file='Smooke/OUTPUT/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='Smooke/OUTPUT/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 0.2d0
  pin = 1.0d+5
  Tin = 1300d0
  dt = tlim/nstep

  neq = nsc + 1
  allocate(Y(neq))
  allocate(sp_Y(nsc))

  sp_Y = 1d-20
  sp_Y(1) = 0.0552d0
  sp_Y(3) = 0.2201d0
  sp_Y(16) = 0.7247d0

  allocate(RT(neq),AT(neq))
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

  write(*,*) 'Smooke Cantera time =', time2-time1
  write(30,*) 'Smooke', time2-time1
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

  write(*,*) 'Smooke explicit time =', time2-time1
  write(20,*) 'Smooke', time2-time1

  !! Native without coded mechanism
  call Assign_Mechanism('nemo')
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,err)
    Tout = y(neq)
    write(100,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'Smooke general time =', time2-time1
  write(10,*) 'Smooke', time2-time1

  close(100); close(200)
# if defined (CANTERA)
  close(300)
# endif

  deallocate(Y); deallocate(sp_Y)
  deallocate(AT); deallocate(RT)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! CORIA
  !-------------------------------------------------------------------------------------------------

  err = read_idealgas_properties('CORIA/INPUT')
  err = read_chemistry_file( folder='CORIA/INPUT', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, 'CORIA/INPUT/coria.yaml')
# endif

  open(100, file='CORIA/OUTPUT/batch-general.dat', status='replace', form='formatted')
  open(200, file='CORIA/OUTPUT/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='CORIA/OUTPUT/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 0.005
  pin = 1.0d+5
  Tin = 1300d0
  dt = tlim/nstep

  neq = nsc + 1
  allocate(Y(neq))
  allocate(sp_Y(nsc))

  sp_Y = 1d-20
  sp_Y(9) = 0.2d0
  sp_Y(4) = 0.8d0

  allocate(RT(neq),AT(neq))
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

  write(*,*) 'CORIA Cantera time =', time2-time1
  write(30,*) 'CORIA', time2-time1
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

  write(*,*) 'CORIA explicit time =', time2-time1
  write(20,*) 'CORIA', time2-time1

  !! Native without coded mechanism
  call Assign_Mechanism('nemo')
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,err)
    Tout = y(neq)
    write(100,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'CORIA general time =', time2-time1
  write(10,*) 'CORIA', time2-time1

  close(100); close(200)
# if defined (CANTERA)
  close(300)
# endif

  deallocate(Y); deallocate(sp_Y)
  deallocate(AT); deallocate(RT)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! TSR-CDF-13
  !-------------------------------------------------------------------------------------------------

  err = read_idealgas_properties('TSR-CDF-13/INPUT')
  err = read_chemistry_file( folder='TSR-CDF-13/INPUT', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, 'TSR-CDF-13/INPUT/TSR-CDF-13.yaml')
# endif

  open(100, file='TSR-CDF-13/OUTPUT/batch-general.dat', status='replace', form='formatted')
  open(200, file='TSR-CDF-13/OUTPUT/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='TSR-CDF-13/OUTPUT/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 0.005
  pin = 5.0d+5
  Tin = 1300d0
  dt = tlim/nstep

  neq = nsc + 1
  allocate(Y(neq))
  allocate(sp_Y(nsc))

  sp_Y = 1d-20
  sp_Y(7) = 0.2d0
  sp_Y(10) = 0.8d0

  allocate(RT(neq),AT(neq))
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

  write(*,*) 'TSR-CDF-13 Cantera time =', time2-time1
  write(30,*) 'TSR-CDF-13', time2-time1
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

  write(*,*) 'TSR-CDF-13 explicit time =', time2-time1
  write(20,*) 'TSR-CDF-13', time2-time1

  !! Native without coded mechanism
  call Assign_Mechanism('nemo')
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,err)
    Tout = y(neq)
    write(100,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'TSR-CDF-13 general time =', time2-time1
  write(10,*) 'TSR-CDF-13', time2-time1

  close(100); close(200)
# if defined (CANTERA)
  close(300)
# endif

  deallocate(Y); deallocate(sp_Y)
  deallocate(AT); deallocate(RT)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! Pelucchi
  !-------------------------------------------------------------------------------------------------

  err = read_idealgas_properties('Pelucchi/INPUT')
  !err = read_chemistry_file( folder='Pelucchi/INPUT', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, 'Pelucchi/INPUT/pelucchi.yaml')
# endif

  open(100, file='Pelucchi/OUTPUT/batch-general.dat', status='replace', form='formatted')
  open(200, file='Pelucchi/OUTPUT/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='Pelucchi/OUTPUT/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 5d-2
  pin = 1.0d+5
  Tin = 1250d0
  dt = tlim/nstep

  neq = nsc + 1
  allocate(Y(neq))
  allocate(sp_Y(nsc))

  sp_Y = 1d-20
  sp_Y(18) = 0.00859d0
  sp_Y(14) = 0.00606d0
  sp_Y(16) = 0.00365d0
  sp_Y(1) = 0.00025d0
  sp_Y(13) = 0.98044d0

  allocate(RT(neq),AT(neq))
  RT(1:nsc)=1d-12
  RT(neq)=1d-12
  AT(1:nsc)=1d-15
  AT(neq)=1d-15
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

  write(*,*) 'Pelucchi Cantera time =', time2-time1
  write(30,*) 'Pelucchi', time2-time1
# endif

  ! !! Native with coded mechanism
  ! call Assign_Mechanism(mech_name)
  ! call initialize
  ! call cpu_time(time1)
  ! do n = 1, nstep
  !   timein = timeout; timeout = timeout+dt
  !   call run_odesolver(neq,timein,timeout,Y,rhs_native,err)
  !   Tout = y(neq)
  !   write(200,*) timeout, Tout
  ! enddo
  ! call cpu_time(time2)

  ! write(*,*) 'Pelucchi explicit time =', time2-time1
  ! write(20,*) 'Pelucchi', time2-time1

  ! !! Native without coded mechanism
  ! call Assign_Mechanism('nemo')
  ! call initialize
  ! call cpu_time(time1)
  ! do n = 1, nstep
  !   timein = timeout; timeout = timeout+dt
  !   call run_odesolver(neq,timein,timeout,Y,rhs_native,err)
  !   Tout = y(neq)
  !   write(100,*) timeout, Tout
  ! enddo
  ! call cpu_time(time2)

  write(*,*) 'Pelucchi general time =', time2-time1
  write(10,*) 'Pelucchi', time2-time1

  close(100); close(200)
# if defined (CANTERA)
  close(300)
# endif

  deallocate(Y); deallocate(sp_Y)
  deallocate(AT); deallocate(RT)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  call free_chemistry_data()

contains

  subroutine initialize()
    implicit none
    R = Rtot(sp_Y)
    rho = pin/(R*Tin)
    Y(1:nsc) = rho*sp_Y
    Y(neq) = Tin
    timein  = 0.D0
    timeout = 0.D0
  end subroutine

end program test