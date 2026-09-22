program test
  use oslo
# if defined (CANTERA)
  use FLINT_cantera_load
# endif
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  use FLINT_Lib_Chemistry_data
  use FLINT_Lib_Chemistry_wdot
  use FLINT_Lib_Chemistry_rhs
  use FLINT_Load_chemistry
  implicit none
  real(8)                    :: R, Tout, pin, Tin, rho
  real(8), allocatable       :: sp_Y(:)
  real(8), allocatable       :: Y(:)
  real(8), allocatable       :: RT(:), AT(:)
  integer                    :: err, neq, nstep, n
  real(8)                    :: timein, timeout, dt=0d0, tlim=0d0, time1, time2
  character(32)              :: solver, mech_name
  integer                    :: iopt(3), sim_type

# if defined(SUNDIALS)
  solver = 'cvode'
# else
  solver = 'ros4'
# endif
  write(*,*)'what kind of simulation do you want to run?'
  write(*,*)'1) verification'
  write(*,*)'2) performance'
  read(*,*) sim_type
  if (sim_type==1) then
    nstep = 1000
  elseif (sim_type==2) then
    nstep = 1
  else
    stop "choose 1 or 2!"
  endif
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

  call execute_command_line('mkdir -p WD/')
  err = read_idealgas_thermo('../database/WD/')
  err = read_chemistry( folder='../database/WD/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../database/WD/WD.yaml')
# endif

  open(200, file='WD/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='WD/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 8.d-3
  pin = 1.0d+5
  Tin = 1000
  dt = tlim/nstep

  neq = ns + 1
  allocate(Y(neq))
  allocate(sp_Y(ns))

  sp_Y = 1d-20
  sp_Y(1) = 0.2
  sp_Y(2) = 0.8

  allocate(RT(neq),AT(neq))
  RT(1:ns)=1d-7
  RT(neq)=1d-7
  AT(1:ns)=1d-7
  AT(neq)=1d-7
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  !! Cantera
# if defined (CANTERA)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_cantera,no_jacobian,0,err)
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
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
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
  deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! Troyes
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p Troyes/')
  err = read_idealgas_thermo('../database/Troyes/')
  err = read_chemistry( folder='../database/Troyes/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../database/Troyes/troyes.yaml')
# endif

  open(100, file='Troyes/batch-general.dat', status='replace', form='formatted')
  open(200, file='Troyes/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='Troyes/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 5d-3
  pin = 1.0d+5
  Tin = 1000
  dt = tlim/nstep

  neq = ns + 1
  allocate(Y(neq))
  allocate(sp_Y(ns))

  sp_Y = 1d-20
  sp_Y(6) = 0.00534
  sp_Y(10) = 0.18796
  sp_Y(12) = 0.80670

  allocate(RT(neq),AT(neq))
  RT(1:ns)=1d-7
  RT(neq)=1d-7
  AT(1:ns)=1d-7
  AT(neq)=1d-7
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  !! Cantera
# if defined (CANTERA)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_cantera,no_jacobian,0,err)
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
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
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
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
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
  deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! Ecker
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p Ecker/')
  err = read_idealgas_thermo('../database/Ecker/')
  err = read_chemistry( folder='../database/Ecker/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../database/Ecker/ecker.yaml')
# endif

  open(100, file='Ecker/batch-general.dat', status='replace', form='formatted')
  open(200, file='Ecker/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='Ecker/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 5d-3
  pin = 1.0d+5
  Tin = 1000d0
  dt = tlim/nstep

  neq = ns + 1
  allocate(Y(neq))
  allocate(sp_Y(ns))

  sp_Y = 1d-20
  sp_Y(12) = 0.18798856d0
  sp_Y(1) = 0.00534534d0
  sp_Y(14) = 0.8066661d0

  allocate(RT(neq),AT(neq))
  RT(1:ns)=1d-7
  RT(neq)=1d-7
  AT(1:ns)=1d-7
  AT(neq)=1d-7
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  !! Cantera
# if defined (CANTERA)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_cantera,no_jacobian,0,err)
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
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
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
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
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
  deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! Cross
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p Cross/')
  err = read_idealgas_thermo('../database/Cross/')
  err = read_chemistry( folder='../database/Cross/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../database/Cross/cross.yaml')
# endif

  open(100, file='Cross/batch-general.dat', status='replace', form='formatted')
  open(200, file='Cross/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='Cross/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 1d-2
  pin = 1.0d+5
  Tin = 1010
  dt = tlim/nstep

  neq = ns + 1
  allocate(Y(neq))
  allocate(sp_Y(ns))

  sp_Y = 1d-20
  sp_Y(9) = 0.18798856d0
  sp_Y(12) = 0.00534534d0
  sp_Y(17) = 0.8066661d0

  allocate(RT(neq),AT(neq))
  RT(1:ns)=1d-7
  RT(neq)=1d-7
  AT(1:ns)=1d-7
  AT(neq)=1d-7
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  !! Cantera
# if defined (CANTERA)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_cantera,no_jacobian,0,err)
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
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
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
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
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
  deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! Smooke
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p Smooke/')
  err = read_idealgas_thermo('../database/Smooke/')
  err = read_chemistry( folder='../database/Smooke/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../database/Smooke/smooke.yaml')
# endif

  open(100, file='Smooke/batch-general.dat', status='replace', form='formatted')
  open(200, file='Smooke/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='Smooke/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 0.2d0
  pin = 1.0d+5
  Tin = 1300d0
  dt = tlim/nstep

  neq = ns + 1
  allocate(Y(neq))
  allocate(sp_Y(ns))

  sp_Y = 1d-20
  sp_Y(1) = 0.0552d0
  sp_Y(3) = 0.2201d0
  sp_Y(16) = 0.7247d0

  allocate(RT(neq),AT(neq))
  RT(1:ns)=1d-7
  RT(neq)=1d-7
  AT(1:ns)=1d-7
  AT(neq)=1d-7
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  !! Cantera
# if defined (CANTERA)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_cantera,no_jacobian,0,err)
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
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
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
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
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
  deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! CORIA-CNRS
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p CORIA/')
  err = read_idealgas_thermo('../database/CORIA/')
  err = read_chemistry( folder='../database/CORIA/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../database/CORIA/coria.yaml')
# endif

  open(100, file='CORIA/batch-general.dat', status='replace', form='formatted')
  open(200, file='CORIA/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='CORIA/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 0.005
  pin = 1.0d+5
  Tin = 1300d0
  dt = tlim/nstep

  neq = ns + 1
  allocate(Y(neq))
  allocate(sp_Y(ns))

  sp_Y = 1d-20
  sp_Y(9) = 0.2d0
  sp_Y(4) = 0.8d0

  allocate(RT(neq),AT(neq))
  RT(1:ns)=1d-7
  RT(neq)=1d-7
  AT(1:ns)=1d-7
  AT(neq)=1d-7
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  !! Cantera
# if defined (CANTERA)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_cantera,no_jacobian,0,err)
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
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
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
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
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
  deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! TSR-CDF-13
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p TSR-CDF-13')
  err = read_idealgas_thermo('../database/TSR-CDF-13/')
  err = read_chemistry( folder='../database/TSR-CDF-13/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../database/TSR-CDF-13/TSR-CDF-13.yaml')
# endif

  open(100, file='TSR-CDF-13/batch-general.dat', status='replace', form='formatted')
  open(200, file='TSR-CDF-13/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='TSR-CDF-13/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 0.005
  pin = 5.0d+5
  Tin = 1300d0
  dt = tlim/nstep

  neq = ns + 1
  allocate(Y(neq))
  allocate(sp_Y(ns))

  sp_Y = 1d-20
  sp_Y(7) = 0.2d0
  sp_Y(10) = 0.8d0

  allocate(RT(neq),AT(neq))
  RT(1:ns)=1d-7
  RT(neq)=1d-7
  AT(1:ns)=1d-7
  AT(neq)=1d-7
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  !! Cantera
# if defined (CANTERA)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_cantera,no_jacobian,0,err)
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
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
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
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
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
  deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! Pelucchi
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p Pelucchi')
  err = read_idealgas_thermo('../database/Pelucchi/')
  err = read_chemistry( folder='../database/Pelucchi/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../database/Pelucchi/pelucchi.yaml')
# endif

  open(100, file='Pelucchi/batch-general.dat', status='replace', form='formatted')
  open(200, file='Pelucchi/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='Pelucchi/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 5d-2
  pin = 1.0d+5
  Tin = 1250d0
  dt = tlim/nstep

  neq = ns + 1
  allocate(Y(neq))
  allocate(sp_Y(ns))

  sp_Y = 1d-20
  sp_Y(18) = 0.00859d0
  sp_Y(14) = 0.00606d0
  sp_Y(16) = 0.00365d0
  sp_Y(1) = 0.00025d0
  sp_Y(13) = 0.98044d0

  allocate(RT(neq),AT(neq))
  RT(1:ns)=1d-12
  RT(neq)=1d-12
  AT(1:ns)=1d-15
  AT(neq)=1d-15
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  !! Cantera
# if defined (CANTERA)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_cantera,no_jacobian,0,err)
    Tout = y(neq)
    write(300,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'Pelucchi Cantera time =', time2-time1
  write(30,*) 'Pelucchi', time2-time1
# endif

  !! Native with coded mechanism
  call Assign_Mechanism(mech_name)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
    Tout = y(neq)
    write(200,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'Pelucchi explicit time =', time2-time1
  write(20,*) 'Pelucchi', time2-time1

  !! Native without coded mechanism
  call Assign_Mechanism('nemo')
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
    Tout = y(neq)
    write(100,*) timeout, Tout
  enddo
  call cpu_time(time2)

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
  deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! ZK
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p ZK/')
  err = read_idealgas_thermo('../database/ZK/')
  err = read_chemistry( folder='../database/ZK/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../database/ZK/ZK.yaml')
# endif

  open(100, file='ZK/batch-general.dat', status='replace', form='formatted')
  open(200, file='ZK/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='ZK/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 0.002
  pin = 5.0d+5
  Tin = 1300d0
  dt = tlim/nstep

  neq = ns + 1
  allocate(Y(neq))
  allocate(sp_Y(ns))

  sp_Y = 1d-20
  sp_Y(17) = 0.2d0
  sp_Y(6) = 0.8d0

  allocate(RT(neq),AT(neq))
  RT(1:ns)=1d-7
  RT(neq)=1d-7
  AT(1:ns)=1d-7
  AT(neq)=1d-7
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  !! Cantera
# if defined (CANTERA)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_cantera,no_jacobian,0,err)
    Tout = y(neq)
    write(300,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'ZK Cantera time =', time2-time1
  write(30,*) 'ZK', time2-time1
# endif

  !! Native with coded mechanism
  call Assign_Mechanism(mech_name)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
    Tout = y(neq)
    write(200,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'ZK explicit time =', time2-time1
  write(20,*) 'ZK', time2-time1

  !! Native without coded mechanism
  call Assign_Mechanism('nemo')
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
    Tout = y(neq)
    write(100,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'ZK general time =', time2-time1
  write(10,*) 'ZK', time2-time1

  close(100); close(200)
# if defined (CANTERA)
  close(300)
# endif

  deallocate(Y); deallocate(sp_Y)
  deallocate(AT); deallocate(RT)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! TSR-GP-24
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p TSR-GP-24/')
  err = read_idealgas_thermo('../database/TSR-GP-24/')
  err = read_chemistry( folder='../database/TSR-GP-24/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../database/TSR-GP-24/TSR-GP-24.yaml')
# endif

  open(100, file='TSR-GP-24/batch-general.dat', status='replace', form='formatted')
  open(200, file='TSR-GP-24/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='TSR-GP-24/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 0.002
  pin = 5.0d+5
  Tin = 1300d0
  dt = tlim/nstep

  neq = ns + 1
  allocate(Y(neq))
  allocate(sp_Y(ns))

  sp_Y = 1d-20
  sp_Y(22) = 0.2d0
  sp_Y(19) = 0.8d0

  allocate(RT(neq),AT(neq))
  RT(1:ns)=1d-7
  RT(neq)=1d-7
  AT(1:ns)=1d-7
  AT(neq)=1d-7
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  !! Cantera
# if defined (CANTERA)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_cantera,no_jacobian,0,err)
    Tout = y(neq)
    write(300,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'TSR-GP-24 Cantera time =', time2-time1
  write(30,*) 'TSR-GP-24', time2-time1
# endif

  !! Native with coded mechanism
  call Assign_Mechanism(mech_name)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
    Tout = y(neq)
    write(200,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'TSR-GP-24 explicit time =', time2-time1
  write(20,*) 'TSR-GP-24', time2-time1

  !! Native without coded mechanism
  call Assign_Mechanism('nemo')
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
    Tout = y(neq)
    write(100,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'TSR-GP-24 general time =', time2-time1
  write(10,*) 'TSR-GP-24', time2-time1

  close(100); close(200)
# if defined (CANTERA)
  close(300)
# endif

  deallocate(Y); deallocate(sp_Y)
  deallocate(AT); deallocate(RT)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! TSR-Rich-31
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p TSR-Rich-31/')
  err = read_idealgas_thermo('../database/TSR-Rich-31/')
  err = read_chemistry( folder='../database/TSR-Rich-31/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../database/TSR-Rich-31/TSR-Rich-31.yaml')
# endif

  open(100, file='TSR-Rich-31/batch-general.dat', status='replace', form='formatted')
  open(200, file='TSR-Rich-31/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='TSR-Rich-31/batch-cantera.dat', status='replace', form='formatted')
# endif

  tlim = 0.005
  pin = 5.0d+5
  Tin = 1300d0
  dt = tlim/nstep

  neq = ns + 1
  allocate(Y(neq))
  allocate(sp_Y(ns))

  sp_Y = 1d-20
  sp_Y(4) = 0.2d0
  sp_Y(16) = 0.8d0

  allocate(RT(neq),AT(neq))
  RT(1:ns)=1d-7
  RT(neq)=1d-7
  AT(1:ns)=1d-7
  AT(neq)=1d-7
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  !! Cantera
# if defined (CANTERA)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_cantera,no_jacobian,0,err)
    Tout = y(neq)
    write(300,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'TSR-Rich-31 Cantera time =', time2-time1
  write(30,*) 'TSR-Rich-31', time2-time1
# endif

  !! Native with coded mechanism
  call Assign_Mechanism(mech_name)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
    Tout = y(neq)
    write(200,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'TSR-Rich-31 explicit time =', time2-time1
  write(20,*) 'TSR-Rich-31', time2-time1

  !! Native without coded mechanism
  call Assign_Mechanism('nemo')
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
    Tout = y(neq)
    write(100,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'TSR-Rich-31 general time =', time2-time1
  write(10,*) 'TSR-Rich-31', time2-time1

  close(100); close(200)
# if defined (CANTERA)
  close(300)
# endif

  deallocate(Y); deallocate(sp_Y)
  deallocate(AT); deallocate(RT)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! Gerlinger
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p Gerlinger/')
  err = read_idealgas_thermo('../database/Gerlinger/')
  err = read_chemistry( folder='../database/Gerlinger/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../database/Gerlinger/Gerlinger-9.yaml')
# endif

  open(100, file='Gerlinger/batch-general.dat', status='replace', form='formatted')
  open(200, file='Gerlinger/batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
  open(300, file='Gerlinger/batch-cantera.dat', status='replace', form='formatted')
# endif

  ! p = 1 bar, T = 1200 K, stoichiometric H2/air (2 H2 + O2 + 3.76 N2)
  tlim = 2.0d-4
  pin = 1.0d+5
  Tin = 1200d0
  dt = tlim/nstep

  neq = ns + 1
  allocate(Y(neq))
  allocate(sp_Y(ns))

  sp_Y = 1d-20
  sp_Y(1) = 0.745124d0   ! N2
  sp_Y(2) = 0.226354d0   ! O2
  sp_Y(3) = 0.028522d0   ! H2

  allocate(RT(neq),AT(neq))
  RT(1:ns)=1d-7
  RT(neq)=1d-7
  AT(1:ns)=1d-7
  AT(neq)=1d-7
  call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)

  !! Cantera
# if defined (CANTERA)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_cantera,no_jacobian,0,err)
    Tout = y(neq)
    write(300,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'Gerlinger Cantera time =', time2-time1
  write(30,*) 'Gerlinger', time2-time1
# endif

  !! Native with coded mechanism
  call Assign_Mechanism(mech_name)
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
    Tout = y(neq)
    write(200,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'Gerlinger explicit time =', time2-time1
  write(20,*) 'Gerlinger', time2-time1

  !! Native without coded mechanism
  call Assign_Mechanism('nemo')
  call initialize
  call cpu_time(time1)
  do n = 1, nstep
    timein = timeout; timeout = timeout+dt
    call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
    Tout = y(neq)
    write(100,*) timeout, Tout
  enddo
  call cpu_time(time2)

  write(*,*) 'Gerlinger general time =', time2-time1
  write(10,*) 'Gerlinger', time2-time1

  close(100); close(200)
# if defined (CANTERA)
  close(300)
# endif

  deallocate(Y); deallocate(sp_Y)
  deallocate(AT); deallocate(RT)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

contains

  subroutine initialize()
    implicit none
    R = f_Rtot(sp_Y)
    rho = pin/(R*Tin)
    Y(1:ns) = rho*sp_Y
    Y(neq) = Tin
    timein  = 0.D0
    timeout = 0.D0
  end subroutine

end program test