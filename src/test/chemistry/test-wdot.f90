! Test for wdot with thrid-body reactions
program test
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  use FLINT_Load_chemistry
  use FLINT_Lib_Chemistry_data
  use FLINT_Lib_Chemistry_wdot
  use FLINT_Lib_Chemistry_data, only: T_tab_min, T_tab_max
# if defined (CANTERA)
  use FLINT_Lib_Chemistry_rhs, only: gas
  use cantera
  use FLINT_cantera_load
# endif
  implicit none
  integer, parameter :: Tend=2000, Tstart=100
  real(8) :: T, rho, R
  real(8), allocatable :: droic(:), rhoi(:)
  real(8), allocatable :: wdot_explicit(:,:), wdot_cantera(:,:)
  real(8) :: time1, time2
  integer :: i, j, err
  character(32) :: mech_name

  !-------------------------------------------------------------------------------------------------
  ! WD
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p wdot/WD')
  err = read_idealgas_thermo('../../database/WD')
  err = read_chemistry( folder='../../database/WD', mech_name=mech_name )
# if defined(CANTERA)
  call load_phase(gas, '../../database/WD/WD.yaml')
# endif
  call Assign_Mechanism(mech_name)

  allocate(rhoi(1:ns))
  allocate(droic(1:ns))
  allocate(wdot_explicit(ns,Tstart:Tend))
  allocate(wdot_cantera(ns,Tstart:Tend))

  rhoi = 1d-20
  rhoi(1) = 0.2
  rhoi(2) = 0.8
  R = f_Rtot(rhoi)
  rho = sum(rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      droic = 0.d0   ! rows outside the rate tables stay zero (TSR-Rich-31's tables start at 500 K)
      if (i >= T_tab_min .and. i < T_tab_max) call Chemistry_Source ( rhoi, T, droic )
      wdot_explicit(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'WD explicit time =', time2-time1

# if defined(CANTERA)

  call setState_TRY(gas, T, rho, rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      call setTemperature(gas, T)
      call getNetProductionRates(gas, droic)
      droic = droic*wm_tab
      wdot_cantera(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'WD Cantera time =', time2-time1

# endif

  open(100, file='wdot/WD/wdot-explicit.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
  enddo
  close(100)
# if defined(CANTERA)
  open(200, file='wdot/WD/wdot-cantera.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
  enddo
  close(200)
# endif

  deallocate(wdot_cantera); deallocate(wdot_explicit)
  deallocate(droic)
  deallocate(rhoi)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(species_names)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! TROYES
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p wdot/Troyes')
  err = read_idealgas_thermo('../../database/Troyes/')
  err = read_chemistry( folder='../../database/Troyes/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../../database/Troyes/troyes.yaml')
# endif
  call Assign_Mechanism(mech_name)

  allocate(rhoi(1:ns))
  allocate(droic(1:ns))
  allocate(wdot_explicit(ns,Tstart:Tend))
  allocate(wdot_cantera(ns,Tstart:Tend))

  rhoi = 1d0/ns
  R = f_Rtot(rhoi)
  rho = sum(rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      droic = 0.d0   ! rows outside the rate tables stay zero (TSR-Rich-31's tables start at 500 K)
      if (i >= T_tab_min .and. i < T_tab_max) call Chemistry_Source ( rhoi, T, droic )
      wdot_explicit(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'Troyes explicit time =', time2-time1

# if defined (CANTERA)

  call setState_TRY(gas, T, rho, rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      call setTemperature(gas, T)
      call getNetProductionRates(gas, droic)
      droic = droic*wm_tab
      wdot_cantera(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'Troyes Cantera time =', time2-time1

# endif

  open(100, file='wdot/Troyes/wdot-explicit.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
  enddo
  close(100)
# if defined (CANTERA)
  open(200, file='wdot/Troyes/wdot-cantera.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
  enddo
  close(200)
# endif

  deallocate(wdot_cantera); deallocate(wdot_explicit)
  deallocate(droic)
  deallocate(rhoi)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! ECKER
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p wdot/Ecker')
  err = read_idealgas_thermo('../../database/Ecker/')
  err = read_chemistry( folder='../../database/Ecker/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../../database/Ecker/ecker.yaml')
# endif
  call Assign_Mechanism(mech_name)

  allocate(rhoi(1:ns))
  allocate(droic(1:ns))
  allocate(wdot_explicit(ns,Tstart:Tend))
  allocate(wdot_cantera(ns,Tstart:Tend))

  rhoi = 1d-20
  rhoi(2) = 0.00606
  rhoi(7) = 0.00365
  rhoi(9) = 0.00861
  rhoi(11) = 0.00025
  rhoi(14) = 0.98143
  R = f_Rtot(rhoi)
  rho = sum(rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      droic = 0.d0   ! rows outside the rate tables stay zero (TSR-Rich-31's tables start at 500 K)
      if (i >= T_tab_min .and. i < T_tab_max) call Chemistry_Source ( rhoi, T, droic )
      wdot_explicit(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'Ecker explicit time =', time2-time1

# if defined (CANTERA)

  call setState_TRY(gas, T, rho, rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      call setTemperature(gas, T)
      call getNetProductionRates(gas, droic)
      droic = droic*wm_tab
      wdot_cantera(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'Ecker Cantera time =', time2-time1

# endif

  open(100, file='wdot/Ecker/wdot-explicit.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
  enddo
  close(100)
# if defined (CANTERA)
  open(200, file='wdot/Ecker/wdot-cantera.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
  enddo
  close(200)
# endif

  deallocate(wdot_cantera); deallocate(wdot_explicit)
  deallocate(droic)
  deallocate(rhoi)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  if (allocated(species_names)) deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! CROSS
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p wdot/Cross')
  err = read_idealgas_thermo('../../database/Cross/')
  err = read_chemistry( folder='../../database/Cross/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../../database/Cross/cross.yaml')
# endif
  call Assign_Mechanism(mech_name)

  allocate(rhoi(1:ns))
  allocate(droic(1:ns))
  allocate(wdot_explicit(ns,Tstart:Tend))
  allocate(wdot_cantera(ns,Tstart:Tend))

  rhoi = 1d0/ns
  R = f_Rtot(rhoi)
  rho = sum(rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      droic = 0.d0   ! rows outside the rate tables stay zero (TSR-Rich-31's tables start at 500 K)
      if (i >= T_tab_min .and. i < T_tab_max) call Chemistry_Source ( rhoi, T, droic )
      wdot_explicit(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'Cross explicit time =', time2-time1

# if defined (CANTERA)

  call setState_TRY(gas, T, rho, rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      call setTemperature(gas, T)
      call getNetProductionRates(gas, droic)
      droic = droic*wm_tab
      wdot_cantera(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'Cross Cantera time =', time2-time1

# endif

  open(100, file='wdot/Cross/wdot-explicit.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
  enddo
  close(100)
# if defined (CANTERA)
  open(200, file='wdot/Cross/wdot-cantera.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
  enddo
  close(200)
# endif

  deallocate(wdot_cantera); deallocate(wdot_explicit)
  deallocate(droic)
  deallocate(rhoi)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  if (allocated(species_names)) deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! SMOOKE
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p wdot/Smooke')
  err = read_idealgas_thermo('../../database/Smooke/')
  err = read_chemistry( folder='../../database/Smooke/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../../database/Smooke/smooke.yaml')
# endif
  call Assign_Mechanism(mech_name)

  allocate(rhoi(1:ns))
  allocate(droic(1:ns))
  allocate(wdot_explicit(ns,Tstart:Tend))
  allocate(wdot_cantera(ns,Tstart:Tend))

  rhoi = 1d0/ns
  R = f_Rtot(rhoi)
  rho = sum(rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      droic = 0.d0   ! rows outside the rate tables stay zero (TSR-Rich-31's tables start at 500 K)
      if (i >= T_tab_min .and. i < T_tab_max) call Chemistry_Source ( rhoi, T, droic )
      wdot_explicit(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'Smooke explicit time =', time2-time1

# if defined (CANTERA)

  call setState_TRY(gas, T, rho, rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      call setTemperature(gas, T)
      call getNetProductionRates(gas, droic)
      droic = droic*wm_tab
      wdot_cantera(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'Smooke Cantera time =', time2-time1

# endif

  open(100, file='wdot/Smooke/wdot-explicit.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
  enddo
  close(100)
# if defined (CANTERA)
  open(200, file='wdot/Smooke/wdot-cantera.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
  enddo
  close(200)
# endif

  deallocate(wdot_cantera); deallocate(wdot_explicit)
  deallocate(droic)
  deallocate(rhoi)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  if (allocated(species_names)) deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! CORIA-CNRS
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p wdot/CORIA')
  err = read_idealgas_thermo('../../database/CORIA/')
  err = read_chemistry( folder='../../database/CORIA/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../../database/CORIA/coria.yaml')
# endif
  call Assign_Mechanism(mech_name)

  allocate(rhoi(1:ns))
  allocate(droic(1:ns))
  allocate(wdot_explicit(ns,Tstart:Tend))
  allocate(wdot_cantera(ns,Tstart:Tend))

  rhoi = 1d0/ns
  R = f_Rtot(rhoi)
  rho = sum(rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      droic = 0.d0   ! rows outside the rate tables stay zero (TSR-Rich-31's tables start at 500 K)
      if (i >= T_tab_min .and. i < T_tab_max) call Chemistry_Source ( rhoi, T, droic )
      wdot_explicit(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'CORIA explicit time =', time2-time1

# if defined (CANTERA)

  call setState_TRY(gas, T, rho, rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      call setTemperature(gas, T)
      call getNetProductionRates(gas, droic)
      droic = droic*wm_tab
      wdot_cantera(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'CORIA Cantera time =', time2-time1

# endif

  open(100, file='wdot/CORIA/wdot-explicit.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
  enddo
  close(100)
# if defined (CANTERA)
  open(200, file='wdot/CORIA/wdot-cantera.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
  enddo
  close(200)
# endif

  deallocate(wdot_cantera); deallocate(wdot_explicit)
  deallocate(droic)
  deallocate(rhoi)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  if (allocated(species_names)) deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! TSR-CDF-13
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p wdot/TSR-CDF-13')
  err = read_idealgas_thermo('../../database/TSR-CDF-13/')
  err = read_chemistry( folder='../../database/TSR-CDF-13/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../../database/TSR-CDF-13/TSR-CDF-13.yaml')
# endif
  call Assign_Mechanism(mech_name)

  allocate(rhoi(1:ns))
  allocate(droic(1:ns))
  allocate(wdot_explicit(ns,Tstart:Tend))
  allocate(wdot_cantera(ns,Tstart:Tend))

  rhoi = 1d0/ns
  R = f_Rtot(rhoi)
  rho = sum(rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      droic = 0.d0   ! rows outside the rate tables stay zero (TSR-Rich-31's tables start at 500 K)
      if (i >= T_tab_min .and. i < T_tab_max) call Chemistry_Source ( rhoi, T, droic )
      wdot_explicit(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'TSR-CDF-13 explicit time =', time2-time1

# if defined (CANTERA)

  call setState_TRY(gas, T, rho, rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      call setTemperature(gas, T)
      call getNetProductionRates(gas, droic)
      droic = droic*wm_tab
      wdot_cantera(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'TSR-CDF-13 Cantera time =', time2-time1

# endif

  open(100, file='wdot/TSR-CDF-13/wdot-explicit.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
  enddo
  close(100)
# if defined (CANTERA)
  open(200, file='wdot/TSR-CDF-13/wdot-cantera.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
  enddo
  close(200)
# endif

  deallocate(wdot_cantera); deallocate(wdot_explicit)
  deallocate(droic)
  deallocate(rhoi)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  if (allocated(species_names)) deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! TSR-GP-24
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p wdot/TSR-GP-24')
  err = read_idealgas_thermo('../../database/TSR-GP-24/')
  err = read_chemistry( folder='../../database/TSR-GP-24/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../../database/TSR-GP-24/TSR-GP-24.yaml')
# endif
  call Assign_Mechanism(mech_name)

  allocate(rhoi(1:ns))
  allocate(droic(1:ns))
  allocate(wdot_explicit(ns,Tstart:Tend))
  allocate(wdot_cantera(ns,Tstart:Tend))

  rhoi = 1d0/ns
  R = f_Rtot(rhoi)
  rho = sum(rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      droic = 0.d0   ! rows outside the rate tables stay zero (TSR-Rich-31's tables start at 500 K)
      if (i >= T_tab_min .and. i < T_tab_max) call Chemistry_Source ( rhoi, T, droic )
      wdot_explicit(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'TSR-GP-24 explicit time =', time2-time1

# if defined (CANTERA)

  call setState_TRY(gas, T, rho, rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      call setTemperature(gas, T)
      call getNetProductionRates(gas, droic)
      droic = droic*wm_tab
      wdot_cantera(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'TSR-GP-24 Cantera time =', time2-time1

# endif

  open(100, file='wdot/TSR-GP-24/wdot-explicit.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
  enddo
  close(100)
# if defined (CANTERA)
  open(200, file='wdot/TSR-GP-24/wdot-cantera.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
  enddo
  close(200)
# endif

  deallocate(wdot_cantera); deallocate(wdot_explicit)
  deallocate(droic)
  deallocate(rhoi)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  if (allocated(species_names)) deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! TSR-Rich-31
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p wdot/TSR-Rich-31')
  err = read_idealgas_thermo('../../database/TSR-Rich-31/')
  err = read_chemistry( folder='../../database/TSR-Rich-31/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../../database/TSR-Rich-31/TSR-Rich-31.yaml')
# endif
  call Assign_Mechanism(mech_name)

  allocate(rhoi(1:ns))
  allocate(droic(1:ns))
  allocate(wdot_explicit(ns,Tstart:Tend))
  allocate(wdot_cantera(ns,Tstart:Tend))

  rhoi = 1d0/ns
  R = f_Rtot(rhoi)
  rho = sum(rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      droic = 0.d0   ! rows outside the rate tables stay zero (TSR-Rich-31's tables start at 500 K)
      if (i >= T_tab_min .and. i < T_tab_max) call Chemistry_Source ( rhoi, T, droic )
      wdot_explicit(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'TSR-Rich-31 explicit time =', time2-time1

# if defined (CANTERA)

  call setState_TRY(gas, T, rho, rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      call setTemperature(gas, T)
      call getNetProductionRates(gas, droic)
      droic = droic*wm_tab
      wdot_cantera(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'TSR-Rich-31 Cantera time =', time2-time1

# endif

  open(100, file='wdot/TSR-Rich-31/wdot-explicit.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
  enddo
  close(100)
# if defined (CANTERA)
  open(200, file='wdot/TSR-Rich-31/wdot-cantera.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
  enddo
  close(200)
# endif

  deallocate(wdot_cantera); deallocate(wdot_explicit)
  deallocate(droic)
  deallocate(rhoi)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  if (allocated(species_names)) deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

  !-------------------------------------------------------------------------------------------------
  ! Pelucchi
  !-------------------------------------------------------------------------------------------------

  call execute_command_line('mkdir -p wdot/Pelucchi')
  err = read_idealgas_thermo('../../database/Pelucchi/')
  err = read_chemistry( folder='../../database/Pelucchi/', mech_name=mech_name )
# if defined (CANTERA)
  call load_phase(gas, '../../database/Pelucchi/pelucchi.yaml')
# endif
  call Assign_Mechanism(mech_name)

  allocate(rhoi(1:ns))
  allocate(droic(1:ns))
  allocate(wdot_explicit(ns,Tstart:Tend))
  allocate(wdot_cantera(ns,Tstart:Tend))

  rhoi = 1d0/ns
  R = f_Rtot(rhoi)
  rho = sum(rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      droic = 0.d0   ! rows outside the rate tables stay zero (TSR-Rich-31's tables start at 500 K)
      if (i >= T_tab_min .and. i < T_tab_max) call Chemistry_Source ( rhoi, T, droic )
      wdot_explicit(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'Pelucchi explicit time =', time2-time1

# if defined (CANTERA)

  call setState_TRY(gas, T, rho, rhoi)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      call setTemperature(gas, T)
      call getNetProductionRates(gas, droic)
      droic = droic*wm_tab
      wdot_cantera(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'Pelucchi Cantera time =', time2-time1

# endif

  open(100, file='wdot/Pelucchi/wdot-explicit.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(100,*) dble(i), (wdot_explicit(j,i),j=1,ns)
  enddo
  close(100)
# if defined (CANTERA)
  open(200, file='wdot/Pelucchi/wdot-cantera.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(200,*) dble(i), (wdot_cantera(j,i),j=1,ns)
  enddo
  close(200)
# endif

  deallocate(wdot_cantera); deallocate(wdot_explicit)
  deallocate(droic)
  deallocate(rhoi)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  if (allocated(species_names)) deallocate(species_names)
  if (allocated(elements_names)) deallocate(elements_names)
  if (allocated(species_composition)) deallocate(species_composition)
  call free_chemistry_data()

end program test