! Test for wdot with thrid-body reactions
program test
  use U_Lib_Thermodynamic
  use U_IO_Table
  use U_IO_chemistry
  use U_Lib_Chemistry_data
  use U_Lib_Chemistry_rhs, only: gas
  use U_Lib_Chemistry_wdot
  use cantera
  use U_cantera_load
  implicit none
  integer, parameter :: Tend=2000, Tstart=100
  real(8) :: T, rho, R
  real(8), allocatable :: droic(:), rhoi(:)
  real(8), allocatable :: wdot_native(:,:), wdot_cantera(:,:)
  real(8) :: time1, time2
  integer :: i, j, err
  character(32) :: mech_name

  !-------------------------------------------------------------------------------------------------
  ! WD
  !-------------------------------------------------------------------------------------------------

  call Read_IdealGas_Properties('WD')
  err = read_chemistry_file( folder='WD', mech_name=mech_name )
  call load_phase(gas, 'WD/WD.yaml')
  call Assign_Mechanism(mech_name)

  allocate(rhoi(1:nsc))
  allocate(droic(1:nsc))
  allocate(wdot_native(nsc,Tstart:Tend))
  allocate(wdot_cantera(nsc,Tstart:Tend))

  rhoi = 1d-20
  rhoi(1) = 0.2
  rhoi(2) = 0.8
  call co_rotot_Rtot(rhoi,rho,R)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      call Chemistry_Source ( rhoi, T, droic, rho )
      wdot_native(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'WD native time =', time2-time1

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

  open(100, file='WD-wdot-native.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(100,*) dble(i), (wdot_native(j,i),j=1,nsc)
  enddo
  close(100)
  open(200, file='WD-wdot-cantera.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(200,*) dble(i), (wdot_cantera(j,i),j=1,nsc)
  enddo
  close(200)

  deallocate(wdot_cantera); deallocate(wdot_native)
  deallocate(droic)
  deallocate(rhoi)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(mi_tab); deallocate(k_tab)
  deallocate(kf_tab); deallocate(kb_tab); deallocate(ni1_tab); deallocate(ni2_tab); deallocate(epsch_tab)

  !-------------------------------------------------------------------------------------------------
  ! TROYES
  !-------------------------------------------------------------------------------------------------

  call Read_IdealGas_Properties('Troyes')
  err = read_chemistry_file( folder='Troyes', mech_name=mech_name )
  call load_phase(gas, 'Troyes/troyes.yaml')
  call Assign_Mechanism(mech_name)

  allocate(rhoi(1:nsc))
  allocate(droic(1:nsc))
  allocate(wdot_native(nsc,Tstart:Tend))
  allocate(wdot_cantera(nsc,Tstart:Tend))

  rhoi = 1d-20
  rhoi(2) = 0.00606
  rhoi(5) = 0.00365
  rhoi(7) = 0.00861
  rhoi(11) = 0.00025
  rhoi(12) = 0.98143
  call co_rotot_Rtot(rhoi,rho,R)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      call Chemistry_Source ( rhoi, T, droic, rho )
      wdot_native(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'Troyes native time =', time2-time1

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

  open(100, file='Troyes-wdot-native.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(100,*) dble(i), (wdot_native(j,i),j=1,nsc)
  enddo
  close(100)
  open(200, file='Troyes-wdot-cantera.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(200,*) dble(i), (wdot_cantera(j,i),j=1,nsc)
  enddo
  close(200)

  deallocate(wdot_cantera); deallocate(wdot_native)
  deallocate(droic)
  deallocate(rhoi)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(mi_tab); deallocate(k_tab)
  deallocate(kf_tab); deallocate(kb_tab); deallocate(ni1_tab); deallocate(ni2_tab); deallocate(epsch_tab)

  !-------------------------------------------------------------------------------------------------
  ! ECKER
  !-------------------------------------------------------------------------------------------------

  call Read_IdealGas_Properties('Ecker')
  err = read_chemistry_file( folder='Ecker', mech_name=mech_name )
  call load_phase(gas, 'Ecker/ecker.yaml')
  call Assign_Mechanism(mech_name)

  allocate(rhoi(1:nsc))
  allocate(droic(1:nsc))
  allocate(wdot_native(nsc,Tstart:Tend))
  allocate(wdot_cantera(nsc,Tstart:Tend))

  rhoi = 1d-20
  rhoi(2) = 0.00606
  rhoi(7) = 0.00365
  rhoi(9) = 0.00861
  rhoi(11) = 0.00025
  rhoi(14) = 0.98143
  call co_rotot_Rtot(rhoi,rho,R)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      call Chemistry_Source ( rhoi, T, droic, rho )
      wdot_native(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'Ecker native time =', time2-time1

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

  open(100, file='Ecker-wdot-native.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(100,*) dble(i), (wdot_native(j,i),j=1,nsc)
  enddo
  close(100)
  open(200, file='Ecker-wdot-cantera.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(200,*) dble(i), (wdot_cantera(j,i),j=1,nsc)
  enddo
  close(200)

  deallocate(wdot_cantera); deallocate(wdot_native)
  deallocate(droic)
  deallocate(rhoi)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(mi_tab); deallocate(k_tab)
  deallocate(kf_tab); deallocate(kb_tab); deallocate(ni1_tab); deallocate(ni2_tab); deallocate(epsch_tab)

  !-------------------------------------------------------------------------------------------------
  ! CROSS
  !-------------------------------------------------------------------------------------------------

  call Read_IdealGas_Properties('Cross')
  err = read_chemistry_file( folder='Cross', mech_name=mech_name )
  call load_phase(gas, 'Cross/cross.yaml')
  call Assign_Mechanism(mech_name)

  allocate(rhoi(1:nsc))
  allocate(droic(1:nsc))
  allocate(wdot_native(nsc,Tstart:Tend))
  allocate(wdot_cantera(nsc,Tstart:Tend))

  rhoi = 1d0/nsc
  call co_rotot_Rtot(rhoi,rho,R)

  call cpu_time(time1)
  do j = 1, 1
    do i = Tstart, Tend
      T = dble(i)
      call Chemistry_Source ( rhoi, T, droic, rho )
      wdot_native(:,i) = droic
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) 'Cross native time =', time2-time1

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

  open(100, file='Cross-wdot-native.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(100,*) dble(i), (wdot_native(j,i),j=1,nsc)
  enddo
  close(100)
  open(200, file='Cross-wdot-cantera.dat', status='replace', form='formatted')
  do i = Tstart, Tend
    write(200,*) dble(i), (wdot_cantera(j,i),j=1,nsc)
  enddo
  close(200)

  deallocate(wdot_cantera); deallocate(wdot_native)
  deallocate(droic)
  deallocate(rhoi)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(mi_tab); deallocate(k_tab)
  deallocate(kf_tab); deallocate(kb_tab); deallocate(ni1_tab); deallocate(ni2_tab); deallocate(epsch_tab)

end program test