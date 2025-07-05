program testN0
  use U_Lib_Thermodynamic
  use U_IO_Table
  use cantera
  use U_cantera_load
  implicit none
  type(phase_t) :: gas
  real(8) :: rhoi(9), T, rho, R
  real(8) :: cp_native, cp_cantera
  real(8) :: lam_cantera, lam_native, mi_native, nu_cantera
  real(8) :: time1, time2
  integer :: i, j

  call Read_IdealGas_Properties()
  call load_phase(gas, 'WD.yaml')

  rhoi = 0.d0
  rhoi(1) = 0.25d0 ! CH4
  rhoi(2) = 0.15d0 ! O2
  rhoi(3) = 0.30d0 ! CO2
  rhoi(5) = 0.40d0 ! CO
  call co_rotot_Rtot(rhoi,rho,R)

  call cpu_time(time1)
  do j = 1, 1000
    do i = 1,2000
      T = dble(i)
      cp_native = co_cp(rhoi,T,rho)
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) time2-time1

  call setState_TRY(gas, T, rho, rhoi)

  call cpu_time(time1)
  do j = 1, 1000
    do i = 1,2000
      T = dble(i)
      call setTemperature(gas, T)
      cp_cantera = cp_mass(gas)
    enddo
  enddo
  call cpu_time(time2)

  write(*,*) time2-time1

  write(*,*) abs(cp_cantera-cp_native)/cp_cantera*100d0

end program testN0