
!> \author Marco Grossi and Andrea Giacomi
!> \date 2026
program test
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  implicit none
  real(8) :: p_, h_
  integer :: err

  FLINT_phase_prefix = ''
  err = read_realfluid_thermo('real-fluid/INPUT/')
  write(*,*) 'Error code from reading real fluid thermo: ', err

  ! Write tables size and range
  write(*,*) 'Real fluid thermo tables loaded with size: ', size(rho_tab,1), size(rho_tab,2)
  write(*,*) 'Pressure range: ', pmin*1e-5, ' bar, to ', pmin*1e-5 + (size(rho_tab,1)-1)*deltap*1e-5, ' bar, delta ', deltap*1e-5, ' bar'
  write(*,*) 'Enthalpy range: ', hmin, ' J/kg, to ', hmin + (size(rho_tab,2)-1)*deltah, ' J/kg, delta ', deltah, ' J/kg'

  p_ = 100.0d5
  h_ = 1.0d5
  write(*,*) 'Testing real fluid thermo properties at p = ', p_, ' Pa and h = ', h_, ' J/kg'
  write(*,*) 'Density: ', ph2vars(p_,h_,rho_tab), ' kg/m^3'
  write(*,*) 'Temperature: ', ph2vars(p_, h_, T_tab), ' K'


end program test