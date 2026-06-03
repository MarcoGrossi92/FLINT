
!> \author Marco Grossi and Andrea Giacomi
!> \date 2026
program test
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  implicit none
  real(8) :: p_, h_
  real(8) :: p_test, T_test, h_test, T_back
  real(8) :: T_min_global, T_max_global
  real(8) :: err_max, err_curr
  integer :: err, i_test, Np_tab, Nh_tab, i_valid, ii

  ! Load data
  FLINT_phase_prefix = ''
  err = read_realfluid_thermo('real-fluid/INPUT/')
  write(*,*) 'Error code from reading real fluid thermo:    ', err
  err = read_realfluid_transport('real-fluid/INPUT/')
  write(*,*) 'Error code from reading real fluid transport: ', err

  ! Write tables size and range
  write(*,*) 'Real fluid thermo tables loaded with size:    ', size(rho_tab,1), size(rho_tab,2)
  write(*,*) 'Real fluid transport tables loaded with size: ', size(mi_tab2D,1), size(mi_tab2D,2)
  write(*,*) ''
  write(*,*) 'Pressure range: ', pmin*1e-5, ' bar, to ', pmin*1e-5 + (size(rho_tab,1)-1)*deltap*1e-5, ' bar, delta ', deltap*1e-5, ' bar'
  write(*,*) 'Enthalpy range: ', hmin, ' J/kg, to ', hmin + (size(rho_tab,2)-1)*deltah, ' J/kg, delta ', deltah, ' J/kg'
  write(*,*) ''

  ! Test interpolation function works
  p_ = 100.0d5
  h_ = 1.0d5
  write(*,*) ''
  write(*,*) 'Testing real fluid thermo properties at p = ', p_, ' Pa and h = ', h_, ' J/kg'
  write(*,*) 'Density: ', ph2vars(p_,h_,rho_tab), ' kg/m^3'
  write(*,*) 'Temperature: ', ph2vars(p_, h_, T_tab), ' K'
  write(*,*) 'Viscosity:    ', ph2vars(p_, h_, mi_tab2D), ' Pa*s'
  write(*,*) 'Conductivity: ', ph2vars(p_, h_, k_tab2D), ' W/(m*K)'


  ! Test ph2pT: build inverse table (p, T) -> h
  write(*,*) ''
  err = ph2pT()
  write(*,*) 'Error code from ph2pT: ', err
  if (err >= 1) then
    write(*,*) 'Fatal error from ph2pT, skipping validation tests'
    stop 1
  end if

  Np_tab = ubound(T_tab, 1)
  Nh_tab = ubound(T_tab, 2)

  T_min_global = minval( max(Tmin2(0:Np_tab-1), Tmin2(1:Np_tab)) )
  T_max_global = maxval( min(Tmin2(0:Np_tab-1) + Nh_tab*deltaT(0:Np_tab-1), &
                             Tmin2(1:Np_tab)   + Nh_tab*deltaT(1:Np_tab)) )

  write(*,'(A4,A12,A12,A18,A14,A14)') '  k ', '  T_in [K] ', '  p [bar] ', '  h = pT2h [J/kg]', '  T_back [K]', '  err [K]'

  err_max = 0.0d0
  do i_test = 0, 10
    T_test = T_min_global + (T_max_global - T_min_global) * dble(i_test) / 10.0d0

    ! find a row pair (i_valid, i_valid+1) whose T-ranges both contain T_test
    i_valid = -1
    do ii = 0, Np_tab - 1
      if (T_test >= Tmin2(ii)   .and. T_test <= Tmin2(ii)   + Nh_tab*deltaT(ii)   .and. &
          T_test >= Tmin2(ii+1) .and. T_test <= Tmin2(ii+1) + Nh_tab*deltaT(ii+1)) then
        i_valid = ii
        exit
      end if
    end do

    p_test = pmin + (dble(i_valid) + 0.5d0) * deltap   ! midpoint between p_{i_valid} and p_{i_valid+1}
    h_test = pT2h(p_test, T_test)
    T_back = ph2vars(p_test, h_test, T_tab)
    err_curr = T_back - T_test
    if (abs(err_curr) > err_max) err_max = abs(err_curr)
    write(*,'(I4,F12.2,F12.2,ES18.5,F14.4,ES14.3)') i_test, T_test, p_test*1.d-5, h_test, T_back, err_curr
  end do
  write(*,'(A,ES12.3,A)') ' Max round-trip error: ', err_max, ' K'


end program test