module FLINT_Lib_Chemistry_rhs
  use OSLo
  use FLINT_Lib_Chemistry_wdot
# if defined (CANTERA)
  use cantera
# endif
  implicit none
  logical :: reactive
# if defined(CANTERA)
  type(phase_t) :: gas
# endif

contains
    
  subroutine rhs_native ( nz, time, Z, F )
    use FLINT_Lib_Thermodynamic
    implicit none
    integer, intent(in)  :: nz
    real(8), intent(in)  :: time
    real(8), intent(in)  :: Z(nz)
    real(8), intent(out) :: F(nz)
    ! Local
    real(8) :: roi(nz-1)
    real(8) :: T
    real(8) :: droic(nz-1)
    real(8) :: eiroi,rho_cv
    integer :: s

    roi(1:nsc) = Z(1:nsc)
    T = Z(nz)
    if (T < 0.d0 .or. T > 10000d0) then
      F(:) = -1.d0
      return 1
    end if

    ! Avoid negative rho_i
    do s = 1, nsc
      roi(s) = max(roi(s), 0.d0)
    end do

    call Chemistry_Source ( roi, T, droic )

    eiroi = 0.d0; rho_cv = 0.d0
    do s = 1, nsc
      eiroi = eiroi + ( comp_ms_tabT(T,s,h_tab) - Ri_tab(s)*T ) * droic(s)
      rho_cv = rho_cv + roi(s)*( comp_ms_tabT(T,s,cp_tab) - Ri_tab(s) )
    enddo

    F(1:nsc) = droic
    F(nz) = -eiroi / rho_cv

  end subroutine rhs_native

# if defined (CANTERA)
  subroutine rhs_cantera ( nz, time, Z, F )
    use FLINT_Lib_Thermodynamic
    use cantera
    implicit none
    integer, intent(in)  :: nz
    real(8), intent(in)  :: time
    real(8), intent(in)  :: Z(nz)
    real(8), intent(out) :: F(nz)
    ! Local
    real(8) :: roi(nz-1), Y(nz-1)
    real(8) :: T, rho
    real(8) :: droic(nz-1),wdot(nz-1)
    real(8) :: eiroi,rho_cv, e_s, cp_s
    integer :: s

    roi(1:nsc) = Z(1:nsc)
    T = Z(nz)
    if (T < 0.d0 .or. T > 10000d0) then
      F(:) = -1.d0
      return 1
    end if

    ! Avoid negative rho_i
    do s = 1, nsc
      roi(s) = max(roi(s), 0.d0)
    end do
    rho = sum(roi)
    Y = roi / rho

    call setState_TRY(gas, T, rho, Y)
    call getNetProductionRates(gas, wdot)
    droic = wdot * wm_tab

    eiroi = 0.d0
    rho_cv = 0.d0
    do s = 1, nsc
      e_s = comp_ms_tabT(T, s, h_tab) - Ri_tab(s) * T
      cp_s = comp_ms_tabT(T, s, cp_tab) - Ri_tab(s)
      eiroi = eiroi + e_s * droic(s)
      rho_cv = rho_cv + roi(s) * cp_s
    end do

    F(1:nsc) = droic
    F(nz) = -eiroi / rho_cv

  end subroutine rhs_cantera
# endif

end module FLINT_Lib_Chemistry_rhs
