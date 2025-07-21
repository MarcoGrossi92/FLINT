module U_Lib_Chemistry_rhs
  use OSLo
  use U_Lib_Chemistry_wdot
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
    use U_Lib_Thermodynamic
    implicit none
    integer, intent(in)  :: nz
    real(8), intent(in)  :: time
    real(8), intent(in)  :: Z(nz)
    real(8), intent(out) :: F(nz)
    ! Local
    real(8) :: roi(nz-1)
    real(8) :: rotot,Rgas
    real(8) :: T,p
    real(8) :: dpc,dTc,droic(nz-1),wdot(nz-1)
    real(8) :: eiroi,rho_cv
    integer :: s

    roi(1:nsc) = Z(1:nsc)
    call co_rotot_Rtot ( roi, rotot, Rgas )
    T = Z(nz)

    ! Concentration source terms due to chemistry droic(i) = omegadot(i)
    call Chemistry_Source ( roi, T, droic, rotot )

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
    use U_Lib_Thermodynamic
    use cantera
    implicit none
    integer, intent(in)  :: nz
    real(8), intent(in)  :: time
    real(8), intent(in)  :: Z(nz)
    real(8), intent(out) :: F(nz)
    ! Local
    real(8) :: roi(nz-1)
    real(8) :: rotot,Rgas
    real(8) :: T,p
    real(8) :: dpc,dTc,droic(nz-1),wdot(nz-1)
    real(8) :: eiroi,rho_cv
    integer :: s

    roi(1:nsc) = Z(1:nsc)
    call co_rotot_Rtot ( roi, rotot, Rgas )
    T = Z(nz)

    ! Concentration source terms due to chemistry droic(i) = omegadot(i)
    call setState_TRY(gas, t, rotot, roi)
    call getNetProductionRates(gas, wdot)
    droic = wdot*wm_tab

    eiroi = 0.d0; rho_cv = 0.d0
    do s = 1, nsc
      eiroi = eiroi + ( comp_ms_tabT(T,s,h_tab) - Ri_tab(s)*T ) * droic(s)
      rho_cv = rho_cv + roi(s)*( comp_ms_tabT(T,s,cp_tab) - Ri_tab(s) )
    enddo

    F(1:nsc) = droic
    F(nz) = -eiroi / rho_cv

  end subroutine rhs_cantera
# endif

end module U_Lib_Chemistry_rhs
