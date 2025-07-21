module globH2_mod
  implicit none
contains

  ! Frolov: hydrogen/air with 3 species, one reaction 
  subroutine Frolov(roi,temp,omegadot,rotot)
    use U_Lib_Thermodynamic
    use U_Lib_Chemistry_data
    implicit none
    integer :: is, T_i, Tint(2)
    real(8), intent(in)    :: roi(nsc)
    real(8), intent(in)    :: temp 
    real(8), intent(out)   :: omegadot(nsc)
    real(8), intent(in)    :: rotot
    ! Local
    real(8) :: coi(nsc+1), Tdiff
    real(8) :: p
    real(8) :: prod1

    !--------------------------------------------------------------

    do is = 1, nsc
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
      if (coi(is).lt.1d-12) coi(is) = 0d0
    enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    p = sum(roi) * sum(coi * Ri_tab) * temp

    ! 2 H2 + O2 --> 2 H2O
    prod1 = -0.5d0 * 8.d11 * ((p / 101325d0)**(-1.15d0)) * coi(2) **2 * coi(5) * exp(-10000d0/temp)

    ! Chemical Source Terms
    omegadot = 0d0
    omegadot(2)= Wm_tab(2) * 2 * prod1        !H2    [kg/(m3*s)]
    omegadot(5)= Wm_tab(5)* prod1             !O2    [kg/(m3*s)]
    omegadot(3)= Wm_tab(3)*(-2 * prod1)       !H20   [kg/(m3*s)]

  end subroutine Frolov


  ! Nassini 
  subroutine Nassini_4(roi,temp,omegadot,rotot)
    use U_Lib_Thermodynamic
    use U_Lib_Chemistry_data
    implicit none
    integer :: is, T_i, Tint(2)
    real(8), intent(in) :: roi(nsc)
    real(8), intent(in)    :: temp 
    real(8), intent(out)   :: omegadot(nsc)
    real(8), intent(in)    :: rotot
    ! Local
    real(8) :: coi(nsc+1), Tdiff
    real(8) :: prod1

    !--------------------------------------------------------------

    do is = 1, nsc
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
      if (coi(is).lt.1d-12) coi(is) = 0d0
    enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! H2 + 0.5 O2 <--> H2O
    prod1 = comp_ch_tabT(1,kf_tab,Tint,Tdiff)*coi(3)*coi(1)-comp_ch_tabT(1,kb_tab,Tint,Tdiff)*coi(2)

    ! Chemical Source Terms
    omegadot = 0d0
    omegadot(1)= Wm_tab(1) * (-0.5d0 * prod1)        !O2    [kg/(m3*s)]
    omegadot(2)= Wm_tab(2)* 1.d0 * prod1             !H2O   [kg/(m3*s)]
    omegadot(3)= Wm_tab(3)*(-1.d0 * prod1)           !H2    [kg/(m3*s)]

  end subroutine Nassini_4
end module globH2_mod