module globH2_mod
  implicit none
contains

  ! Frolov: hydrogen/air with 3 species, one reaction 
  subroutine Frolov(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    integer :: is, T_i, Tint(2)
    real(8), intent(inout)    :: roi(ns)
    real(8), intent(in)    :: temp 
    real(8), intent(out)   :: omegadot(ns)
    ! Local
    real(8) :: coi(ns), Tdiff
    real(8) :: p
    real(8) :: prod1

    !--------------------------------------------------------------

    do is = 1, ns
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
      if (coi(is).lt.1d-12) coi(is) = 0d0
    enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    p = sum(roi * Ri_tab) * temp

    ! 2 H2 + O2 --> 2 H2O
    prod1 = -0.5d0 * 8.d11 * ((p / 101325d0)**(-1.15d0)) * coi(3) **2 * coi(1) * exp(-10000d0/temp)

    ! Chemical Source Terms
    omegadot = 0d0
    omegadot(1)= Wm_tab(1)* prod1             !O2    [kg/(m3*s)]
    omegadot(2)= Wm_tab(2)*(-2 * prod1)       !H20   [kg/(m3*s)]
    omegadot(3)= Wm_tab(3) * 2 * prod1        !H2    [kg/(m3*s)]
    

  end subroutine Frolov


  ! Frolov_nopressure: the same single-step H2/O2 reaction as Frolov above, but
  ! without the (p/p0)^-1.15 pressure correction, and with the reverse rate taken
  ! from kb_tab instead of being omitted.
  !
  !     2 H2 + O2  <=>  2 H2O
  !     forward Arrhenius: A = 8e11, b = 0, Ea = 8.314e7 J/kmol
  !
  ! Reverse rate comes from the equilibrium constant computed by the chemistry
  ! loader and supplied through kb_tab (same convention as ONERA-7, general, ...).
  !
  ! Species ordering (matches the YAML phase `species: [O2,H2O,H2,N2]`):
  !     1: O2   2: H2O   3: H2   4: N2 (inert; row is all zero)
  !
  ! Hand-rolled to bypass the data-driven `general` dispatcher in
  ! FLINT_Lib_Chemistry_wdot. On the THOR Scalability case `general` was the #2
  ! hotspot at ~9000 s of CPU per rank -- dominated by dispatch overhead (outer
  ! reaction loop, branch on `ni1 /= 0`, integer exponent via `**nint(...)`).
  ! This straight-line version eliminates all of that for the single reaction.
  subroutine Frolov_nopressure(roi, temp, omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in)    :: temp
    real(8), intent(out)   :: omegadot(ns)
    ! Local
    real(8) :: coi(ns), Tdiff
    integer :: is, T_i, Tint(2)
    real(8) :: kf, kb, net_rate

    do is = 1, ns
      coi(is) = roi(is) / Wm_tab(is)   ! kmol/m^3
    end do

    T_i     = int(temp)
    Tdiff   = temp - T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    kf = f_kf(1, Tint, Tdiff)
    kb = f_kb(1, Tint, Tdiff)

    ! Reaction 1:  2 H2 + O2 <=> 2 H2O
    ! prodf = kf * c(3)^2 * c(1) ;  prodb = kb * c(2)^2
    net_rate = kf * coi(3) * coi(3) * coi(1) - kb * coi(2) * coi(2)

    ! Stoichiometric coefficients (signed, product - reactant):
    !   O2 (1): -1     H2O (2): +2     H2 (3): -2     N2 (4): 0
    omegadot(1) = -      Wm_tab(1) * net_rate
    omegadot(2) =  2.d0 * Wm_tab(2) * net_rate
    omegadot(3) = -2.d0 * Wm_tab(3) * net_rate
    omegadot(4) = 0.d0
  end subroutine Frolov_nopressure


  ! Analytical Jacobian of Frolov_nopressure.
  !
  !   dwdr(i,j) = d omegadot(i) / d roi(j)     [1/s]
  !   dwdT(i)   = d omegadot(i) / d T          [kg/(m^3 s K)]
  !
  ! Same table-derivative trick as ONERA_7_jac: kf and kb are linearly
  ! interpolated in T from kf_tab / kb_tab, so dk/dT is a free table lookup.
  subroutine Frolov_nopressure_jac(roi, temp, dwdr, dwdT)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(in)  :: roi(ns), temp
    real(8), intent(out) :: dwdr(ns, ns)
    real(8), intent(out) :: dwdT(ns)
    ! Local
    real(8) :: coi(ns), Tdiff
    integer :: is, T_i, j
    real(8) :: kf, kb, dkf_dT, dkb_dT
    real(8) :: dnet_dc(ns), dnet_dT_r
    real(8) :: dwdr_c(ns, ns)

    do is = 1, ns
      coi(is) = roi(is) / Wm_tab(is)
    end do

    T_i   = int(temp)
    Tdiff = temp - T_i

    kf     = kf_tab(T_i  , 1) + Tdiff * (kf_tab(T_i+1, 1) - kf_tab(T_i, 1))
    kb     = kb_tab(T_i  , 1) + Tdiff * (kb_tab(T_i+1, 1) - kb_tab(T_i, 1))
    dkf_dT = kf_tab(T_i+1, 1) - kf_tab(T_i, 1)
    dkb_dT = kb_tab(T_i+1, 1) - kb_tab(T_i, 1)

    ! d(net_rate)/d coi(j) for reaction 1:
    !   net_rate = kf * c(3)^2 * c(1) - kb * c(2)^2
    !   d/dc(1) =        kf * c(3)^2
    !   d/dc(2) = - 2 * kb * c(2)
    !   d/dc(3) =   2 * kf * c(3) * c(1)
    !   d/dc(4) = 0  (N2 inert)
    dnet_dc(:) = 0.d0
    dnet_dc(1) =          kf * coi(3) * coi(3)
    dnet_dc(2) = -2.d0 *  kb * coi(2)
    dnet_dc(3) =  2.d0 *  kf * coi(3) * coi(1)
    dnet_dT_r  = dkf_dT * coi(3) * coi(3) * coi(1) &
               - dkb_dT * coi(2) * coi(2)

    ! Apply stoichiometry to get d omegadot/d coi(j) and d omegadot/d T.
    dwdr_c = 0.d0
    dwdr_c(1, :) = -      Wm_tab(1) * dnet_dc
    dwdr_c(2, :) =  2.d0 * Wm_tab(2) * dnet_dc
    dwdr_c(3, :) = -2.d0 * Wm_tab(3) * dnet_dc
    ! row 4 (N2) stays zero by construction.

    dwdT(:) = 0.d0
    dwdT(1) = -      Wm_tab(1) * dnet_dT_r
    dwdT(2) =  2.d0 * Wm_tab(2) * dnet_dT_r
    dwdT(3) = -2.d0 * Wm_tab(3) * dnet_dT_r

    ! Convert d/d coi(j) -> d/d roi(j): divide column j by Wm_tab(j).
    do j = 1, ns
      dwdr(:, j) = dwdr_c(:, j) / Wm_tab(j)
    end do
  end subroutine Frolov_nopressure_jac


  ! Nassini
  subroutine Nassini_4(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    integer :: is, T_i, Tint(2)
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in)    :: temp 
    real(8), intent(out)   :: omegadot(ns)
    ! Local
    real(8) :: coi(ns), Tdiff
    real(8) :: prod1, prod2

    !--------------------------------------------------------------

    do is = 1, ns
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
      if (coi(is).lt.1d-12) coi(is) = 0d0
    enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! H2 + 0.5 O2 --> H2O
    prod1 = comp_ch_tabT(1,kf_tab,Tint,Tdiff)*coi(3)*coi(1)

    ! H2O --> H2 + 0.5 O2
    prod2 = comp_ch_tabT(2,kf_tab,Tint,Tdiff)*coi(2)

    ! Chemical Source Terms
    omegadot = 0d0
    omegadot(1)= Wm_tab(1) * (-0.5d0 * prod1 + 0.5 * prod2)        !O2    [kg/(m3*s)]
    omegadot(2)= Wm_tab(2)* (1.d0 * prod1 - 1 * prod2)             !H2O   [kg/(m3*s)]
    omegadot(3)= Wm_tab(3)*(-1.d0 * prod1 + 1 * prod2)             !H2    [kg/(m3*s)]

  end subroutine Nassini_4
end module globH2_mod
