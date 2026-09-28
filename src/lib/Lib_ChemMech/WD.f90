  ! WD: Global Westbrook-Dryer mechanism
  ! 5 species & 3 reactions
  module WD_mod
    implicit none
    contains
  subroutine WD(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in)    :: temp 
    real(8), intent(out)   :: omegadot(ns)
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(ns), Tdiff
    real(8) :: prod1,prod2,prod3

    do is = 1, ns
      roi(is) = max(roi(is), 0.d0)
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
    enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! species: [CH4, O2, CO2, H2O, CO]

    ! CH4 + 1.5 O2 => CO + 2 H2O
    prod1 = f_kf(1,Tint,Tdiff)*(coi(1)**0.70)*(coi(2)**0.80)

    ! CO + 0.5 O2 + H2O => CO2 + H2O
    prod2 = f_kf(2,Tint,Tdiff)*(coi(4)**0.5)*coi(5)*(coi(2)**0.25)

    ! CO2 => CO + 0.5 O2
    prod3 = f_kf(3,Tint,Tdiff)*coi(3)
     
    ! Chemical Source Terms
    omegadot = 0d0
    omegadot(1)=Wm_tab(1)*(-prod1)
    omegadot(2)=Wm_tab(2)*(-1.5*prod1-0.5*prod2+0.5*prod3)
    omegadot(3)=Wm_tab(3)*(prod2-prod3)
    omegadot(4)=Wm_tab(4)*(2*prod1)
    omegadot(5)=Wm_tab(5)*(prod1-prod2+prod3)

  end subroutine WD

  !> WD-Andersen: the Westbrook-Dryer steps 1-2 with the CO2 dissociation step (3) written as the
  !> explicit inverse of step 2, concentration orders [CO2]^1 [H2O]^0.5 [O2]^-0.25 (Andersen et al.,
  !> 2009: the pair 2/3 then reaches the equilibrium of CO + 0.5 O2 <-> CO2, which the former law
  !> [CO2]^1.25 did not). Zero-concentration convention (one rule for every FLINT site with a
  !> negative order: pow_order of the general procedure, the guarded steps of Coronetti and JLR):
  !> a negative-order species at zero concentration gives a zero rate, as Cantera does; the power
  !> is evaluated only when the concentration is positive (0**(-0.25) is +Infinity and traps under
  !> -fpe0 / -ffpe-trap). Literals in double precision (the routine is compared with Cantera to 1e-12).
  subroutine Andersen(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in)    :: temp 
    real(8), intent(out)   :: omegadot(ns)
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(ns), Tdiff
    real(8) :: prod1,prod2,prod3,prod4,prod5,prod6

    do is = 1, ns
      roi(is) = max(roi(is), 0.d0)
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
    enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1
    ! species: [CH4, O2, CO2, H2O, CO]

    ! CH4 + 1.5 O2 => CO + 2 H2O
    prod1 = f_kf(1,Tint,Tdiff)*(coi(1)**0.70d0)*(coi(2)**0.80d0)

    ! CO + 0.5 O2 + H2O => CO2 + H2O
    prod2 = f_kf(2,Tint,Tdiff)*(coi(4)**0.5d0)*coi(5)*(coi(2)**0.25d0)

    ! CO2 => CO + 0.5 O2, rate = kf(3) [CO2] [H2O]^0.5 [O2]^-0.25 (Andersen et al. 2009);
    ! zero rate at [O2] = 0 (Cantera's convention for a negative order), the power is not evaluated there
    if (coi(2) > 0d0) then
      prod3 = f_kf(3,Tint,Tdiff)*coi(3)*(coi(4)**0.5d0)*(coi(2)**(-0.25d0))
    else
      prod3 = 0d0
    endif
     
    ! Chemical Source Terms
    omegadot = 0d0
    omegadot(1)=Wm_tab(1)*(-prod1)
    omegadot(2)=Wm_tab(2)*(-1.5*prod1-0.5*prod2+0.5*prod3)
    omegadot(3)=Wm_tab(3)*(prod2-prod3)
    omegadot(4)=Wm_tab(4)*(2*prod1)
    omegadot(5)=Wm_tab(5)*(prod1-prod2+prod3)

  end subroutine Andersen

  subroutine OSK(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in)    :: temp 
    real(8), intent(out)   :: omegadot(ns)
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(ns), Tdiff
    real(8) :: prod1,prod2,prod3,prod4,prod5,prod6

    do is = 1, ns
      roi(is) = max(roi(is), 0.d0)
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
    enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1
    ! species: [O2, CH4, H2O, CO2]

    ! CH4 + 2 O2 => CO2 + 2 H2O
    prod1 = f_kf(1,Tint,Tdiff)*(coi(2)**0.7)*(coi(1)**0.8)

     
    ! Chemical Source Terms
    omegadot = 0d0
    omegadot(1)=Wm_tab(1)*(-2*prod1)
    omegadot(2)=Wm_tab(2)*(-prod1)
    omegadot(3)=Wm_tab(3)*(2*prod1)
    omegadot(4)=Wm_tab(4)*(prod1)

    !omegadot(5)=Wm_tab(5)*(prod1-prod2+prod3)

  end subroutine OSK
  


end module WD_mod
