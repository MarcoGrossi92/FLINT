module JLRs_mod
  implicit none
contains
  ! John Lindstedt with Recombination Reaction Mechanism
  ! 9 species & 7 reactions
  subroutine JLR(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in)    :: temp 
    real(8), intent(out)   :: omegadot(ns)
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(ns), Tdiff
    real(8) :: prod1,prod2,prod3,prod4,prod5,prod6,prod7

   do is = 1, ns
     ! no cut-off: zeroing coi below 1e-10 made the rates jump there (see fpow)
     coi(is)=max(roi(is),0d0)/Wm_tab(is)  ! kmol/m^3
   enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! 0.5 CH4 + 1.25 O2 --> CO + 2 H2 - 0.5 CH4 + 0.75 O2
    prod1 = f_kf(1,Tint,Tdiff)*fpow(coi(2),0.50d0)*(coi(1)**1.25)

    ! CH4 + H2O --> CO + 3 H2
    prod2 = f_kf(2,Tint,Tdiff)*(coi(2)*coi(3))

    ! CO + H2O <--> CO2 + H2
    prod3 = f_kf(3,Tint,Tdiff)*coi(4)*coi(3)- &
            f_kb(3,Tint,Tdiff)*coi(5)*coi(6)

    ! 1/4 H2 + 3/2 O2 <--> H2O + O2 - 3/4 H2
    ! reverse term regularised continuously: an on/off switch at [H2]=1e-10
    ! pins H2 at the switch and collapses the ODE step size
    prod4 = f_kf(4,Tint,Tdiff)*fpow(coi(6),0.25d0)*(coi(1)**1.50) - &
            f_kb(4,Tint,Tdiff)*(coi(3))*(coi(1))*(max(coi(6),1.d-10)**(-0.75))

    ! O2 <--> 2O
    prod5 = f_kf(5,Tint,Tdiff)*coi(1)-f_kb(5,Tint,Tdiff)*coi(8)**2

    ! H2O <--> H + OH
    prod6 = f_kf(6,Tint,Tdiff)*coi(3)-f_kb(6,Tint,Tdiff)*coi(7)*coi(9)

    ! OH + H2 <--> H + H2O
    prod7 = f_kf(7,Tint,Tdiff)*coi(9)*coi(6)- &
            f_kb(7,Tint,Tdiff)*coi(7)*coi(3)
     
    ! Chemical Source Terms
    omegadot = 0d0
    omegadot(1)=Wm_tab(1)*(-0.5*prod1-0.5*prod4-prod5)        !O2    [kg/(m3*s)]
    omegadot(2)=Wm_tab(2)*(-prod1-prod2)                      !CH4   [kg/(m3*s)]
    omegadot(3)=Wm_tab(3)*(-prod2-prod3+prod4-prod6+prod7)    !H2O   [kg/(m3*s)]
    omegadot(4)=Wm_tab(4)*(-prod3+prod1+prod2)                !CO    [kg/(m3*s)]
    omegadot(5)=Wm_tab(5)*(prod3)                             !CO2   [kg/(m3*s)]
    omegadot(6)=Wm_tab(6)*(2*prod1+3*prod2+prod3-prod4-prod7) !H2    [kg/(m3*s)]
    omegadot(7)=Wm_tab(7)*(prod6+prod7)                       !H     [kg/(m3*s)]
    omegadot(8)=Wm_tab(8)*(2*prod5)                           !O     [kg/(m3*s)]
    omegadot(9)=Wm_tab(9)*(prod6-prod7)                       !OH    [kg/(m3*s)]

  end subroutine JLR

  ! Frassoldati: John Lindstedt with Recombination Reaction Mechanism
  ! 9 species & 6 reactions
  subroutine Frassoldati(roi,temp,omegadot)
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
     ! no cut-off: zeroing coi below 1e-10 made the rates jump there (see fpow)
     coi(is)=max(roi(is),0d0)/Wm_tab(is)  ! kmol/m^3
   enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! 0.5 CH4 + 1.25 O2 --> CO + 2 H2 - 0.5 CH4 + 0.75 O2
    prod1 = f_kf(1,Tint,Tdiff)*fpow(coi(2),0.50d0)*(coi(1)**1.30)

    ! CH4 + H2O --> CO + 3 H2
    prod2 = f_kf(2,Tint,Tdiff)*(coi(2)*coi(3))

    ! CO + H2O <--> CO2 + H2
    prod3 = f_kf(3,Tint,Tdiff)*coi(4)*coi(3)- &
            f_kb(3,Tint,Tdiff)*coi(5)*coi(6)

    ! 1/4 H2 + 3/2 O2 <--> H2O + O2 - 3/4 H2
    ! reverse term regularised continuously: an on/off switch at [H2]=1e-10
    ! pins H2 at the switch and collapses the ODE step size
    prod4 = f_kf(4,Tint,Tdiff)*fpow(coi(6),0.30d0)*(coi(1)**1.55) - &
            f_kb(4,Tint,Tdiff)*(coi(3))*(coi(1))*(max(coi(6),1.d-10)**(-0.75))

    ! O2 <--> 2O
    prod5 = f_kf(5,Tint,Tdiff)*coi(1)-f_kb(5,Tint,Tdiff)*coi(8)**2

    ! H2O <--> H + OH
    prod6 = f_kf(6,Tint,Tdiff)*coi(3)-f_kb(6,Tint,Tdiff)*coi(7)*coi(9)
     
    ! Chemical Source Terms
    omegadot = 0d0
    omegadot(1)=Wm_tab(1)*(-0.5*prod1-0.5*prod4-prod5)        !O2    [kg/(m3*s)]
    omegadot(2)=Wm_tab(2)*(-prod1-prod2)                      !CH4   [kg/(m3*s)]
    omegadot(3)=Wm_tab(3)*(-prod2-prod3+prod4-prod6)          !H2O   [kg/(m3*s)]
    omegadot(4)=Wm_tab(4)*(-prod3+prod1+prod2)                !CO    [kg/(m3*s)]
    omegadot(5)=Wm_tab(5)*(prod3)                             !CO2   [kg/(m3*s)]
    omegadot(6)=Wm_tab(6)*(2*prod1+3*prod2+prod3-prod4)       !H2    [kg/(m3*s)]
    omegadot(7)=Wm_tab(7)*(prod6)                             !H     [kg/(m3*s)]
    omegadot(8)=Wm_tab(8)*(2*prod5)                           !O     [kg/(m3*s)]
    omegadot(9)=Wm_tab(9)*(prod6)                             !OH    [kg/(m3*s)]

  end subroutine Frassoldati

  !> John Lindstedt with Recombination Reaction Mechanism adapted for
  !> kerosene decomposition products combustion 10 species & 8 reactions
  !> SPECIES: O2, C2H4, H2O, CO, CO2, H2, H, O, OH, RP-1 surrogate (C12H24)
  subroutine CKJLR10sp(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    real(8), intent(inout)  :: roi(ns)
    real(8), intent(in)  :: temp 
    real(8), intent(out) :: omegadot(ns)
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(ns), Tdiff
    real(8) :: prod1, prod2, prod3, prod4, prod5, prod6, prod7, prod8

    do is = 1, ns
      ! no cut-off: zeroing coi below 1e-10 made the rates jump there (see fpow)
      coi(is)=max(roi(is),0d0)/Wm_tab(is)  ! kmol/m^3
    enddo

    ! Preliminary
    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    !> 0.5 C2H4 + 1.25 O2 --> 2 CO + 2 H2 - 0.5 C2H4 + 0.25 O2
    prod1 = f_kf(1,Tint,Tdiff)*fpow(coi(2),0.50d0)*(coi(1)**1.25)

    !> C2H4 + H2O --> 2 CO + 4 H2 - H2O
    prod2 = f_kf(2,Tint,Tdiff)*coi(2)*coi(3)

    !> CO + H2O <--> CO2 + H2
    prod3 = f_kf(3,Tint,Tdiff)*coi(4)*coi(3)-&
            f_kb(3,Tint,Tdiff)*coi(5)*coi(6)

    !> 1/4 H2 + 3/2 O2 <--> H2O + O2 - 3/4 H2
    ! reverse term regularised continuously: an on/off switch at [H2]=1e-10
    ! pins H2 at the switch and collapses the ODE step size
    prod4 = f_kf(4,Tint,Tdiff)*fpow(coi(6),0.25d0)*(coi(1)**1.50) - &
            f_kb(4,Tint,Tdiff)*coi(3)*coi(1)*(max(coi(6),1.d-10)**(-0.75))

    !> O2 <--> 2O
    prod5 = f_kf(5,Tint,Tdiff)*coi(1)-f_kb(5,Tint,Tdiff)*coi(8)**2

    !> H2O <--> H + OH
    prod6 = f_kf(6,Tint,Tdiff)*coi(3)-f_kb(6,Tint,Tdiff)*coi(7)*coi(9)

    !> OH + H2 <--> H + H2O
    prod7 = f_kf(7,Tint,Tdiff)*coi(9)*coi(6)- &
            f_kb(7,Tint,Tdiff)*coi(7)*coi(3)

    !> C12H24 --> 6C2H4 
    prod8 = f_kf(8,Tint,Tdiff)*coi(10)

    !> Chemical Source Terms
    omegadot = 0d0
    omegadot(1) =Wm_tab(1)*(-prod1-0.5*prod4-prod5)             !O2     [kg/(m3*s)]  
    omegadot(2) =Wm_tab(2)*(-prod1-prod2+6*prod8)               !C2H4   [kg/(m3*s)]  
    omegadot(3) =Wm_tab(3)*(-2*prod2-prod3+prod4-prod6+prod7)   !H2O    [kg/(m3*s)]  
    omegadot(4) =Wm_tab(4)*(2*prod1+2*prod2-prod3)              !CO     [kg/(m3*s)]  
    omegadot(5) =Wm_tab(5)*(prod3)                              !CO2    [kg/(m3*s)]  
    omegadot(6) =Wm_tab(6)*(2*prod1+4*prod2+prod3-prod4-prod7)  !H2     [kg/(m3*s)]  
    omegadot(7) =Wm_tab(7)*(prod6+prod7)                        !H      [kg/(m3*s)]  
    omegadot(8) =Wm_tab(8)*(2*prod5)                            !O      [kg/(m3*s)]  
    omegadot(9) =Wm_tab(9)*(prod6-prod7)                        !OH     [kg/(m3*s)]      
    omegadot(10)=Wm_tab(10)*(-prod8)                            !C12H24 [kg/(m3*s)]
    
  end subroutine CKJLR10sp


  pure function fpow(c,a) result(f)
    ! c**a for a fractional order 0<a<1, linear below c0. The plain power has an infinite slope at c=0,
    ! and the former cut-off (coi set to 0 below 1e-10) made the rates jump at c0; either one pins a
    ! trace species there and collapses the step size of the stiff ODE solvers (sdirk4b -6/-7)
    real(8), intent(in) :: c, a
    real(8) :: f
    real(8), parameter :: c0 = 1d-10 ! kmol/m^3
    if (c >= c0) then
      f = c**a
    else
      f = c0**a*(c/c0)
    end if
  end function fpow

end module JLRs_mod