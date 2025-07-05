
module U_Lib_ChemMech
  use U_Lib_Thermodynamic

  implicit none
  integer :: nrc
  real(8), dimension(:,:), allocatable, public :: ni1_tab
  real(8), dimension(:,:), allocatable, public :: ni2_tab
  real(8), dimension(:,:), allocatable, public :: kf_tab
  real(8), dimension(:,:), allocatable, public :: kb_tab
  real(8), dimension(:,:), allocatable, public :: epsch_tab

  !> Concrete procedure pointing to one of the subroutine realizations
  procedure(chemsource_if), pointer, public :: chemistry_source

  !> Abstract interface relative to the finite-rate reactions source procedure
  abstract interface
  subroutine chemsource_if(roi,temp,omegadot,rotot)
    use U_Lib_Thermodynamic
    implicit none
    integer :: is, T_i, Tint(2)
    real(8), intent(in) :: roi(nsc)
    real(8), intent(in) :: temp 
    real(8), intent(out) :: omegadot(nsc)
    real(8), intent(in) :: rotot
    real(8) :: coi(nsc+1), Tdiff
    !--------------------------------------------------------------
  end subroutine chemsource_if
  end interface


contains


  subroutine Assign_Mechanism(mad_world)
    implicit none
    character(*), intent(in) :: mad_world

    select case(mad_world)
    case('WD')
      chemistry_source => WD
    case('JLR-Nasuti')
      chemistry_source => JLR
    case('Frassoldati')
      chemistry_source => Frassoldati
    case('CoronettiC4H6')
      chemistry_source => Coronetti
    case('CKJLR-10sp')
      chemistry_source => CKJLR10sp
    case('Singh')
      chemistry_source => Singh
    case('Singh-WC32')
      chemistry_source => Singh_WC32
    case('Frolov')
      chemistry_source => Frolov
    case('Nassini')
      chemistry_source => Nassini_4
    case('ciottoli20')
      chemistry_source => ciottoli20
    case default
      chemistry_source => general
    end select

  end subroutine Assign_Mechanism


  ! General mechanism
  subroutine general(roi,temp,omegadot,rotot)
    use U_Lib_Thermodynamic, only: nsc
    implicit none
    real(8), intent(in) :: roi(nsc)
    real(8), intent(in) :: temp 
    real(8), intent(out) :: omegadot(nsc)
    real(8), intent(in) :: rotot
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(nsc+1), Tdiff
    real(8) :: prod1, prod2, deltani
    integer :: s, ir

    do is = 1, nsc
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
      ! Loop done in order to avoid numerical issues in omegadot evaluation
      ! Very low coi could produce finite prods
      if (coi(is).lt.1d-10) coi(is) = 0.0d0
    enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    do s = 1, nsc
      omegadot(s) = 0.d0
      do ir = 1, nrc
        deltani = ni2_tab(s,ir)-ni1_tab(s,ir)
        if (deltani==0) cycle
        prod1 = 1d0
        prod2 = 1d0
        coi(nsc+1) = 0d0
        do is = 1, nsc+1
          coi(nsc+1) = coi(nsc+1)+coi(is)*epsch_tab(is,ir)  !terzo corpo
          if(ni1_tab(is,ir).ne.0) prod1 = prod1*coi(is)**ni1_tab(is,ir)
          if(ni2_tab(is,ir).ne.0) prod2 = prod2*coi(is)**ni2_tab(is,ir)
        enddo  !is=1,nsc+1
        prod1 = prod1*comp_ch_tabT(ir,kf_tab,Tint,Tdiff)
        prod2 = prod2*comp_ch_tabT(ir,kb_tab,Tint,Tdiff)
        omegadot(s) = omegadot(s)+deltani*(prod1-prod2)
      enddo  !ir=1,nrc
      omegadot(s) = Wm_tab(s)*omegadot(s)  !kg/(m3*s)
    enddo  !s=1,nsc

  end subroutine general


  ! WD: Global Westbrook-Dryer mechanism
  ! 5 species & 3 reactions
  subroutine WD(roi,temp,omegadot,rotot)
    use U_Lib_Thermodynamic, only: nsc
    implicit none
    real(8), intent(in) :: roi(nsc)
    real(8), intent(in)    :: temp 
    real(8), intent(out)   :: omegadot(nsc)
    real(8), intent(in)    :: rotot
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(nsc+1), Tdiff
    real(8) :: prod1,prod2,prod3,prod4,prod5,prod6

   do is = 1, nsc
     coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
     ! Loop done in order to avoid numerical issues in omegadot evaluation
     ! Very low coi could produce finite prods
     if (coi(is).lt.1d-10) coi(is) = 0.0d0
   enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! species: [CH4, O2, CO2, H2O, CO]

    ! CH4 + 1.5 O2 => CO + 2 H2O
    prod1 = comp_ch_tabT(1,kf_tab,Tint,Tdiff)*(coi(1)**0.70)*(coi(2)**0.80)

    ! CO + 0.5 O2 + H2O => CO2 + H2O
    prod2 = comp_ch_tabT(2,kf_tab,Tint,Tdiff)*(coi(4)*coi(5)*coi(2)**0.50)

    ! CO2 => CO + 0.5 O2
    prod3 = comp_ch_tabT(3,kf_tab,Tint,Tdiff)*coi(3)
     
    ! Chemical Source Terms
    omegadot = 0d0
    omegadot(1)=Wm_tab(1)*(-prod1)
    omegadot(2)=Wm_tab(2)*(-1.5*prod1-0.5*prod2+0.5*prod3)
    omegadot(3)=Wm_tab(3)*(prod2-prod3)
    omegadot(4)=Wm_tab(4)*(2*prod1)
    omegadot(5)=Wm_tab(5)*(prod1-prod2+prod3)

  end subroutine WD


  ! John Lindstedt with Recombination Reaction Mechanism
  ! 9 species & 7 reactions
  subroutine JLR(roi,temp,omegadot,rotot)
    use U_Lib_Thermodynamic, only: nsc
    implicit none
    real(8), intent(in) :: roi(nsc)
    real(8), intent(in)    :: temp 
    real(8), intent(out)   :: omegadot(nsc)
    real(8), intent(in)    :: rotot
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(nsc+1), Tdiff
    real(8) :: prod1,prod2,prod3,prod4,prod5,prod6,prod7

   do is = 1, nsc
     coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
     ! Loop done in order to avoid numerical issues in omegadot evaluation
     ! Very low coi could produce finite prods
     if (coi(is).lt.1d-10) coi(is) = 0.0d0
   enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! 0.5 CH4 + 1.25 O2 --> CO + 2 H2 - 0.5 CH4 + 0.75 O2
    prod1 = comp_ch_tabT(1,kf_tab,Tint,Tdiff)*(coi(2)**0.50)*(coi(1)**1.25)

    ! CH4 + H2O --> CO + 3 H2
    prod2 = comp_ch_tabT(2,kf_tab,Tint,Tdiff)*(coi(2)*coi(3))

    ! CO + H2O <--> CO2 + H2
    prod3 = comp_ch_tabT(3,kf_tab,Tint,Tdiff)*coi(4)*coi(3)- &
            comp_ch_tabT(3,kb_tab,Tint,Tdiff)*coi(5)*coi(6)

    ! 1/4 H2 + 3/2 O2 <--> H2O + O2 - 3/4 H2
    if (coi(6) < 1.d-10) then
      prod4 = (comp_ch_tabT(4,kf_tab,Tint,Tdiff)*(coi(6)**0.25)*(coi(1)**1.50))
    else
      prod4 = comp_ch_tabT(4,kf_tab,Tint,Tdiff)*(coi(6)**0.25)* & 
              (coi(1)**1.50)-comp_ch_tabT(4,kb_tab,Tint,Tdiff)*(coi(3))*(coi(1))*(coi(6)**(-0.75))
    endif

    ! O2 <--> 2O
    prod5 = comp_ch_tabT(5,kf_tab,Tint,Tdiff)*coi(1)-comp_ch_tabT(5,kb_tab,Tint,Tdiff)*coi(8)**2

    ! H2O <--> H + OH
    prod6 = comp_ch_tabT(6,kf_tab,Tint,Tdiff)*coi(3)-comp_ch_tabT(6,kb_tab,Tint,Tdiff)*coi(7)*coi(9)

    ! OH + H2 <--> H + H2O
    prod7 = comp_ch_tabT(7,kf_tab,Tint,Tdiff)*coi(9)*coi(6)- &
            comp_ch_tabT(7,kb_tab,Tint,Tdiff)*coi(7)*coi(3)
     
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
  subroutine Frassoldati(roi,temp,omegadot,rotot)
    use U_Lib_Thermodynamic, only: nsc
    implicit none
    real(8), intent(in) :: roi(nsc)
    real(8), intent(in)    :: temp 
    real(8), intent(out)   :: omegadot(nsc)
    real(8), intent(in)    :: rotot
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(nsc+1), Tdiff
    real(8) :: prod1,prod2,prod3,prod4,prod5,prod6

   do is = 1, nsc
     coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
     ! Loop done in order to avoid numerical issues in omegadot evaluation
     ! Very low coi could produce finite prods
     if (coi(is).lt.1d-10) coi(is) = 0.0d0
   enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! 0.5 CH4 + 1.25 O2 --> CO + 2 H2 - 0.5 CH4 + 0.75 O2
    prod1 = comp_ch_tabT(1,kf_tab,Tint,Tdiff)*(coi(2)**0.50)*(coi(1)**1.30)

    ! CH4 + H2O --> CO + 3 H2
    prod2 = comp_ch_tabT(2,kf_tab,Tint,Tdiff)*(coi(2)*coi(3))

    ! CO + H2O <--> CO2 + H2
    prod3 = comp_ch_tabT(3,kf_tab,Tint,Tdiff)*coi(4)*coi(3)- &
            comp_ch_tabT(3,kb_tab,Tint,Tdiff)*coi(5)*coi(6)

    ! 1/4 H2 + 3/2 O2 <--> H2O + O2 - 3/4 H2
    if (coi(6) < 1.d-10) then
      prod4 = (comp_ch_tabT(4,kf_tab,Tint,Tdiff)*(coi(6)**0.3)*(coi(1)**1.55))
    else
      prod4 = comp_ch_tabT(4,kf_tab,Tint,Tdiff)*(coi(6)**0.3)* & 
              (coi(1)**1.55)-comp_ch_tabT(4,kb_tab,Tint,Tdiff)*(coi(3))*(coi(1))*(coi(6)**(-0.75))
    endif

    ! O2 <--> 2O
    prod5 = comp_ch_tabT(5,kf_tab,Tint,Tdiff)*coi(1)-comp_ch_tabT(5,kb_tab,Tint,Tdiff)*coi(8)**2

    ! H2O <--> H + OH
    prod6 = comp_ch_tabT(6,kf_tab,Tint,Tdiff)*coi(3)-comp_ch_tabT(6,kb_tab,Tint,Tdiff)*coi(7)*coi(9)
     
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

  ! Coronetti for butadiene combustion (DOI:10.2514/1.B34760)
  ! Order of species: O2, C4H6, H2O, CO, CO2, H2, O, H, OH
  subroutine Coronetti(roi,temp,omegadot,rotot)
    use U_Lib_Thermodynamic, only: nsc
    implicit none
    real(8), intent(in) :: roi(nsc)
    real(8), intent(in) :: temp 
    real(8), intent(out) :: omegadot(nsc)
    real(8), intent(in) :: rotot
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(nsc+1), Tdiff
    real(8) :: prod1, prod2, prod3, prod4, prod5, prod6
    real(8) :: kf, kb, CP1, CP2, xi, lin1, lin2
    real(8), parameter :: limitH2 = 1d-10, limitO2 = 1d-10, sigma = 23d0, tau = 17d0

    do is = 1, nsc
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
    enddo
    
    ! Preliminary
    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! irreversible reaction: C4H6 + 2 O2 -> 4 CO + 3 H2 
    kf = comp_ch_tabT(1,kf_tab,Tint,Tdiff)
    CP1 = (coi(2)**0.5)*(coi(1)**1.25)
    if (( coi(1) < limitO2 ).or.(coi(2)<limitO2) )then  
      prod1 = 0.d0
    else
      prod1 = kf*CP1
    end if

    ! irreversible reaction: C4H6 + 4 H2O -> 4 CO + 7 H2
    prod2 = comp_ch_tabT(2,kf_tab,Tint,Tdiff)*coi(2)*coi(3)
      
    ! reversible reaction: CO + H2O <--> CO2 + H2
    kf = comp_ch_tabT(3,kf_tab,Tint,Tdiff)
    kb = comp_ch_tabT(3,kb_tab,Tint,Tdiff)
    CP1 = coi(4)*coi(3)
    CP2 = coi(5)*coi(6)
    prod3 = kf*CP1 - kb*CP2

    ! reversible reaction: H2 + 1/2 O2 <--> H2O
    kf = comp_ch_tabT(4,kf_tab,Tint,Tdiff)
    kb = comp_ch_tabT(4,kb_tab,Tint,Tdiff)
    CP1 = (coi(6)**0.25)*(coi(1)**1.50)
    CP2 = (coi(3))*(coi(1))*(coi(6)**(-0.75))
    if (( coi(6) < limitH2 ).or.(coi(1)< limitH2)) then
      CP1 = 0.d0 
    end if
    if ((coi(3)<limitH2).or.(coi(1)<limitH2).or.(coi(6)<limitH2)) then
      CP2 = 0.d0 
    endif
    prod4 = kf*CP1 - kb*CP2
   
    ! reversible reaction: O2 <--> 2 O
    prod5 = comp_ch_tabT(5,kf_tab,Tint,Tdiff)*coi(1) &      
          - comp_ch_tabT(5,kb_tab,Tint,Tdiff)*coi(7)**2
      
    ! reversible reaction: H2O <--> OH + H
    prod6 = comp_ch_tabT(6,kf_tab,Tint,Tdiff)*coi(3) &
          - comp_ch_tabT(6,kb_tab,Tint,Tdiff)*coi(9)*coi(8)
      
    omegadot(1)=Wm_tab(1)*(-2.d0*prod1-0.5d0*prod4-prod5)         ! O2 
    omegadot(2)=Wm_tab(2)*(-prod1-prod2)                          ! C4H6
    omegadot(3)=Wm_tab(3)*(-4.d0*prod2-prod3+prod4-prod6)         ! H2O 
    omegadot(4)=Wm_tab(4)*(4.d0*prod1+4.d0*prod2-prod3)           ! CO
    omegadot(5)=Wm_tab(5)*(prod3)                                 ! CO2
    omegadot(6)=Wm_tab(6)*(3.d0*prod1+7.d0*prod2+prod3-prod4)     ! H2
    omegadot(7)=Wm_tab(7)*(2.d0*prod5)                            ! O
    omegadot(8)=Wm_tab(8)*(prod6)                                 ! H 
    omegadot(9)=Wm_tab(9)*(prod6)                                 ! OH

  end subroutine Coronetti

  !> John Lindstedt with Recombination Reaction Mechanism adapted for
  !> kerosene decomposition products combustion 10 species & 8 reactions
  !> SPECIES: O2, C2H4, H2O, CO, CO2, H2, H, O, OH, RP-1 (C12H24)
  subroutine CKJLR10sp(roi,temp,omegadot,rotot)
    use U_Lib_Thermodynamic, only: nsc
    implicit none
    real(8), intent(in)  :: roi(nsc)
    real(8), intent(in)  :: temp 
    real(8), intent(in)  :: rotot
    real(8), intent(out) :: omegadot(nsc)
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(nsc+1), Tdiff
    real(8) :: prod1, prod2, prod3, prod4, prod5, prod6, prod7, prod8
    real(8), parameter :: limitH2=1.d-10

    do is = 1, nsc
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
      ! Loop done in order to avoid numerical issues in omegadot evaluation
      ! Very low coi could produce finite prods
      if (coi(is).lt.1d-20) coi(is) = 0.0d0
    enddo

    ! Preliminary
    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    !> 0.5 C2H4 + 1.25 O2 --> 2 CO + 2 H2 - 0.5 C2H4 + 0.25 O2
    prod1=comp_ch_tabT(1,kf_tab,Tint,Tdiff)*(coi(2)**0.50)*(coi(1)**1.25)

    !> C2H4 + H2O --> 2 CO + 4 H2 - H2O
    prod2=comp_ch_tabT(2,kf_tab,Tint,Tdiff)*coi(2)*coi(3)

    !> CO + H2O <--> CO2 + H2
    prod3=comp_ch_tabT(3,kf_tab,Tint,Tdiff)*coi(4)*coi(3)-comp_ch_tabT(3,kb_tab,Tint,Tdiff)*coi(5)*coi(6)

    !> 1/4 H2 + 3/2 O2 <--> H2O + O2 - 3/4 H2
    if (coi(6) < limitH2) then
      prod4=comp_ch_tabT(4,kf_tab,Tint,Tdiff)*(coi(6)**0.25)*(coi(1)**1.50)
    else
      prod4=comp_ch_tabT(4,kf_tab,Tint,Tdiff)*(coi(6)**0.25)*(coi(1)**1.50) - &
            comp_ch_tabT(4,kb_tab,Tint,Tdiff)*coi(3)*coi(1)*(coi(6)**(-0.75))
    endif

    !> O2 <--> 2O
    prod5=comp_ch_tabT(5,kf_tab,Tint,Tdiff)*coi(1)-comp_ch_tabT(5,kb_tab,Tint,Tdiff)*coi(8)**2

    !> H2O <--> H + OH
    prod6=comp_ch_tabT(6,kf_tab,Tint,Tdiff)*coi(3)-comp_ch_tabT(6,kb_tab,Tint,Tdiff)*coi(7)*coi(9)

    !> OH + H2 <--> H + H2O
    prod7=comp_ch_tabT(7,kf_tab,Tint,Tdiff)*coi(9)*coi(6)-comp_ch_tabT(7,kb_tab,Tint,Tdiff)*coi(7)*coi(3)

    !> C12H24 --> 6C2H4 
    prod8=comp_ch_tabT(8,kf_tab,Tint,Tdiff)*coi(10)

    !> Chemical Source Terms
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

  ! Singh (1994) global mechanism for ethylene combustion.
  ! 10 reversible reactions for 9 species.
  subroutine Singh(roi,temp,omegadot,rotot)
    use U_Lib_Thermodynamic, only: nsc
    implicit none
    real(8), intent(in) :: roi(nsc)
    real(8), intent(in) :: temp 
    real(8), intent(out) :: omegadot(nsc)
    real(8), intent(in) :: rotot
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(nsc+1), Tdiff
    real(8) :: M, prod1f, prod1b, prod1, prod2f, prod2b, prod2
    real(8) :: prod3f, prod3b, prod3, prod4f, prod4b, prod4
    real(8) :: prod5f, prod5b, prod5, prod6f, prod6b, prod6
    real(8) :: prod7f, prod7b, prod7, prod8f, prod8b, prod8
    real(8) :: prod9f, prod9b, prod9, prod10f, prod10b, prod10

    do is = 1, nsc
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
      if (coi(is).lt.1d-12) coi(is) = 0d0
    enddo
    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! 3rd body eff: 2.5 H2, 16.0 H2O, 1.0 else
    M = sum(coi(1:5))+2.5d0*coi(6)+16d0*coi(7)+coi(8)+coi(9)
    ! C2H4 + O2 <-> 2CO+ 2H2
    prod1f=comp_ch_tabT(1,kf_tab,Tint,Tdiff)*coi(2)*coi(8)
    prod1b=comp_ch_tabT(1,kb_tab,Tint,Tdiff)*coi(4)**2*coi(6)**2
    prod1=prod1f-prod1b
    ! CO + O <-> CO2 + M
    prod2f=comp_ch_tabT(2,kf_tab,Tint,Tdiff)*coi(4)*coi(9)*M
    prod2b=comp_ch_tabT(2,kb_tab,Tint,Tdiff)*coi(5)*M
    prod2=prod2f-prod2b
    !CO + OH <-> CO2 + H
    prod3f=comp_ch_tabT(3,kf_tab,Tint,Tdiff)*coi(4)*coi(3)
    prod3b=comp_ch_tabT(3,kb_tab,Tint,Tdiff)*coi(5)*coi(1)
    prod3=prod3f-prod3b
    !H2 + O2 <-> OH+OH
    prod4f=comp_ch_tabT(4,kf_tab,Tint,Tdiff)*coi(6)*coi(8)
    prod4b=comp_ch_tabT(4,kb_tab,Tint,Tdiff)*coi(3)*coi(3)
    prod4=prod4f-prod4b
    !H + O2 <-> OH + O
    prod5f=comp_ch_tabT(5,kf_tab,Tint,Tdiff)*coi(1)*coi(8)
    prod5b=comp_ch_tabT(5,kb_tab,Tint,Tdiff)*coi(3)*coi(9)
    prod5=prod5f-prod5b
    !OH +H2 <-> H2O + H
    prod6f=comp_ch_tabT(6,kf_tab,Tint,Tdiff)*coi(3)*coi(6)
    prod6b=comp_ch_tabT(6,kb_tab,Tint,Tdiff)*coi(7)*coi(1)
    prod6=prod6f-prod6b
    !O + H2 <-> OH + H
    prod7f=comp_ch_tabT(7,kf_tab,Tint,Tdiff)*coi(9)*coi(6)
    prod7b=comp_ch_tabT(7,kb_tab,Tint,Tdiff)*coi(3)*coi(1)
    prod7=prod7f-prod7b
    !OH + OH <-> H2O + O
    prod8f=comp_ch_tabT(8,kf_tab,Tint,Tdiff)*coi(3)*coi(3)
    prod8b=comp_ch_tabT(8,kb_tab,Tint,Tdiff)*coi(7)*coi(9)
    prod8=prod8f-prod8b
    !H + H <-> H2 + M
    prod9f=comp_ch_tabT(9,kf_tab,Tint,Tdiff)*coi(1)*coi(1)*M
    prod9b=comp_ch_tabT(9,kb_tab,Tint,Tdiff)*coi(6)*M
    prod9=prod9f-prod9b
    !H + OH <-> H2O + M
    prod10f=comp_ch_tabT(10,kf_tab,Tint,Tdiff)*coi(1)*coi(3)*M
    prod10b=comp_ch_tabT(10,kb_tab,Tint,Tdiff)*coi(7)*M
    prod10=prod10f-prod10b

    ! H
    omegadot(1)=wm_tab(1)*(prod3+prod6+prod7-2*prod9-prod10)
    ! C2H4
    omegadot(2)=wm_tab(2)*(-prod1)
    ! OH
    omegadot(3)=wm_tab(3)*(-prod3+2*prod4+prod5-prod6+ &
                            prod7-2*prod8-prod10)
    ! CO
    omegadot(4)=wm_tab(4)*(2*prod1-prod2-prod3)
    ! CO2
    omegadot(5)=wm_tab(5)*(prod2+prod3)
    ! H2
    omegadot(6)=wm_tab(6)*(2*prod1-prod4-prod6-prod7+prod9)
    ! H2O
    omegadot(7)=wm_tab(7)*(prod6+prod8+prod10)
    ! O2
    omegadot(8)=wm_tab(8)*(-prod1-prod4-prod5)
    ! O
    omegadot(9)=wm_tab(9)*(-prod2+prod5-prod7+prod8)

  end subroutine Singh


  ! Singh (1994) + paraffin cracking
  subroutine Singh_WC32(roi,temp,omegadot,rotot)
    use U_Lib_Thermodynamic, only: nsc
    implicit none
    real(8), intent(in) :: roi(nsc)
    real(8), intent(in) :: temp 
    real(8), intent(out) :: omegadot(nsc)
    real(8), intent(in) :: rotot
    ! Local
    integer :: is, T_i, Tint(2)
    integer, parameter :: iH=1, iC2H4=2, iOH=3, iCO=4, iCO2=5
    integer, parameter :: iH2=6, iH2O=7, iO2=8, iO=9, iC32H66=10
    real(8) :: coi(nsc+1), Tdiff
    real(8) :: M, prodf(11), prodb(11), prod(11)

    do is = 1, nsc
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
      if (coi(is).lt.1d-12) coi(is) = 0d0
    enddo
    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! 3rd body eff: 2.5 H2, 16.0 H2O, 1.0 else
    M = sum(coi(iH:iCO2)) + 2.5d0*coi(iH2) + 16d0*coi(iH2O) + coi(iO2) + coi(iO)
    ! C2H4 + O2 <-> 2CO+ 2H2
    prodf(1)=comp_ch_tabT(1,kf_tab,Tint,Tdiff)*coi(iC2H4)*coi(iO2)
    prodb(1)=comp_ch_tabT(1,kb_tab,Tint,Tdiff)*coi(iCO)**2*coi(iH2)**2
    prod(1)=prodf(1)-prodb(1)
    ! CO + O + M <-> CO2 + M
    prodf(2)=comp_ch_tabT(2,kf_tab,Tint,Tdiff)*coi(iCO)*coi(iO)*M
    prodb(2)=comp_ch_tabT(2,kb_tab,Tint,Tdiff)*coi(iCO2)*M
    prod(2)=prodf(2)-prodb(2)
    !CO + OH <-> CO2 + H
    prodf(3)=comp_ch_tabT(3,kf_tab,Tint,Tdiff)*coi(iCO)*coi(iOH)
    prodb(3)=comp_ch_tabT(3,kb_tab,Tint,Tdiff)*coi(iCO2)*coi(iH)
    prod(3)=prodf(3)-prodb(3)
    !H2 + O2 <-> OH+OH
    prodf(4)=comp_ch_tabT(4,kf_tab,Tint,Tdiff)*coi(iH2)*coi(iO2)
    prodb(4)=comp_ch_tabT(4,kb_tab,Tint,Tdiff)*coi(iOH)**2
    prod(4)=prodf(4)-prodb(4)
    !H + O2 <-> OH + O
    prodf(5)=comp_ch_tabT(5,kf_tab,Tint,Tdiff)*coi(iH)*coi(iO2)
    prodb(5)=comp_ch_tabT(5,kb_tab,Tint,Tdiff)*coi(iOH)*coi(iO)
    prod(5)=prodf(5)-prodb(5)
    !OH +H2 <-> H2O + H
    prodf(6)=comp_ch_tabT(6,kf_tab,Tint,Tdiff)*coi(iOH)*coi(iH2)
    prodb(6)=comp_ch_tabT(6,kb_tab,Tint,Tdiff)*coi(iH2O)*coi(iH)
    prod(6)=prodf(6)-prodb(6)
    !O + H2 <-> OH + H
    prodf(7)=comp_ch_tabT(7,kf_tab,Tint,Tdiff)*coi(iO)*coi(iH2)
    prodb(7)=comp_ch_tabT(7,kb_tab,Tint,Tdiff)*coi(iOH)*coi(iH)
    prod(7)=prodf(7)-prodb(7)
    !OH + OH <-> H2O + O
    prodf(8)=comp_ch_tabT(8,kf_tab,Tint,Tdiff)*coi(iOH)**2
    prodb(8)=comp_ch_tabT(8,kb_tab,Tint,Tdiff)*coi(iH2O)*coi(iO)
    prod(8)=prodf(8)-prodb(8)
    !H + H + M <-> H2 + M
    prodf(9)=comp_ch_tabT(9,kf_tab,Tint,Tdiff)*coi(iH)**2*M
    prodb(9)=comp_ch_tabT(9,kb_tab,Tint,Tdiff)*coi(iH2)*M
    prod(9)=prodf(9)-prodb(9)
    !H + OH <-> H2O + M
    prodf(10)=comp_ch_tabT(10,kf_tab,Tint,Tdiff)*coi(iH)*coi(iOH)*M
    prodb(10)=comp_ch_tabT(10,kb_tab,Tint,Tdiff)*coi(iH2O)*M
    prod(10)=prodf(10)-prodb(10)
    !C32H66 --> 16C2H4 + H2
    if (coi(iC32H66)>1d-10) then
      prodf(11)=comp_ch_tabT(11,kf_tab,Tint,Tdiff)*(coi(iC32H66))
    else
      prodf(11)=0d0
    endif
    prodb(11)=0d0
    prod(11)=prodf(11)

    ! H
    omegadot(iH)=wm_tab(iH)*(prod(3)-prod(5)+prod(6)+prod(7)-2*prod(9)-prod(10))
    ! C2H4
    omegadot(iC2H4)=wm_tab(iC2H4)*(-prod(1)+16*prod(11))
    ! OH
    omegadot(iOH)=wm_tab(iOH)*(-prod(3)+2*prod(4)+prod(5)-prod(6)+ &
                                prod(7)-2*prod(8)-prod(10))
    ! CO
    omegadot(iCO)=wm_tab(iCO)*(2*prod(1)-prod(2)-prod(3))
    ! CO2
    omegadot(iCO2)=wm_tab(iCO2)*(prod(2)+prod(3))
    ! H2
    omegadot(iH2)=wm_tab(iH2)*(2*prod(1)-prod(4)-prod(6)-prod(7)+prod(9)+prod(11))
    ! H2O
    omegadot(iH2O)=wm_tab(iH2O)*(prod(6)+prod(8)+prod(10))
    ! O2
    omegadot(iO2)=wm_tab(iO2)*(-prod(1)-prod(4)-prod(5))
    ! O
    omegadot(iO)=wm_tab(iO)*(-prod(2)+prod(5)-prod(7)+prod(8))
    ! C32H66
    omegadot(iC32H66)=wm_tab(iC32H66)*(-prod(11))

  end subroutine Singh_WC32


  ! Frolov: hydrogen/air with 3 species, one reaction 
  subroutine Frolov(roi,temp,omegadot,rotot)
    use U_Lib_Thermodynamic, only: nsc
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

!    do is = 1, nsc
!      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
!      if (coi(is).lt.1d-12) coi(is) = 0d0
!    enddo

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
    omegadot(3)= Wm_tab(3)*(-1.d0 * prod1)            !H2    [kg/(m3*s)]

  end subroutine Nassini_4

  subroutine ciottoli20(roi,temp,omegadot,rotot)
    use U_Lib_Thermodynamic, only: nsc
    implicit none
    integer :: is, T_i, Tint(2)
    real(8), intent(in)    :: roi(nsc)
    real(8), intent(in)    :: temp 
    real(8), intent(out)   :: omegadot(nsc)
    real(8), intent(in)    :: rotot
    ! Local
    real(8) :: coi(nsc+1), Tdiff
    real(8) :: p, M
    real(8) :: prodf(1:104),prodb(1:104)

    !--------------------------------------------------------------

    do is = 1, nsc
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
      if (coi(is).lt.1d-10) coi(is) = 0d0
    enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    p = sum(roi) * sum(coi * Ri_tab) * temp

    ! reac n. 1: HCO + OH => CO + H2O
  prodf(1)=comp_ch_tabT(1,kf_tab,Tint,Tdiff)*(coi(7)**1.0)*(coi(1)**1.0)
  prodb(1)=comp_ch_tabT(1,kb_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(9)**1.0)
  ! reac n. 2: CO + H2O => HCO + OH
  prodf(2)=comp_ch_tabT(2,kf_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(9)**1.0)
  prodb(2)=comp_ch_tabT(2,kb_tab,Tint,Tdiff)*(coi(7)**1.0)*(coi(1)**1.0)
  ! reac n. 3: CO + OH => CO2 + H
  prodf(3)=comp_ch_tabT(3,kf_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(1)**1.0)
  prodb(3)=comp_ch_tabT(3,kb_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(3)**1.0)
  ! reac n. 4: CO2 + H => CO + OH
  prodf(4)=comp_ch_tabT(4,kf_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(3)**1.0)
  prodb(4)=comp_ch_tabT(4,kb_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(1)**1.0)
  ! reac n. 5: H + O2 => O + OH
  prodf(5)=comp_ch_tabT(5,kf_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(4)**1.0)
  prodb(5)=comp_ch_tabT(5,kb_tab,Tint,Tdiff)*(coi(5)**1.0)*(coi(1)**1.0)
  ! reac n. 6: O + OH => H + O2
  prodf(6)=comp_ch_tabT(6,kf_tab,Tint,Tdiff)*(coi(5)**1.0)*(coi(1)**1.0)
  prodb(6)=comp_ch_tabT(6,kb_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(4)**1.0)
  ! reac n. 7: H2 + O => H + OH
  prodf(7)=comp_ch_tabT(7,kf_tab,Tint,Tdiff)*(coi(6)**1.0)*(coi(5)**1.0)
  prodb(7)=comp_ch_tabT(7,kb_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(1)**1.0)
  ! reac n. 8: H + OH => H2 + O
  prodf(8)=comp_ch_tabT(8,kf_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(1)**1.0)
  prodb(8)=comp_ch_tabT(8,kb_tab,Tint,Tdiff)*(coi(6)**1.0)*(coi(5)**1.0)
  ! reac n. 9: H2O + O => 2 OH
  prodf(9)=comp_ch_tabT(9,kf_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(5)**1.0)
  prodb(9)=comp_ch_tabT(9,kb_tab,Tint,Tdiff)*(coi(1)**2.0)
  ! reac n. 10: 2 OH => H2O + O
  prodf(10)=comp_ch_tabT(10,kf_tab,Tint,Tdiff)*(coi(1)**2.0)
  prodb(10)=comp_ch_tabT(10,kb_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(5)**1.0)
  ! reac n. 11: H2 + OH => H + H2O
  prodf(11)=comp_ch_tabT(11,kf_tab,Tint,Tdiff)*(coi(6)**1.0)*(coi(1)**1.0)
  prodb(11)=comp_ch_tabT(11,kb_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(9)**1.0)
  ! reac n. 12: H + H2O => H2 + OH
  prodf(12)=comp_ch_tabT(12,kf_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(9)**1.0)
  prodb(12)=comp_ch_tabT(12,kb_tab,Tint,Tdiff)*(coi(6)**1.0)*(coi(1)**1.0)
  ! reac n. 13: HCO + M => CO + H + M
  M=0.d0+coi(1)*1.0+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*2.5+coi(7)*1.0+coi(8)*1.9+coi(9)*12.0+coi(10)*3.8+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.0+coi(16)*1.0+coi(17)*1.0+coi(18)*1.0+coi(19)*1.0+coi(20)*1.0
  prodf(13)=comp_ch_tabT(13,kf_tab,Tint,Tdiff)*(coi(7)**1.0)*M
  prodb(13)=comp_ch_tabT(13,kb_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(3)**1.0)*M
  ! reac n. 14: CO + H + M => HCO + M
  M=0.d0+coi(1)*1.0+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*2.5+coi(7)*1.0+coi(8)*1.9+coi(9)*12.0+coi(10)*3.8+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.0+coi(16)*1.0+coi(17)*1.0+coi(18)*1.0+coi(19)*1.0+coi(20)*1.0
  prodf(14)=comp_ch_tabT(14,kf_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(3)**1.0)*M
  prodb(14)=comp_ch_tabT(14,kb_tab,Tint,Tdiff)*(coi(7)**1.0)*M
  ! reac n. 15: CO + HO2 => CO2 + OH
  prodf(15)=comp_ch_tabT(15,kf_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(11)**1.0)
  prodb(15)=comp_ch_tabT(15,kb_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(1)**1.0)
  ! reac n. 16: CO2 + OH => CO + HO2
  prodf(16)=comp_ch_tabT(16,kf_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(1)**1.0)
  prodb(16)=comp_ch_tabT(16,kb_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(11)**1.0)
  ! reac n. 17: H2O + M => H + OH + M
  M=0.d0+coi(1)*1.0+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*2.5+coi(7)*1.0+coi(8)*1.9+coi(9)*12.0+coi(10)*3.8+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.0+coi(16)*1.0+coi(17)*1.0+coi(18)*1.0+coi(19)*1.0+coi(20)*1.0
  prodf(17)=comp_ch_tabT(17,kf_tab,Tint,Tdiff)*(coi(9)**1.0)*M
  prodb(17)=comp_ch_tabT(17,kb_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(1)**1.0)*M
  ! reac n. 18: H + OH + M => H2O + M
  M=0.d0+coi(1)*1.0+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*2.5+coi(7)*1.0+coi(8)*1.9+coi(9)*12.0+coi(10)*3.8+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.0+coi(16)*1.0+coi(17)*1.0+coi(18)*1.0+coi(19)*1.0+coi(20)*1.0
  prodf(18)=comp_ch_tabT(18,kf_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(1)**1.0)*M
  prodb(18)=comp_ch_tabT(18,kb_tab,Tint,Tdiff)*(coi(9)**1.0)*M
  ! reac n. 19: H + O2 (+M) <=> HO2 (+M)
  M=0.d0+coi(1)*1.0+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*2.5+coi(7)*1.0+coi(8)*1.9+coi(9)*12.0+coi(10)*3.8+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.0+coi(16)*1.0+coi(17)*1.0+coi(18)*1.0+coi(19)*1.0+coi(20)*1.0
  prodf(19)=comp_ch_tabT(19,kf_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(4)**1.0)*M
  prodb(19)=comp_ch_tabT(19,kb_tab,Tint,Tdiff)*(coi(11)**1.0)*M
  ! reac n. 20: CO + O (+M) <=> CO2 (+M)
  M=0.d0+coi(1)*1.0+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*2.5+coi(7)*1.0+coi(8)*1.9+coi(9)*12.0+coi(10)*3.8+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.0+coi(16)*1.0+coi(17)*1.0+coi(18)*1.0+coi(19)*1.0+coi(20)*1.0
  prodf(20)=comp_ch_tabT(20,kf_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(5)**1.0)*M 
  prodb(20)=comp_ch_tabT(20,kb_tab,Tint,Tdiff)*(coi(10)**1.0)*M
  ! reac n. 21: CO + O2 => CO2 + O
  prodf(21)=comp_ch_tabT(21,kf_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(4)**1.0)
  prodb(21)=comp_ch_tabT(21,kb_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(5)**1.0)
  ! reac n. 22: CO2 + O => CO + O2
  prodf(22)=comp_ch_tabT(22,kf_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(5)**1.0)
  prodb(22)=comp_ch_tabT(22,kb_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(4)**1.0)
  ! reac n. 23: H + HCO => CO + H2
  prodf(23)=comp_ch_tabT(23,kf_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(7)**1.0)
  prodb(23)=comp_ch_tabT(23,kb_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(6)**1.0)
  ! reac n. 24: CO + H2 => H + HCO
  prodf(24)=comp_ch_tabT(24,kf_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(6)**1.0)
  prodb(24)=comp_ch_tabT(24,kb_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(7)**1.0)
  ! reac n. 25: HCO + O => CO + OH
  prodf(25)=comp_ch_tabT(25,kf_tab,Tint,Tdiff)*(coi(7)**1.0)*(coi(5)**1.0)
  prodb(25)=comp_ch_tabT(25,kb_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(1)**1.0)
  ! reac n. 26: CO + OH => HCO + O
  prodf(26)=comp_ch_tabT(26,kf_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(1)**1.0)
  prodb(26)=comp_ch_tabT(26,kb_tab,Tint,Tdiff)*(coi(7)**1.0)*(coi(5)**1.0)
  ! reac n. 27: C2H4 (+M) <=> C2H2 + H2 (+M)
  M = sum(coi(1:nsc))
  prodf(27)=comp_ch_tabT(27,kf_tab,Tint,Tdiff)*(coi(12)**1.0)*M
  prodb(27)=comp_ch_tabT(27,kb_tab,Tint,Tdiff)*(coi(14)**1.0)*(coi(6)**1.0)*M
  ! reac n. 28: HO2 + O => O2 + OH
  prodf(28)=comp_ch_tabT(28,kf_tab,Tint,Tdiff)*(coi(11)**1.0)*(coi(5)**1.0)
  prodb(28)=comp_ch_tabT(28,kb_tab,Tint,Tdiff)*(coi(4)**1.0)*(coi(1)**1.0)
  ! reac n. 29: O2 + OH => HO2 + O
  prodf(29)=comp_ch_tabT(29,kf_tab,Tint,Tdiff)*(coi(4)**1.0)*(coi(1)**1.0)
  prodb(29)=comp_ch_tabT(29,kb_tab,Tint,Tdiff)*(coi(11)**1.0)*(coi(5)**1.0)
  ! reac n. 30: HCO + O2 => CO + HO2
  prodf(30)=comp_ch_tabT(30,kf_tab,Tint,Tdiff)*(coi(7)**1.0)*(coi(4)**1.0)
  prodb(30)=comp_ch_tabT(30,kb_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(11)**1.0)
  ! reac n. 31: CO + HO2 => HCO + O2
  prodf(31)=comp_ch_tabT(31,kf_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(11)**1.0)
  prodb(31)=comp_ch_tabT(31,kb_tab,Tint,Tdiff)*(coi(7)**1.0)*(coi(4)**1.0)
  ! reac n. 32: H + HO2 => 2 OH
  prodf(32)=comp_ch_tabT(32,kf_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(11)**1.0)
  prodb(32)=comp_ch_tabT(32,kb_tab,Tint,Tdiff)*(coi(1)**2.0)
  ! reac n. 33: 2 OH => H + HO2
  prodf(33)=comp_ch_tabT(33,kf_tab,Tint,Tdiff)*(coi(1)**2.0)
  prodb(33)=comp_ch_tabT(33,kb_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(11)**1.0)
  ! reac n. 34: H + HO2 => H2 + O2
  prodf(34)=comp_ch_tabT(34,kf_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(11)**1.0)
  prodb(34)=comp_ch_tabT(34,kb_tab,Tint,Tdiff)*(coi(6)**1.0)*(coi(4)**1.0)
  ! reac n. 35: H2 + O2 => H + HO2
  prodf(35)=comp_ch_tabT(35,kf_tab,Tint,Tdiff)*(coi(6)**1.0)*(coi(4)**1.0)
  prodb(35)=comp_ch_tabT(35,kb_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(11)**1.0)
  ! reac n. 36: HO2 + OH => H2O + O2
  prodf(36)=comp_ch_tabT(36,kf_tab,Tint,Tdiff)*(coi(11)**1.0)*(coi(1)**1.0)
  prodb(36)=comp_ch_tabT(36,kb_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(4)**1.0)
  ! reac n. 37: H2O + O2 => HO2 + OH
  prodf(37)=comp_ch_tabT(37,kf_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(4)**1.0)
  prodb(37)=comp_ch_tabT(37,kb_tab,Tint,Tdiff)*(coi(11)**1.0)*(coi(1)**1.0)
  ! reac n. 38: OH + M => H + O + M
  M=0.d0+coi(1)*1.0+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*2.5+coi(7)*1.0+coi(8)*1.9+coi(9)*12.0+coi(10)*3.8+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.0+coi(16)*1.0+coi(17)*1.0+coi(18)*1.0+coi(19)*1.0+coi(20)*1.0
  prodf(38)=comp_ch_tabT(38,kf_tab,Tint,Tdiff)*(coi(1)**1.0)*M
  prodb(38)=comp_ch_tabT(38,kb_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(5)**1.0)*M
  ! reac n. 39: H + O + M => OH + M
  M=0.d0+coi(1)*1.0+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*2.5+coi(7)*1.0+coi(8)*1.9+coi(9)*12.0+coi(10)*3.8+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.0+coi(16)*1.0+coi(17)*1.0+coi(18)*1.0+coi(19)*1.0+coi(20)*1.0
  prodf(39)=comp_ch_tabT(39,kf_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(5)**1.0)*M
  prodb(39)=comp_ch_tabT(39,kb_tab,Tint,Tdiff)*(coi(1)**1.0)*M
  ! reac n. 40: O2 + M => 2 O + M
  M=0.d0+coi(1)*1.0+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*2.5+coi(7)*1.0+coi(8)*1.9+coi(9)*12.0+coi(10)*3.8+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.0+coi(16)*1.0+coi(17)*1.0+coi(18)*1.0+coi(19)*1.0+coi(20)*1.0
  prodf(40)=comp_ch_tabT(40,kf_tab,Tint,Tdiff)*(coi(4)**1.0)*M
  prodb(40)=comp_ch_tabT(40,kb_tab,Tint,Tdiff)*(coi(5)**2.0)*M

  ! reac n. 41: 2 O + M => O2 + M
  M=0.d0+coi(1)*1.0+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*2.5+coi(7)*1.0+coi(8)*1.9+coi(9)*12.0+coi(10)*3.8+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.0+coi(16)*1.0+coi(17)*1.0+coi(18)*1.0+coi(19)*1.0+coi(20)*1.0
  prodf(41)=comp_ch_tabT(41,kf_tab,Tint,Tdiff)*(coi(5)**2.0)*M
  prodb(41)=comp_ch_tabT(41,kb_tab,Tint,Tdiff)*(coi(4)**1.0)*M
  ! reac n. 42: H2 + M => 2 H + M
  M=0.d0+coi(1)*1.0+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*2.5+coi(7)*1.0+coi(8)*1.9+coi(9)*12.0+coi(10)*3.8+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.0+coi(16)*1.0+coi(17)*1.0+coi(18)*1.0+coi(19)*1.0+coi(20)*1.0
  prodf(42)=comp_ch_tabT(42,kf_tab,Tint,Tdiff)*(coi(6)**1.0)*M
  prodb(42)=comp_ch_tabT(42,kb_tab,Tint,Tdiff)*(coi(3)**2.0)*M
  ! reac n. 43: 2 H + M => H2 + M
  M=0.d0+coi(1)*1.0+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*2.5+coi(7)*1.0+coi(8)*1.9+coi(9)*12.0+coi(10)*3.8+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.0+coi(16)*1.0+coi(17)*1.0+coi(18)*1.0+coi(19)*1.0+coi(20)*1.0
  prodf(43)=comp_ch_tabT(43,kf_tab,Tint,Tdiff)*(coi(3)**2.0)*M
  prodb(43)=comp_ch_tabT(43,kb_tab,Tint,Tdiff)*(coi(6)**1.0)*M
  ! reac n. 44: C2H3 + H (+M) <=> C2H4 (+M)
  M=sum(coi(1:20))
  prodf(44)=comp_ch_tabT(44,kf_tab,Tint,Tdiff)*(coi(15)**1.0)*(coi(3)**1.0)*M
  prodb(44)=comp_ch_tabT(44,kb_tab,Tint,Tdiff)*(coi(12)**1.0)*M
  ! reac n. 45: C2H2 + H (+M) <=> C2H3 (+M)
  M=0.d0+coi(1)*1.0+coi(2)*1.0+coi(3)*1.0+coi(4)*1.0+coi(5)*1.0+coi(6)*2.0+coi(7)*1.0+coi(8)*2.0+coi(9)*5.0+coi(10)*3.0+coi(11)*1.0+coi(12)*1.0+coi(13)*1.0+coi(14)*1.0+coi(15)*1.0+coi(16)*1.0+coi(17)*1.0+coi(18)*1.0+coi(19)*1.0+coi(20)*1.0
  prodf(45)=comp_ch_tabT(45,kf_tab,Tint,Tdiff)*(coi(14)**1.0)*(coi(3)**1.0)*M
  prodb(45)=comp_ch_tabT(45,kb_tab,Tint,Tdiff)*(coi(15)**1.0)*M
  ! reac n. 46: C2H4 + H => C2H3 + H2
  prodf(46)=comp_ch_tabT(46,kf_tab,Tint,Tdiff)*(coi(12)**1.0)*(coi(3)**1.0)
  prodb(46)=comp_ch_tabT(46,kb_tab,Tint,Tdiff)*(coi(15)**1.0)*(coi(6)**1.0)
  ! reac n. 47: C2H3 + H2 => C2H4 + H
  prodf(47)=comp_ch_tabT(47,kf_tab,Tint,Tdiff)*(coi(15)**1.0)*(coi(6)**1.0)
  prodb(47)=comp_ch_tabT(47,kb_tab,Tint,Tdiff)*(coi(12)**1.0)*(coi(3)**1.0)
  ! reac n. 48: C2H4 + OH => C2H3 + H2O
  prodf(48)=comp_ch_tabT(48,kf_tab,Tint,Tdiff)*(coi(12)**1.0)*(coi(1)**1.0)
  prodb(48)=comp_ch_tabT(48,kb_tab,Tint,Tdiff)*(coi(15)**1.0)*(coi(9)**1.0)
  ! reac n. 49: C2H3 + H2O => C2H4 + OH
  prodf(49)=comp_ch_tabT(49,kf_tab,Tint,Tdiff)*(coi(15)**1.0)*(coi(9)**1.0)
  prodb(49)=comp_ch_tabT(49,kb_tab,Tint,Tdiff)*(coi(12)**1.0)*(coi(1)**1.0)
  ! reac n. 50: C2H2 + O2 => HCCO + OH
  prodf(50)=comp_ch_tabT(50,kf_tab,Tint,Tdiff)*(coi(14)**1.0)*(coi(4)**1.0)
  prodb(50)=comp_ch_tabT(50,kb_tab,Tint,Tdiff)*(coi(16)**1.0)*(coi(1)**1.0)
  ! reac n. 51: HCCO + OH => C2H2 + O2
  prodf(51)=comp_ch_tabT(51,kf_tab,Tint,Tdiff)*(coi(16)**1.0)*(coi(1)**1.0)
  prodb(51)=comp_ch_tabT(51,kb_tab,Tint,Tdiff)*(coi(14)**1.0)*(coi(4)**1.0)
  ! reac n. 52: CH2 + O2 => CO + H2O
  prodf(52)=comp_ch_tabT(52,kf_tab,Tint,Tdiff)*(coi(13)**1.0)*(coi(4)**1.0)
  prodb(52)=comp_ch_tabT(52,kb_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(9)**1.0)
  ! reac n. 53: CO + H2O => CH2 + O2
  prodf(53)=comp_ch_tabT(53,kf_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(9)**1.0)
  prodb(53)=comp_ch_tabT(53,kb_tab,Tint,Tdiff)*(coi(13)**1.0)*(coi(4)**1.0)
  ! reac n. 54: C2H2 + O => CH2 + CO
  prodf(54)=comp_ch_tabT(54,kf_tab,Tint,Tdiff)*(coi(14)**1.0)*(coi(5)**1.0)
  prodb(54)=comp_ch_tabT(54,kb_tab,Tint,Tdiff)*(coi(13)**1.0)*(coi(8)**1.0)
  ! reac n. 55: CH2 + CO => C2H2 + O
  prodf(55)=comp_ch_tabT(55,kf_tab,Tint,Tdiff)*(coi(13)**1.0)*(coi(8)**1.0)
  prodb(55)=comp_ch_tabT(55,kb_tab,Tint,Tdiff)*(coi(14)**1.0)*(coi(5)**1.0)
  ! reac n. 56: CH2 + O2 => HCO + OH
  prodf(56)=comp_ch_tabT(56,kf_tab,Tint,Tdiff)*(coi(13)**1.0)*(coi(4)**1.0)
  prodb(56)=comp_ch_tabT(56,kb_tab,Tint,Tdiff)*(coi(7)**1.0)*(coi(1)**1.0)
  ! reac n. 57: HCO + OH => CH2 + O2
  prodf(57)=comp_ch_tabT(57,kf_tab,Tint,Tdiff)*(coi(7)**1.0)*(coi(1)**1.0)
  prodb(57)=comp_ch_tabT(57,kb_tab,Tint,Tdiff)*(coi(13)**1.0)*(coi(4)**1.0)
  ! reac n. 58: CH2 + O => CO + 2 H
  prodf(58)=comp_ch_tabT(58,kf_tab,Tint,Tdiff)*(coi(13)**1.0)*(coi(5)**1.0)
  prodb(58)=comp_ch_tabT(58,kb_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(3)**2.0)
  ! reac n. 59: CH2 + O2 => CO2 + 2 H
  prodf(59)=comp_ch_tabT(59,kf_tab,Tint,Tdiff)*(coi(13)**1.0)*(coi(4)**1.0)
  prodb(59)=comp_ch_tabT(59,kb_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(3)**2.0)
  ! reac n. 60: C2H3 + O2 => C2H2 + HO2
  prodf(60)=comp_ch_tabT(60,kf_tab,Tint,Tdiff)*(coi(15)**1.0)*(coi(4)**1.0)
  prodb(60)=comp_ch_tabT(60,kb_tab,Tint,Tdiff)*(coi(14)**1.0)*(coi(11)**1.0)
  ! reac n. 61: C2H2 + HO2 => C2H3 + O2
  prodf(61)=comp_ch_tabT(61,kf_tab,Tint,Tdiff)*(coi(14)**1.0)*(coi(11)**1.0)
  prodb(61)=comp_ch_tabT(61,kb_tab,Tint,Tdiff)*(coi(15)**1.0)*(coi(4)**1.0)
  ! reac n. 62: C2H2 + O => H + HCCO
  prodf(62)=comp_ch_tabT(62,kf_tab,Tint,Tdiff)*(coi(14)**1.0)*(coi(5)**1.0)
  prodb(62)=comp_ch_tabT(62,kb_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(16)**1.0)
  ! reac n. 63: H + HCCO => C2H2 + O
  prodf(63)=comp_ch_tabT(63,kf_tab,Tint,Tdiff)*(coi(3)**1.0)*(coi(16)**1.0)
  prodb(63)=comp_ch_tabT(63,kb_tab,Tint,Tdiff)*(coi(14)**1.0)*(coi(5)**1.0)
  ! reac n. 64: C2H2 + OH => CH2CO + H
  prodf(64)=comp_ch_tabT(64,kf_tab,Tint,Tdiff)*(coi(14)**1.0)*(coi(1)**1.0)
  prodb(64)=comp_ch_tabT(64,kb_tab,Tint,Tdiff)*(coi(17)**1.0)*(coi(3)**1.0)
  ! reac n. 65: CH2CO + H => C2H2 + OH
  prodf(65)=comp_ch_tabT(65,kf_tab,Tint,Tdiff)*(coi(17)**1.0)*(coi(3)**1.0)
  prodb(65)=comp_ch_tabT(65,kb_tab,Tint,Tdiff)*(coi(14)**1.0)*(coi(1)**1.0)
  ! reac n. 66: CH2CO + O => CH2 + CO2
  prodf(66)=comp_ch_tabT(66,kf_tab,Tint,Tdiff)*(coi(17)**1.0)*(coi(5)**1.0)
  prodb(66)=comp_ch_tabT(66,kb_tab,Tint,Tdiff)*(coi(13)**1.0)*(coi(10)**1.0)
  ! reac n. 67: CH2 + CO2 => CH2CO + O
  prodf(67)=comp_ch_tabT(67,kf_tab,Tint,Tdiff)*(coi(13)**1.0)*(coi(10)**1.0)
  prodb(67)=comp_ch_tabT(67,kb_tab,Tint,Tdiff)*(coi(17)**1.0)*(coi(5)**1.0)
  ! reac n. 68: CH2CO (+M) <=> CH2 + CO (+M)
  M = sum(coi(1:nsc))
  prodf(68)=comp_ch_tabT(68,kf_tab,Tint,Tdiff)*(coi(17)**1.0)*M
  prodb(68)=comp_ch_tabT(68,kb_tab,Tint,Tdiff)*(coi(13)**1.0)*(coi(8)**1.0)*M
  ! reac n. 69: CH2CO + O => HCCO + OH
  prodf(69)=comp_ch_tabT(69,kf_tab,Tint,Tdiff)*(coi(17)**1.0)*(coi(5)**1.0)
  prodb(69)=comp_ch_tabT(69,kb_tab,Tint,Tdiff)*(coi(16)**1.0)*(coi(1)**1.0)
  ! reac n. 70: HCCO + OH => CH2CO + O
  prodf(70)=comp_ch_tabT(70,kf_tab,Tint,Tdiff)*(coi(16)**1.0)*(coi(1)**1.0)
  prodb(70)=comp_ch_tabT(70,kb_tab,Tint,Tdiff)*(coi(17)**1.0)*(coi(5)**1.0)
  ! reac n. 71: CH2CO + OH => H2O + HCCO
  prodf(71)=comp_ch_tabT(71,kf_tab,Tint,Tdiff)*(coi(17)**1.0)*(coi(1)**1.0)
  prodb(71)=comp_ch_tabT(71,kb_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(16)**1.0)
  ! reac n. 72: H2O + HCCO => CH2CO + OH
  prodf(72)=comp_ch_tabT(72,kf_tab,Tint,Tdiff)*(coi(9)**1.0)*(coi(16)**1.0)
  prodb(72)=comp_ch_tabT(72,kb_tab,Tint,Tdiff)*(coi(17)**1.0)*(coi(1)**1.0)
  ! reac n. 73: CH2CO + H => H2 + HCCO
  prodf(73)=comp_ch_tabT(73,kf_tab,Tint,Tdiff)*(coi(17)**1.0)*(coi(3)**1.0)
  prodb(73)=comp_ch_tabT(73,kb_tab,Tint,Tdiff)*(coi(6)**1.0)*(coi(16)**1.0)
  ! reac n. 74: H2 + HCCO => CH2CO + H
  prodf(74)=comp_ch_tabT(74,kf_tab,Tint,Tdiff)*(coi(6)**1.0)*(coi(16)**1.0)
  prodb(74)=comp_ch_tabT(74,kb_tab,Tint,Tdiff)*(coi(17)**1.0)*(coi(3)**1.0)
  ! reac n. 75: HCCO + OH => 2 HCO
  prodf(75)=comp_ch_tabT(75,kf_tab,Tint,Tdiff)*(coi(16)**1.0)*(coi(1)**1.0)
  prodb(75)=comp_ch_tabT(75,kb_tab,Tint,Tdiff)*(coi(7)**2.0)
  ! reac n. 76: 2 HCO => HCCO + OH
  prodf(76)=comp_ch_tabT(76,kf_tab,Tint,Tdiff)*(coi(7)**2.0)
  prodb(76)=comp_ch_tabT(76,kb_tab,Tint,Tdiff)*(coi(16)**1.0)*(coi(1)**1.0)
  ! reac n. 77: HCCO + O => 2 CO + H
  prodf(77)=comp_ch_tabT(77,kf_tab,Tint,Tdiff)*(coi(16)**1.0)*(coi(5)**1.0)
  prodb(77)=comp_ch_tabT(77,kb_tab,Tint,Tdiff)*(coi(8)**2.0)*(coi(3)**1.0)
  ! reac n. 78: CH2 + O2 => CO2 + H2
  prodf(78)=comp_ch_tabT(78,kf_tab,Tint,Tdiff)*(coi(13)**1.0)*(coi(4)**1.0)
  prodb(78)=comp_ch_tabT(78,kb_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(6)**1.0)
  ! reac n. 79: CO2 + H2 => CH2 + O2
  prodf(79)=comp_ch_tabT(79,kf_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(6)**1.0)
  prodb(79)=comp_ch_tabT(79,kb_tab,Tint,Tdiff)*(coi(13)**1.0)*(coi(4)**1.0)
  ! reac n. 80: C2H3 + H => C2H2 + H2
  prodf(80)=comp_ch_tabT(80,kf_tab,Tint,Tdiff)*(coi(15)**1.0)*(coi(3)**1.0)
  prodb(80)=comp_ch_tabT(80,kb_tab,Tint,Tdiff)*(coi(14)**1.0)*(coi(6)**1.0)
  ! reac n. 81: C2H2 + H2 => C2H3 + H
  prodf(81)=comp_ch_tabT(81,kf_tab,Tint,Tdiff)*(coi(14)**1.0)*(coi(6)**1.0)
  prodb(81)=comp_ch_tabT(81,kb_tab,Tint,Tdiff)*(coi(15)**1.0)*(coi(3)**1.0)
  ! reac n. 82: C4H6 => 2 C2H3
  prodf(82)=comp_ch_tabT(82,kf_tab,Tint,Tdiff)*(coi(2)**1.0)
  prodb(82)=comp_ch_tabT(82,kb_tab,Tint,Tdiff)*(coi(15)**2.0)
  ! reac n. 83: 2 C2H3 => C4H6
  prodf(83)=comp_ch_tabT(83,kf_tab,Tint,Tdiff)*(coi(15)**2.0)
  prodb(83)=comp_ch_tabT(83,kb_tab,Tint,Tdiff)*(coi(2)**1.0)
  ! reac n. 84: C4H6 + O => C2H4 + CH2CO
  prodf(84)=comp_ch_tabT(84,kf_tab,Tint,Tdiff)*(coi(2)**1.0)*(coi(5)**1.0)
  prodb(84)=comp_ch_tabT(84,kb_tab,Tint,Tdiff)*(coi(12)**1.0)*(coi(17)**1.0)
  ! reac n. 85: C2H4 + CH2CO => C4H6 + O
  prodf(85)=comp_ch_tabT(85,kf_tab,Tint,Tdiff)*(coi(12)**1.0)*(coi(17)**1.0)
  prodb(85)=comp_ch_tabT(85,kb_tab,Tint,Tdiff)*(coi(2)**1.0)*(coi(5)**1.0)
  ! reac n. 86: C2H4 + O2 => C2H3 + HO2
  prodf(86)=comp_ch_tabT(86,kf_tab,Tint,Tdiff)*(coi(12)**1.0)*(coi(4)**1.0)
  prodb(86)=comp_ch_tabT(86,kb_tab,Tint,Tdiff)*(coi(15)**1.0)*(coi(11)**1.0)
  ! reac n. 87: C2H3 + HO2 => C2H4 + O2
  prodf(87)=comp_ch_tabT(87,kf_tab,Tint,Tdiff)*(coi(15)**1.0)*(coi(11)**1.0)
  prodb(87)=comp_ch_tabT(87,kb_tab,Tint,Tdiff)*(coi(12)**1.0)*(coi(4)**1.0)
  ! reac n. 88: HCO + O => CO2 + H
  prodf(88)=comp_ch_tabT(88,kf_tab,Tint,Tdiff)*(coi(7)**1.0)*(coi(5)**1.0)
  prodb(88)=comp_ch_tabT(88,kb_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(3)**1.0)
  ! reac n. 89: CO2 + H => HCO + O
  prodf(89)=comp_ch_tabT(89,kf_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(3)**1.0)
  prodb(89)=comp_ch_tabT(89,kb_tab,Tint,Tdiff)*(coi(7)**1.0)*(coi(5)**1.0)
  ! reac n. 90: C2H3 + C2H4 => C4H6 + H
  prodf(90)=comp_ch_tabT(90,kf_tab,Tint,Tdiff)*(coi(15)**1.0)*(coi(12)**1.0)
  prodb(90)=comp_ch_tabT(90,kb_tab,Tint,Tdiff)*(coi(2)**1.0)*(coi(3)**1.0)
  ! reac n. 91: C4H6 + H => C2H3 + C2H4
  prodf(91)=comp_ch_tabT(91,kf_tab,Tint,Tdiff)*(coi(2)**1.0)*(coi(3)**1.0)
  prodb(91)=comp_ch_tabT(91,kb_tab,Tint,Tdiff)*(coi(15)**1.0)*(coi(12)**1.0)
  ! reac n. 92: C3H2 + O2 => CO + H + HCCO
  prodf(92)=comp_ch_tabT(92,kf_tab,Tint,Tdiff)*(coi(19)**1.0)*(coi(4)**1.0)
  prodb(92)=comp_ch_tabT(92,kb_tab,Tint,Tdiff)*(coi(8)**1.0)*(coi(3)**1.0)*(coi(16)**1.0)
  ! reac n. 93: C2H3CO => C2H3 + CO
  prodf(93)=comp_ch_tabT(93,kf_tab,Tint,Tdiff)*(coi(18)**1.0)
  prodb(93)=comp_ch_tabT(93,kb_tab,Tint,Tdiff)*(coi(15)**1.0)*(coi(8)**1.0)
  ! reac n. 94: C2H3 + CO => C2H3CO
  prodf(94)=comp_ch_tabT(94,kf_tab,Tint,Tdiff)*(coi(15)**1.0)*(coi(8)**1.0)
  prodb(94)=comp_ch_tabT(94,kb_tab,Tint,Tdiff)*(coi(18)**1.0)
  ! reac n. 95: C3H3 + H => C3H2 + H2
  prodf(95)=comp_ch_tabT(95,kf_tab,Tint,Tdiff)*(coi(20)**1.0)*(coi(3)**1.0)
  prodb(95)=comp_ch_tabT(95,kb_tab,Tint,Tdiff)*(coi(19)**1.0)*(coi(6)**1.0)
  ! reac n. 96: C3H2 + H2 => C3H3 + H
  prodf(96)=comp_ch_tabT(96,kf_tab,Tint,Tdiff)*(coi(19)**1.0)*(coi(6)**1.0)
  prodb(96)=comp_ch_tabT(96,kb_tab,Tint,Tdiff)*(coi(20)**1.0)*(coi(3)**1.0)
  ! reac n. 97: C3H2 + OH => C2H2 + HCO
  prodf(97)=comp_ch_tabT(97,kf_tab,Tint,Tdiff)*(coi(19)**1.0)*(coi(1)**1.0)
  prodb(97)=comp_ch_tabT(97,kb_tab,Tint,Tdiff)*(coi(14)**1.0)*(coi(7)**1.0)
  ! reac n. 98: C2H2 + HCO => C3H2 + OH
  prodf(98)=comp_ch_tabT(98,kf_tab,Tint,Tdiff)*(coi(14)**1.0)*(coi(7)**1.0)
  prodb(98)=comp_ch_tabT(98,kb_tab,Tint,Tdiff)*(coi(19)**1.0)*(coi(1)**1.0)
  ! reac n. 99: C3H3 + OH => C3H2 + H2O
  prodf(99)=comp_ch_tabT(99,kf_tab,Tint,Tdiff)*(coi(20)**1.0)*(coi(1)**1.0)
  prodb(99)=comp_ch_tabT(99,kb_tab,Tint,Tdiff)*(coi(19)**1.0)*(coi(9)**1.0)
  ! reac n. 100: C3H2 + H2O => C3H3 + OH
  prodf(100)=comp_ch_tabT(100,kf_tab,Tint,Tdiff)*(coi(19)**1.0)*(coi(9)**1.0)
  prodb(100)=comp_ch_tabT(100,kb_tab,Tint,Tdiff)*(coi(20)**1.0)*(coi(1)**1.0)
  ! reac n. 101: C3H3 + O2 => CH2CO + HCO
  prodf(101)=comp_ch_tabT(101,kf_tab,Tint,Tdiff)*(coi(20)**1.0)*(coi(4)**1.0)
  prodb(101)=comp_ch_tabT(101,kb_tab,Tint,Tdiff)*(coi(17)**1.0)*(coi(7)**1.0)
  ! reac n. 102: CH2CO + HCO => C3H3 + O2
  prodf(102)=comp_ch_tabT(102,kf_tab,Tint,Tdiff)*(coi(17)**1.0)*(coi(7)**1.0)
  prodb(102)=comp_ch_tabT(102,kb_tab,Tint,Tdiff)*(coi(20)**1.0)*(coi(4)**1.0)
  ! reac n. 103: HCCO + O2 => CO2 + HCO
  prodf(103)=comp_ch_tabT(103,kf_tab,Tint,Tdiff)*(coi(16)**1.0)*(coi(4)**1.0)
  prodb(103)=comp_ch_tabT(103,kb_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(7)**1.0)
  ! reac n. 104: CO2 + HCO => HCCO + O2
  prodf(104)=comp_ch_tabT(104,kf_tab,Tint,Tdiff)*(coi(10)**1.0)*(coi(7)**1.0)
  prodb(104)=comp_ch_tabT(104,kb_tab,Tint,Tdiff)*(coi(16)**1.0)*(coi(4)**1.0)
  ! species source terms
  omegadot(1)=Wm_tab(1)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(7)-prodb(7))+(0.0-1.0)*(prodf(8)-prodb(8))+(2.0-0.0)*(prodf(9)-prodb(9))+(0.0-2.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(12)-prodb(12))+(1.0-0.0)*(prodf(15)-prodb(15))+(0.0-1.0)*(prodf(16)-prodb(16))+(1.0-0.0)*(prodf(17)-prodb(17))+(0.0-1.0)*(prodf(18)-prodb(18))+(1.0-0.0)*(prodf(25)-prodb(25))+(0.0-1.0)*(prodf(26)-prodb(26))+(1.0-0.0)*(prodf(28)-prodb(28))+(0.0-1.0)*(prodf(29)-prodb(29))+(2.0-0.0)*(prodf(32)-prodb(32))+(0.0-2.0)*(prodf(33)-prodb(33))+(0.0-1.0)*(prodf(36)-prodb(36))+(1.0-0.0)*(prodf(37)-prodb(37))+(0.0-1.0)*(prodf(38)-prodb(38))+(1.0-0.0)*(prodf(39)-prodb(39))+(0.0-1.0)*(prodf(48)-prodb(48))+(1.0-0.0)*(prodf(49)-prodb(49))+(1.0-0.0)*(prodf(50)-prodb(50))+(0.0-1.0)*(prodf(51)-prodb(51))+(1.0-0.0)*(prodf(56)-prodb(56))+(0.0-1.0)*(prodf(57)-prodb(57))+(0.0-1.0)*(prodf(64)-prodb(64))+(1.0-0.0)*(prodf(65)-prodb(65))+(1.0-0.0)*(prodf(69)-prodb(69))+(0.0-1.0)*(prodf(70)-prodb(70))+(0.0-1.0)*(prodf(71)-prodb(71))+(1.0-0.0)*(prodf(72)-prodb(72))+(0.0-1.0)*(prodf(75)-prodb(75))+(1.0-0.0)*(prodf(76)-prodb(76))+(0.0-1.0)*(prodf(97)-prodb(97))+(1.0-0.0)*(prodf(98)-prodb(98))+(0.0-1.0)*(prodf(99)-prodb(99))+(1.0-0.0)*(prodf(100)-prodb(100)))
  omegadot(2)=Wm_tab(2)*(+(0.0-1.0)*(prodf(82)-prodb(82))+(1.0-0.0)*(prodf(83)-prodb(83))+(0.0-1.0)*(prodf(84)-prodb(84))+(1.0-0.0)*(prodf(85)-prodb(85))+(1.0-0.0)*(prodf(90)-prodb(90))+(0.0-1.0)*(prodf(91)-prodb(91)))
  omegadot(3)=Wm_tab(3)*(+(1.0-0.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(4)-prodb(4))+(0.0-1.0)*(prodf(5)-prodb(5))+(1.0-0.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(7)-prodb(7))+(0.0-1.0)*(prodf(8)-prodb(8))+(1.0-0.0)*(prodf(11)-prodb(11))+(0.0-1.0)*(prodf(12)-prodb(12))+(1.0-0.0)*(prodf(13)-prodb(13))+(0.0-1.0)*(prodf(14)-prodb(14))+(1.0-0.0)*(prodf(17)-prodb(17))+(0.0-1.0)*(prodf(18)-prodb(18))+(0.0-1.0)*(prodf(19)-prodb(19))+(0.0-1.0)*(prodf(23)-prodb(23))+(1.0-0.0)*(prodf(24)-prodb(24))+(0.0-1.0)*(prodf(32)-prodb(32))+(1.0-0.0)*(prodf(33)-prodb(33))+(0.0-1.0)*(prodf(34)-prodb(34))+(1.0-0.0)*(prodf(35)-prodb(35))+(1.0-0.0)*(prodf(38)-prodb(38))+(0.0-1.0)*(prodf(39)-prodb(39))+(2.0-0.0)*(prodf(42)-prodb(42))+(0.0-2.0)*(prodf(43)-prodb(43))+(0.0-1.0)*(prodf(44)-prodb(44))+(0.0-1.0)*(prodf(45)-prodb(45))+(0.0-1.0)*(prodf(46)-prodb(46))+(1.0-0.0)*(prodf(47)-prodb(47))+(2.0-0.0)*(prodf(58)-prodb(58))+(2.0-0.0)*(prodf(59)-prodb(59))+(1.0-0.0)*(prodf(62)-prodb(62))+(0.0-1.0)*(prodf(63)-prodb(63))+(1.0-0.0)*(prodf(64)-prodb(64))+(0.0-1.0)*(prodf(65)-prodb(65))+(0.0-1.0)*(prodf(73)-prodb(73))+(1.0-0.0)*(prodf(74)-prodb(74))+(1.0-0.0)*(prodf(77)-prodb(77))+(0.0-1.0)*(prodf(80)-prodb(80))+(1.0-0.0)*(prodf(81)-prodb(81))+(1.0-0.0)*(prodf(88)-prodb(88))+(0.0-1.0)*(prodf(89)-prodb(89))+(1.0-0.0)*(prodf(90)-prodb(90))+(0.0-1.0)*(prodf(91)-prodb(91))+(1.0-0.0)*(prodf(92)-prodb(92))+(0.0-1.0)*(prodf(95)-prodb(95))+(1.0-0.0)*(prodf(96)-prodb(96)))
  omegadot(4)=Wm_tab(4)*(+(0.0-1.0)*(prodf(5)-prodb(5))+(1.0-0.0)*(prodf(6)-prodb(6))+(0.0-1.0)*(prodf(19)-prodb(19))+(0.0-1.0)*(prodf(21)-prodb(21))+(1.0-0.0)*(prodf(22)-prodb(22))+(1.0-0.0)*(prodf(28)-prodb(28))+(0.0-1.0)*(prodf(29)-prodb(29))+(0.0-1.0)*(prodf(30)-prodb(30))+(1.0-0.0)*(prodf(31)-prodb(31))+(1.0-0.0)*(prodf(34)-prodb(34))+(0.0-1.0)*(prodf(35)-prodb(35))+(1.0-0.0)*(prodf(36)-prodb(36))+(0.0-1.0)*(prodf(37)-prodb(37))+(0.0-1.0)*(prodf(40)-prodb(40))+(1.0-0.0)*(prodf(41)-prodb(41))+(0.0-1.0)*(prodf(50)-prodb(50))+(1.0-0.0)*(prodf(51)-prodb(51))+(0.0-1.0)*(prodf(52)-prodb(52))+(1.0-0.0)*(prodf(53)-prodb(53))+(0.0-1.0)*(prodf(56)-prodb(56))+(1.0-0.0)*(prodf(57)-prodb(57))+(0.0-1.0)*(prodf(59)-prodb(59))+(0.0-1.0)*(prodf(60)-prodb(60))+(1.0-0.0)*(prodf(61)-prodb(61))+(0.0-1.0)*(prodf(78)-prodb(78))+(1.0-0.0)*(prodf(79)-prodb(79))+(0.0-1.0)*(prodf(86)-prodb(86))+(1.0-0.0)*(prodf(87)-prodb(87))+(0.0-1.0)*(prodf(92)-prodb(92))+(0.0-1.0)*(prodf(101)-prodb(101))+(1.0-0.0)*(prodf(102)-prodb(102))+(0.0-1.0)*(prodf(103)-prodb(103))+(1.0-0.0)*(prodf(104)-prodb(104)))
  omegadot(5)=Wm_tab(5)*(+(1.0-0.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(6)-prodb(6))+(0.0-1.0)*(prodf(7)-prodb(7))+(1.0-0.0)*(prodf(8)-prodb(8))+(0.0-1.0)*(prodf(9)-prodb(9))+(1.0-0.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(20)-prodb(20))+(1.0-0.0)*(prodf(21)-prodb(21))+(0.0-1.0)*(prodf(22)-prodb(22))+(0.0-1.0)*(prodf(25)-prodb(25))+(1.0-0.0)*(prodf(26)-prodb(26))+(0.0-1.0)*(prodf(28)-prodb(28))+(1.0-0.0)*(prodf(29)-prodb(29))+(1.0-0.0)*(prodf(38)-prodb(38))+(0.0-1.0)*(prodf(39)-prodb(39))+(2.0-0.0)*(prodf(40)-prodb(40))+(0.0-2.0)*(prodf(41)-prodb(41))+(0.0-1.0)*(prodf(54)-prodb(54))+(1.0-0.0)*(prodf(55)-prodb(55))+(0.0-1.0)*(prodf(58)-prodb(58))+(0.0-1.0)*(prodf(62)-prodb(62))+(1.0-0.0)*(prodf(63)-prodb(63))+(0.0-1.0)*(prodf(66)-prodb(66))+(1.0-0.0)*(prodf(67)-prodb(67))+(0.0-1.0)*(prodf(69)-prodb(69))+(1.0-0.0)*(prodf(70)-prodb(70))+(0.0-1.0)*(prodf(77)-prodb(77))+(0.0-1.0)*(prodf(84)-prodb(84))+(1.0-0.0)*(prodf(85)-prodb(85))+(0.0-1.0)*(prodf(88)-prodb(88))+(1.0-0.0)*(prodf(89)-prodb(89)))
  omegadot(6)=Wm_tab(6)*(+(0.0-1.0)*(prodf(7)-prodb(7))+(1.0-0.0)*(prodf(8)-prodb(8))+(0.0-1.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(12)-prodb(12))+(1.0-0.0)*(prodf(23)-prodb(23))+(0.0-1.0)*(prodf(24)-prodb(24))+(1.0-0.0)*(prodf(27)-prodb(27))+(1.0-0.0)*(prodf(34)-prodb(34))+(0.0-1.0)*(prodf(35)-prodb(35))+(0.0-1.0)*(prodf(42)-prodb(42))+(1.0-0.0)*(prodf(43)-prodb(43))+(1.0-0.0)*(prodf(46)-prodb(46))+(0.0-1.0)*(prodf(47)-prodb(47))+(1.0-0.0)*(prodf(73)-prodb(73))+(0.0-1.0)*(prodf(74)-prodb(74))+(1.0-0.0)*(prodf(78)-prodb(78))+(0.0-1.0)*(prodf(79)-prodb(79))+(1.0-0.0)*(prodf(80)-prodb(80))+(0.0-1.0)*(prodf(81)-prodb(81))+(1.0-0.0)*(prodf(95)-prodb(95))+(0.0-1.0)*(prodf(96)-prodb(96)))
  omegadot(7)=Wm_tab(7)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(13)-prodb(13))+(1.0-0.0)*(prodf(14)-prodb(14))+(0.0-1.0)*(prodf(23)-prodb(23))+(1.0-0.0)*(prodf(24)-prodb(24))+(0.0-1.0)*(prodf(25)-prodb(25))+(1.0-0.0)*(prodf(26)-prodb(26))+(0.0-1.0)*(prodf(30)-prodb(30))+(1.0-0.0)*(prodf(31)-prodb(31))+(1.0-0.0)*(prodf(56)-prodb(56))+(0.0-1.0)*(prodf(57)-prodb(57))+(2.0-0.0)*(prodf(75)-prodb(75))+(0.0-2.0)*(prodf(76)-prodb(76))+(0.0-1.0)*(prodf(88)-prodb(88))+(1.0-0.0)*(prodf(89)-prodb(89))+(1.0-0.0)*(prodf(97)-prodb(97))+(0.0-1.0)*(prodf(98)-prodb(98))+(1.0-0.0)*(prodf(101)-prodb(101))+(0.0-1.0)*(prodf(102)-prodb(102))+(1.0-0.0)*(prodf(103)-prodb(103))+(0.0-1.0)*(prodf(104)-prodb(104)))
  omegadot(8)=Wm_tab(8)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(13)-prodb(13))+(0.0-1.0)*(prodf(14)-prodb(14))+(0.0-1.0)*(prodf(15)-prodb(15))+(1.0-0.0)*(prodf(16)-prodb(16))+(0.0-1.0)*(prodf(20)-prodb(20))+(0.0-1.0)*(prodf(21)-prodb(21))+(1.0-0.0)*(prodf(22)-prodb(22))+(1.0-0.0)*(prodf(23)-prodb(23))+(0.0-1.0)*(prodf(24)-prodb(24))+(1.0-0.0)*(prodf(25)-prodb(25))+(0.0-1.0)*(prodf(26)-prodb(26))+(1.0-0.0)*(prodf(30)-prodb(30))+(0.0-1.0)*(prodf(31)-prodb(31))+(1.0-0.0)*(prodf(52)-prodb(52))+(0.0-1.0)*(prodf(53)-prodb(53))+(1.0-0.0)*(prodf(54)-prodb(54))+(0.0-1.0)*(prodf(55)-prodb(55))+(1.0-0.0)*(prodf(58)-prodb(58))+(1.0-0.0)*(prodf(68)-prodb(68))+(2.0-0.0)*(prodf(77)-prodb(77))+(1.0-0.0)*(prodf(92)-prodb(92))+(1.0-0.0)*(prodf(93)-prodb(93))+(0.0-1.0)*(prodf(94)-prodb(94)))
  omegadot(9)=Wm_tab(9)*(+(1.0-0.0)*(prodf(1)-prodb(1))+(0.0-1.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(9)-prodb(9))+(1.0-0.0)*(prodf(10)-prodb(10))+(1.0-0.0)*(prodf(11)-prodb(11))+(0.0-1.0)*(prodf(12)-prodb(12))+(0.0-1.0)*(prodf(17)-prodb(17))+(1.0-0.0)*(prodf(18)-prodb(18))+(1.0-0.0)*(prodf(36)-prodb(36))+(0.0-1.0)*(prodf(37)-prodb(37))+(1.0-0.0)*(prodf(48)-prodb(48))+(0.0-1.0)*(prodf(49)-prodb(49))+(1.0-0.0)*(prodf(52)-prodb(52))+(0.0-1.0)*(prodf(53)-prodb(53))+(1.0-0.0)*(prodf(71)-prodb(71))+(0.0-1.0)*(prodf(72)-prodb(72))+(1.0-0.0)*(prodf(99)-prodb(99))+(0.0-1.0)*(prodf(100)-prodb(100)))
  omegadot(10)=Wm_tab(10)*(+(1.0-0.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(15)-prodb(15))+(0.0-1.0)*(prodf(16)-prodb(16))+(1.0-0.0)*(prodf(20)-prodb(20))+(1.0-0.0)*(prodf(21)-prodb(21))+(0.0-1.0)*(prodf(22)-prodb(22))+(1.0-0.0)*(prodf(59)-prodb(59))+(1.0-0.0)*(prodf(66)-prodb(66))+(0.0-1.0)*(prodf(67)-prodb(67))+(1.0-0.0)*(prodf(78)-prodb(78))+(0.0-1.0)*(prodf(79)-prodb(79))+(1.0-0.0)*(prodf(88)-prodb(88))+(0.0-1.0)*(prodf(89)-prodb(89))+(1.0-0.0)*(prodf(103)-prodb(103))+(0.0-1.0)*(prodf(104)-prodb(104)))
  omegadot(11)=Wm_tab(11)*(+(0.0-1.0)*(prodf(15)-prodb(15))+(1.0-0.0)*(prodf(16)-prodb(16))+(1.0-0.0)*(prodf(19)-prodb(19))+(0.0-1.0)*(prodf(28)-prodb(28))+(1.0-0.0)*(prodf(29)-prodb(29))+(1.0-0.0)*(prodf(30)-prodb(30))+(0.0-1.0)*(prodf(31)-prodb(31))+(0.0-1.0)*(prodf(32)-prodb(32))+(1.0-0.0)*(prodf(33)-prodb(33))+(0.0-1.0)*(prodf(34)-prodb(34))+(1.0-0.0)*(prodf(35)-prodb(35))+(0.0-1.0)*(prodf(36)-prodb(36))+(1.0-0.0)*(prodf(37)-prodb(37))+(1.0-0.0)*(prodf(60)-prodb(60))+(0.0-1.0)*(prodf(61)-prodb(61))+(1.0-0.0)*(prodf(86)-prodb(86))+(0.0-1.0)*(prodf(87)-prodb(87)))
  omegadot(12)=Wm_tab(12)*(+(0.0-1.0)*(prodf(27)-prodb(27))+(1.0-0.0)*(prodf(44)-prodb(44))+(0.0-1.0)*(prodf(46)-prodb(46))+(1.0-0.0)*(prodf(47)-prodb(47))+(0.0-1.0)*(prodf(48)-prodb(48))+(1.0-0.0)*(prodf(49)-prodb(49))+(1.0-0.0)*(prodf(84)-prodb(84))+(0.0-1.0)*(prodf(85)-prodb(85))+(0.0-1.0)*(prodf(86)-prodb(86))+(1.0-0.0)*(prodf(87)-prodb(87))+(0.0-1.0)*(prodf(90)-prodb(90))+(1.0-0.0)*(prodf(91)-prodb(91)))
  omegadot(13)=Wm_tab(13)*(+(0.0-1.0)*(prodf(52)-prodb(52))+(1.0-0.0)*(prodf(53)-prodb(53))+(1.0-0.0)*(prodf(54)-prodb(54))+(0.0-1.0)*(prodf(55)-prodb(55))+(0.0-1.0)*(prodf(56)-prodb(56))+(1.0-0.0)*(prodf(57)-prodb(57))+(0.0-1.0)*(prodf(58)-prodb(58))+(0.0-1.0)*(prodf(59)-prodb(59))+(1.0-0.0)*(prodf(66)-prodb(66))+(0.0-1.0)*(prodf(67)-prodb(67))+(1.0-0.0)*(prodf(68)-prodb(68))+(0.0-1.0)*(prodf(78)-prodb(78))+(1.0-0.0)*(prodf(79)-prodb(79)))
  omegadot(14)=Wm_tab(14)*(+(1.0-0.0)*(prodf(27)-prodb(27))+(0.0-1.0)*(prodf(45)-prodb(45))+(0.0-1.0)*(prodf(50)-prodb(50))+(1.0-0.0)*(prodf(51)-prodb(51))+(0.0-1.0)*(prodf(54)-prodb(54))+(1.0-0.0)*(prodf(55)-prodb(55))+(1.0-0.0)*(prodf(60)-prodb(60))+(0.0-1.0)*(prodf(61)-prodb(61))+(0.0-1.0)*(prodf(62)-prodb(62))+(1.0-0.0)*(prodf(63)-prodb(63))+(0.0-1.0)*(prodf(64)-prodb(64))+(1.0-0.0)*(prodf(65)-prodb(65))+(1.0-0.0)*(prodf(80)-prodb(80))+(0.0-1.0)*(prodf(81)-prodb(81))+(1.0-0.0)*(prodf(97)-prodb(97))+(0.0-1.0)*(prodf(98)-prodb(98)))
  omegadot(15)=Wm_tab(15)*(+(0.0-1.0)*(prodf(44)-prodb(44))+(1.0-0.0)*(prodf(45)-prodb(45))+(1.0-0.0)*(prodf(46)-prodb(46))+(0.0-1.0)*(prodf(47)-prodb(47))+(1.0-0.0)*(prodf(48)-prodb(48))+(0.0-1.0)*(prodf(49)-prodb(49))+(0.0-1.0)*(prodf(60)-prodb(60))+(1.0-0.0)*(prodf(61)-prodb(61))+(0.0-1.0)*(prodf(80)-prodb(80))+(1.0-0.0)*(prodf(81)-prodb(81))+(2.0-0.0)*(prodf(82)-prodb(82))+(0.0-2.0)*(prodf(83)-prodb(83))+(1.0-0.0)*(prodf(86)-prodb(86))+(0.0-1.0)*(prodf(87)-prodb(87))+(0.0-1.0)*(prodf(90)-prodb(90))+(1.0-0.0)*(prodf(91)-prodb(91))+(1.0-0.0)*(prodf(93)-prodb(93))+(0.0-1.0)*(prodf(94)-prodb(94)))
  omegadot(16)=Wm_tab(16)*(+(1.0-0.0)*(prodf(50)-prodb(50))+(0.0-1.0)*(prodf(51)-prodb(51))+(1.0-0.0)*(prodf(62)-prodb(62))+(0.0-1.0)*(prodf(63)-prodb(63))+(1.0-0.0)*(prodf(69)-prodb(69))+(0.0-1.0)*(prodf(70)-prodb(70))+(1.0-0.0)*(prodf(71)-prodb(71))+(0.0-1.0)*(prodf(72)-prodb(72))+(1.0-0.0)*(prodf(73)-prodb(73))+(0.0-1.0)*(prodf(74)-prodb(74))+(0.0-1.0)*(prodf(75)-prodb(75))+(1.0-0.0)*(prodf(76)-prodb(76))+(0.0-1.0)*(prodf(77)-prodb(77))+(1.0-0.0)*(prodf(92)-prodb(92))+(0.0-1.0)*(prodf(103)-prodb(103))+(1.0-0.0)*(prodf(104)-prodb(104)))
  omegadot(17)=Wm_tab(17)*(+(1.0-0.0)*(prodf(64)-prodb(64))+(0.0-1.0)*(prodf(65)-prodb(65))+(0.0-1.0)*(prodf(66)-prodb(66))+(1.0-0.0)*(prodf(67)-prodb(67))+(0.0-1.0)*(prodf(68)-prodb(68))+(0.0-1.0)*(prodf(69)-prodb(69))+(1.0-0.0)*(prodf(70)-prodb(70))+(0.0-1.0)*(prodf(71)-prodb(71))+(1.0-0.0)*(prodf(72)-prodb(72))+(0.0-1.0)*(prodf(73)-prodb(73))+(1.0-0.0)*(prodf(74)-prodb(74))+(1.0-0.0)*(prodf(84)-prodb(84))+(0.0-1.0)*(prodf(85)-prodb(85))+(1.0-0.0)*(prodf(101)-prodb(101))+(0.0-1.0)*(prodf(102)-prodb(102)))
  omegadot(18)=Wm_tab(18)*(+(0.0-1.0)*(prodf(93)-prodb(93))+(1.0-0.0)*(prodf(94)-prodb(94)))
  omegadot(19)=Wm_tab(19)*(+(0.0-1.0)*(prodf(92)-prodb(92))+(1.0-0.0)*(prodf(95)-prodb(95))+(0.0-1.0)*(prodf(96)-prodb(96))+(0.0-1.0)*(prodf(97)-prodb(97))+(1.0-0.0)*(prodf(98)-prodb(98))+(1.0-0.0)*(prodf(99)-prodb(99))+(0.0-1.0)*(prodf(100)-prodb(100)))
  omegadot(20)=Wm_tab(20)*(+(0.0-1.0)*(prodf(95)-prodb(95))+(1.0-0.0)*(prodf(96)-prodb(96))+(0.0-1.0)*(prodf(99)-prodb(99))+(1.0-0.0)*(prodf(100)-prodb(100))+(0.0-1.0)*(prodf(101)-prodb(101))+(1.0-0.0)*(prodf(102)-prodb(102)))
  end subroutine ciottoli20

  function comp_ch_tabT(ireact,tab,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: tab(:,:), Tdiff
    ! Local
    real(8) :: a, b
    real(8) :: result
      
    a = tab(Tint(1),ireact)      ! int(T)   <- Tint(1)
    b = tab(Tint(2),ireact)      ! int(T)+1 <- Tint(2)
    result = a+(b-a)*Tdiff

  end function comp_ch_tabT


end module U_Lib_ChemMech
