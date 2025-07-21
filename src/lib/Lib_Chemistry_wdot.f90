
module U_Lib_Chemistry_wdot
  implicit none

  !> Concrete procedure pointing to one of the subroutine realizations
  procedure(chemsource_if), pointer, public :: chemistry_source

  !> Abstract interface relative to the finite-rate reactions source procedure
  abstract interface
  subroutine chemsource_if(roi,temp,omegadot,rotot)
    use U_Lib_Thermodynamic
    use U_Lib_Chemistry_data
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
    use WD_mod
    use globH2_mod
    use JLRs_mod
    use coronetti_mod
    use singh_mod
    use troyes_mod
    use ciottoli20_mod
    use ecker_mod
    use cross_mod
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
    case('Ciottoli20')
      chemistry_source => ciottoli20
    case('Troyes')
      chemistry_source => troyes
    case('Ecker')
      chemistry_source => ecker
    case('Cross')
      chemistry_source => cross
    case default
      write(*,*) "[WARNING] Explicit procedure for "//trim(mad_world)//" not found, defaulting to the general procedure"
      chemistry_source => general
    end select

  end subroutine Assign_Mechanism


  ! General mechanism
  subroutine general(roi,temp,omegadot,rotot)
    use U_Lib_Thermodynamic
    use U_Lib_Chemistry_data
    implicit none
    real(8), intent(in) :: roi(nsc)
    real(8), intent(in) :: temp 
    real(8), intent(out) :: omegadot(nsc)
    real(8), intent(in) :: rotot
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(nsc+1), Tdiff
    real(8) :: prod_fwd, prod_rev, deltani
    real(8) :: rate_fwd, rate_rev, net_rate, coM
    integer :: s, ir

    do is = 1, nsc
      coi(is)=roi(is)/Wm_tab(is)  ! kmol/m^3
      ! Loop done in order to avoid numerical issues in omegadot evaluation
      ! Very low coi could produce finite prods
      if (coi(is).lt.0d0) coi(is) = 0.0d0
    enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! Initialize mass source terms
    omegadot(:) = 0.d0

    ! Loop over all reactions
    do ir = 1, nrc

      ! Compute third-body effective concentration
      if (ni1_tab(nsc+1, ir)>0d0) then
        coM = 0.d0
        do is = 1, nsc
          coM = coM + coi(is) * epsch_tab(is, ir)
        enddo
      else 
        coM = 1.d0
      endif

      ! Compute forward and reverse rate-of-progress
      prod_fwd = 1.0d0
      prod_rev = 1.0d0
      do is = 1, nsc
        if (ni1_tab(is, ir) /= 0) prod_fwd = prod_fwd * coi(is)**ni1_tab(is, ir)
        if (ni2_tab(is, ir) /= 0) prod_rev = prod_rev * coi(is)**ni2_tab(is, ir)
      enddo
      rate_fwd = comp_ch_tabT(ir,kf_tab,Tint,Tdiff) * prod_fwd * coM
      rate_rev = comp_ch_tabT(ir,kb_tab,Tint,Tdiff) * prod_rev * coM
      net_rate = rate_fwd - rate_rev

      ! Sum up net production rate for each species
      do is = 1, nsc
        deltani = ni2_tab(is, ir) - ni1_tab(is, ir)
        omegadot(is) = omegadot(is) + Wm_tab(is) * deltani * net_rate
      enddo

    enddo

  end subroutine general

end module U_Lib_Chemistry_wdot
