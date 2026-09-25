
module FLINT_Lib_Chemistry_wdot
  use iso_fortran_env, only: error_unit
  implicit none

  !> Fallback policy of Assign_Mechanism for a name that is not hooked: .false.
  !> (default) selects the general procedure and prints a WARNING on stdout and
  !> on the error unit; .true. refuses the name (error stop). The environment
  !> variable FLINT_STRICT_MECHANISM=1|true|yes|on turns it on as well.
  logical, public :: FLINT_strict_mechanism = .false.

  !> Concrete procedure pointing to one of the subroutine realizations
  procedure(chemsource_if), pointer, public :: chemistry_source

  !> Concrete procedure pointing to the analytical Jacobian of `chemistry_source`
  !> (null when the active mechanism has no analytical Jacobian yet; callers
  !> must fall back to finite-difference Jacobian in that case).
  procedure(chemjac_if), pointer, public :: chemistry_jacobian => null()

  !> Abstract interface relative to the finite-rate reactions source procedure
  abstract interface
  subroutine chemsource_if(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    implicit none
    integer :: is, T_i, Tint(2)
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in) :: temp
    real(8), intent(out) :: omegadot(ns)
    real(8) :: coi(ns+1), Tdiff
  end subroutine chemsource_if
  end interface

  !> Abstract interface for the analytical chemistry Jacobian.
  !>   dwdr(i,j) = d omegadot(i) / d roi(j)
  !>   dwdT(i)   = d omegadot(i) / d T
  abstract interface
  subroutine chemjac_if(roi,temp,dwdr,dwdT)
    use FLINT_Lib_Thermodynamic
    implicit none
    real(8), intent(in)  :: roi(ns), temp
    real(8), intent(out) :: dwdr(ns,ns)
    real(8), intent(out) :: dwdT(ns)
  end subroutine chemjac_if
  end interface

contains

  subroutine Assign_Mechanism(mad_world)
    use WD_mod
    use globH2_mod
    use JLRs_mod
    use smooke_mod
    use coria_mod
    use TSRCDF13_mod
    use TSRGP24_mod
    use TSRRich31_mod
    use ZK_mod
    use coronetti_mod
    use singh_mod
    use troyes_mod
    use ecker_mod
    use cross_mod
    use pelucchi_mod
    use ONERA7_mod
    use sandiego_mod
    use FFCMy_12_mod
    use Gerlinger9_mod
    use FLINT_Lib_Chemistry_contract, only: check_mechanism_contract
    implicit none
    character(*), intent(in) :: mad_world
    logical :: hooked, ok

    ! Default: no analytical Jacobian available. Each mechanism that has one
    ! overrides this below.
    chemistry_jacobian => null()
    hooked = .true.

    select case(mad_world)
    case('WD')
      chemistry_source => WD
    case('JLR-Nasuti')
      chemistry_source => JLR
    case('Frassoldati')
      chemistry_source => Frassoldati
    case('Smooke')
      chemistry_source => smooke
    case('CORIA-CNRS')
      chemistry_source => coria
    case('TSR-CDF-13')
      chemistry_source => TSRCDF13
    case('TSR-GP-24')
      chemistry_source => TSRGP24
    case('TSR-Rich-31')
      chemistry_source => TSRRich31
    case('ZK')
      chemistry_source => ZK
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
    case('Troyes')
      chemistry_source => troyes
    case('Ecker')
      chemistry_source => ecker
    case('Cross')
      chemistry_source => cross
    case('Pelucchi')
      chemistry_source => pelucchi
    case('WD-Andersen')
      chemistry_source => Andersen
    case('OSK')
      chemistry_source => OSK
    case('ONERA-7')
      chemistry_source   => ONERA_7
      chemistry_jacobian => ONERA_7_jac
    case('Frolov_nopressure')
      chemistry_source   => Frolov_nopressure
      chemistry_jacobian => Frolov_nopressure_jac
    case('SanDiego')
      chemistry_source => sandiego20161214
    case('FFCMy-12')
      chemistry_source => FFCMy_12
    case('Gerlinger-9')
      chemistry_source => Gerlinger9

    case default
      hooked = .false.
      write(*,*) "[WARNING] Explicit procedure for "//trim(mad_world)//" not found, defaulting to the general procedure"
      write(error_unit,'(A)') "[WARNING] FLINT Assign_Mechanism: explicit procedure for "//trim(mad_world)// &
        " not found, defaulting to the general procedure"
      if (strict_mechanism()) then
        write(*,'(A)') "[ERROR] FLINT Assign_Mechanism: mechanism "//trim(mad_world)// &
          " is not hooked and strict mode is on (FLINT_STRICT_MECHANISM / FLINT_strict_mechanism)"
        write(error_unit,'(A)') "[ERROR] FLINT Assign_Mechanism: mechanism "//trim(mad_world)// &
          " is not hooked and strict mode is on (FLINT_STRICT_MECHANISM / FLINT_strict_mechanism)"
        error stop 1
      endif
      chemistry_source => general
    end select

    ! Mechanism contract: a hooked name selects a compiled routine whose species
    ! slots and reaction tables are fixed; the loaded data must match them
    ! (see FLINT_Lib_Chemistry_contract for the rules). The tables must have been
    ! loaded (read_idealgas_thermo, read_chemistry) before this call.
    if (hooked) then
      call check_mechanism_contract(mad_world, ok)
      if (.not. ok) error stop '[ERROR] FLINT Assign_Mechanism: mechanism contract violated (see the two lists above)'
    endif

  end subroutine Assign_Mechanism


  !> Strict fallback policy: the module flag or the environment variable FLINT_STRICT_MECHANISM,
  !> read case-insensitively: 1/true/yes/on turn strict mode on, 0/false/no/off (or an empty value)
  !> leave the module flag as it is; any other value is reported with a WARNING on both units and
  !> ignored (a value such as 'Yes' or 'On' was silently ignored by an exact-case comparison).
  function strict_mechanism() result(strict)
    implicit none
    logical :: strict
    character(len=64) :: val
    integer :: istat, k, c
    strict = FLINT_strict_mechanism
    val = ''
    call get_environment_variable('FLINT_STRICT_MECHANISM', value=val, status=istat)
    if (istat == 1 .or. istat == 2) return      ! not set / no environment on this processor
    val = adjustl(val)
    do k = 1, len_trim(val)
      c = iachar(val(k:k))
      if (c >= iachar('A') .and. c <= iachar('Z')) val(k:k) = achar(c + 32)
    enddo
    if (istat == 0) then
      select case (trim(val))
      case ('1', 'true', 'yes', 'on')
        strict = .true.
        return
      case ('', '0', 'false', 'no', 'off')
        return
      end select
    endif
    write(*,'(A)') "[WARNING] FLINT Assign_Mechanism: FLINT_STRICT_MECHANISM='"//trim(val)// &
      "' is not one of 1/true/yes/on or 0/false/no/off: ignored"
    write(error_unit,'(A)') "[WARNING] FLINT Assign_Mechanism: FLINT_STRICT_MECHANISM='"//trim(val)// &
      "' is not one of 1/true/yes/on or 0/false/no/off: ignored"
  end function strict_mechanism

  ! General mechanism
  subroutine general(roi,temp,omegadot)
    use FLINT_Lib_Thermodynamic
    use FLINT_Lib_Chemistry_data
    use FLINT_Lib_Chemistry_falloff
    implicit none
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in) :: temp 
    real(8), intent(out) :: omegadot(ns)
    ! Local
    integer :: is, T_i, Tint(2)
    real(8) :: coi(ns), Tdiff
    real(8) :: prod_fwd, prod_rev, deltani, k(2)
    real(8) :: rate_fwd, rate_rev, net_rate, coM
    integer :: ir

    do is = 1, ns
      coi(is) = roi(is)/Wm_tab(is)  ! kmol/m^3
    enddo

    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    ! Initialize mass source terms
    omegadot(:) = 0.d0

    ! Loop over Arrhenius reactions
    do ir = 1, nrc_arrh

      ! Compute third-body effective concentration
      if (ni1_arrh_tab(ns+1, ir)>0d0) then
        coM = 0.d0
        do is = 1, ns
          coM = coM + coi(is) * epsch_arrh_tab(is, ir)
        enddo
      else 
        coM = 1.d0
      endif

      ! Compute forward and reverse rate-of-progress: mass-action law with the real stoichiometric
      ! coefficients (reactants forward, products reverse), as in Cantera; the explicit orders of
      ! the optional block replace the reactant coefficients of the forward rate
      prod_fwd = 1.0d0
      prod_rev = 1.0d0
      do is = 1, ns
        if (have_orders) then
          ! explicit reaction orders (optional block of chemistry-info.txt, this contract)
          if (ord_arrh_tab(is, ir) /= 0) prod_fwd = prod_fwd * pow_order(coi(is), ord_arrh_tab(is, ir))
        else
          if (ni1_arrh_tab(is, ir) /= 0) prod_fwd = prod_fwd * pow_order(coi(is), ni1_arrh_tab(is, ir))
        endif
        if (ni2_arrh_tab(is, ir) /= 0) prod_rev = prod_rev * pow_order(coi(is), ni2_arrh_tab(is, ir))
      enddo
      rate_fwd = f_kf(ir,Tint,Tdiff) * prod_fwd * coM
      rate_rev = f_kb(ir,Tint,Tdiff) * prod_rev * coM
      net_rate = rate_fwd - rate_rev

      ! Sum up net production rate for each species
      do is = 1, ns
        deltani = ni2_arrh_tab(is, ir) - ni1_arrh_tab(is, ir)
        omegadot(is) = omegadot(is) + Wm_tab(is) * deltani * net_rate
      enddo

    enddo

    ! Loop over falloff-Troe reactions
    do ir = 1, nrc_troe

      ! Compute third-body effective concentration
      if (ni1_troe_tab(ns+1, ir)>0d0) then
        coM = 0.d0
        do is = 1, ns
          coM = coM + coi(is) * epsch_troe_tab(is, ir)
        enddo
      else 
        coM = 1.d0
      endif

      ! Compute forward and reverse rate-of-progress
      prod_fwd = 1.0d0
      prod_rev = 1.0d0
      do is = 1, ns
        if (ni1_troe_tab(is, ir) /= 0) prod_fwd = prod_fwd * coi(is)**ni1_troe_tab(is, ir)
        if (ni2_troe_tab(is, ir) /= 0) prod_rev = prod_rev * coi(is)**ni2_troe_tab(is, ir)
      enddo
      k = f_k_troe(ir,Tint,Tdiff,coM)
      rate_fwd = k(1) * prod_fwd
      rate_rev = k(2) * prod_rev
      net_rate = rate_fwd - rate_rev

      ! Sum up net production rate for each species
      do is = 1, ns
        deltani = ni2_troe_tab(is, ir) - ni1_troe_tab(is, ir)
        omegadot(is) = omegadot(is) + Wm_tab(is) * deltani * net_rate
      enddo

    enddo

    ! Loop over falloff-Lindemann reactions
    do ir = 1, nrc_lindemann

      ! Compute third-body effective concentration
      if (ni1_lind_tab(ns+1, ir)>0d0) then
        coM = 0.d0
        do is = 1, ns
          coM = coM + coi(is) * epsch_lind_tab(is, ir)
        enddo
      else 
        coM = 1.d0
      endif

      ! Compute forward and reverse rate-of-progress
      prod_fwd = 1.0d0
      prod_rev = 1.0d0
      do is = 1, ns
        if (ni1_lind_tab(is, ir) /= 0) prod_fwd = prod_fwd * coi(is)**ni1_lind_tab(is, ir)
        if (ni2_lind_tab(is, ir) /= 0) prod_rev = prod_rev * coi(is)**ni2_lind_tab(is, ir)
      enddo
      k = f_k_lindemann(ir,Tint,Tdiff,coM)
      rate_fwd = k(1) * prod_fwd
      rate_rev = k(2) * prod_rev
      net_rate = rate_fwd - rate_rev

      ! Sum up net production rate for each species
      do is = 1, ns
        deltani = ni2_lind_tab(is, ir) - ni1_lind_tab(is, ir)
        omegadot(is) = omegadot(is) + Wm_tab(is) * deltani * net_rate
      enddo

    enddo

  end subroutine general


  !> x**n for the small integer exponents that reaction orders actually take.
  !>
  !> The `**` operator with a run-time integer exponent lowers to a call into
  !> libgcc's __powidf2. Orders of 1, 2 and 3 cover every
  !> reaction in the shipped mechanisms and are expanded inline here; anything
  !> else falls back to the intrinsic, so results are unchanged.
  !> c**o for a reaction order or a stoichiometric coefficient o: an integer-valued o goes through ipow
  !> (bit-identical to the integer powers used before), a non-integer one through the real power of
  !> max(c, 0) (Cantera also gives 0 at c <= 0 for a non-integer exponent); a negative order
  !> at zero concentration gives 0, the convention of Cantera (its forward rate of progress is 0
  !> there, verified with Cantera 3.0.1) so that the same mechanism gives the same rates.
  pure function pow_order(c, o) result(y)
    implicit none
    real(8), intent(in) :: c, o
    real(8) :: y
    if (c <= 0d0 .and. o < 0d0) then
      y = 0d0
    else if (o == dble(nint(o))) then
      y = ipow(c, nint(o))
    else
      y = max(c, 0d0)**o
    endif
  end function pow_order

  pure function ipow(x, n) result(y)
    implicit none
    real(8), intent(in) :: x
    integer, intent(in) :: n
    real(8) :: y

    select case (n)
    case (1)
      y = x
    case (2)
      y = x*x
    case (3)
      y = x*x*x
    case (0)
      y = 1.0d0
    case default
      y = x**n
    end select

  end function ipow

end module FLINT_Lib_Chemistry_wdot
