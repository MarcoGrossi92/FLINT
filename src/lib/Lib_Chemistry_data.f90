module FLINT_Lib_Chemistry_data
  implicit none

  integer                              :: nrc
  integer, dimension(:), allocatable   :: rxn_type ! 0 -> Arrhenius, 1 -> Troe, 2 -> Lindemann
  !> Temperature bounds [K] of the loaded rate tables: every table below is
  !> allocated as tab(T_tab_min:T_tab_max, :) and row T holds the rate at T kelvin
  !> (the first row of chemistry-*.dat gives T_tab_min, e.g. 1 K for the tables of
  !> database/WD). Set by read_chemistry, reset by free_chemistry_data.
  integer                              :: T_tab_min = 1
  integer                              :: T_tab_max = 0
  ! Arrhenius
  integer                              :: nrc_arrh
  real(8), dimension(:,:), allocatable :: kf_tab
  real(8), dimension(:,:), allocatable :: kb_tab
  ! Falloff-Troe
  integer                              :: nrc_troe
  real(8), dimension(:,:), allocatable :: kinf_troe_tab
  real(8), dimension(:,:), allocatable :: k0_troe_tab
  real(8), dimension(:,:), allocatable :: kc_troe_tab
  real(8), dimension(:,:), allocatable :: Fcent_tab
  ! Falloff-Lindemann
  integer                              :: nrc_lindemann
  real(8), dimension(:,:), allocatable :: kinf_lind_tab
  real(8), dimension(:,:), allocatable :: k0_lind_tab
  real(8), dimension(:,:), allocatable :: kc_lind_tab
  ! Arrhenius (for the general loop only)
  !> Forward reaction orders of the Arrhenius reactions from the optional 'Reaction orders' block
  !> of chemistry-info.txt: ord_arrh_tab(is, ir) = order of species is in the
  !> Arrhenius reaction ir (the explicit order where given, the stoichiometric reactant
  !> coefficient elsewhere). have_orders = .false. (no block): the general loop raises the
  !> concentrations to the real stoichiometric reactant coefficients, the same exponents as a block with
  !> no rows. A file without the block comes from an older table writer: when the general procedure is
  !> selected (general_selected, set by Assign_Mechanism) a WARNING says so once (orders_block_warned).
  logical                              :: have_orders = .false.
  logical                              :: general_selected = .false.
  logical                              :: orders_block_warned = .false.
  real(8), dimension(:,:), allocatable :: ord_arrh_tab
  real(8), dimension(:,:), allocatable :: ni1_arrh_tab
  real(8), dimension(:,:), allocatable :: ni2_arrh_tab
  real(8), dimension(:,:), allocatable :: epsch_arrh_tab
  ! Falloff-Troe (for the general loop only)
  real(8), dimension(:,:), allocatable :: ni1_troe_tab
  real(8), dimension(:,:), allocatable :: ni2_troe_tab
  real(8), dimension(:,:), allocatable :: epsch_troe_tab
  ! Falloff-Lindemann (for the general loop only)
  real(8), dimension(:,:), allocatable :: ni1_lind_tab
  real(8), dimension(:,:), allocatable :: ni2_lind_tab
  real(8), dimension(:,:), allocatable :: epsch_lind_tab

contains

  !> NaN (or infinity) test on the bit pattern (exponent bits all set): unlike isnan(x) or x /= x it
  !> is not folded away by -ffast-math / -ffinite-math-only and it raises no floating-point exception.
  pure logical function nan_bits(x)
    real(8), intent(in) :: x
    nan_bits = iand(shiftr(transfer(x, 0_8), 52), 2047_8) == 2047_8
  end function nan_bits

  !> Linear interpolation of the rate of reaction `ireact` in the rate table `tab`
  !> between the rows Tint(1) = int(T) and Tint(2) = int(T)+1. `tab` is one of the
  !> module rate tables (kf_tab, kb_tab, the falloff tables) or any table on the
  !> same temperature grid: the dummy argument takes the lower bound T_tab_min of
  !> the loaded tables, so row T is the rate at T kelvin whatever the first
  !> temperature of the table. (Up to FLINT 2223136 the dummy was tab(:,:), which
  !> renumbers the rows from 1: for a table starting at T0 > 1 K that version
  !> returned the rate at T + T0 - 1; for tables starting at 1 K both versions
  !> return the same numbers.) The routines of FLINT read
  !> the tables through f_kf/f_kb and the other accessors below.
  pure function comp_ch_tabT(ireact,tab,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: tab(T_tab_min:,:), Tdiff
    ! Local
    real(8) :: a, b
    real(8) :: result

    a = tab(Tint(1),ireact)      ! int(T)   <- Tint(1)
    b = tab(Tint(2),ireact)      ! int(T)+1 <- Tint(2)
    result = a+(b-a)*Tdiff

  end function comp_ch_tabT

  !> Linear interpolation of the forward (f_kf) / backward (f_kb) rate of the
  !> Arrhenius reaction `ireact` between the rows Tint(1) = int(T) and
  !> Tint(2) = int(T)+1 of the module tables, whose row T is the rate at T kelvin
  !> whatever the first temperature of the table.

  pure function f_kf(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = kf_tab(Tint(1),ireact)
    b = kf_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_kf

  pure function f_kb(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = kb_tab(Tint(1),ireact)
    b = kb_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_kb

  pure function f_kc_lind(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = kc_lind_tab(Tint(1),ireact)
    b = kc_lind_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_kc_lind

  pure function f_kinf_lind(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = kinf_lind_tab(Tint(1),ireact)
    b = kinf_lind_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_kinf_lind

  pure function f_k0_lind(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = k0_lind_tab(Tint(1),ireact)
    b = k0_lind_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_k0_lind

  pure function f_kc_troe(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = kc_troe_tab(Tint(1),ireact)
    b = kc_troe_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_kc_troe

  pure function f_kinf_troe(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = kinf_troe_tab(Tint(1),ireact)
    b = kinf_troe_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_kinf_troe

  pure function f_k0_troe(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = k0_troe_tab(Tint(1),ireact)
    b = k0_troe_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_k0_troe

  pure function f_Fcent(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = Fcent_tab(Tint(1),ireact)
    b = Fcent_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_Fcent

  ! Free all chemistry data
  subroutine free_chemistry_data()
    implicit none
    if (allocated(rxn_type)) deallocate(rxn_type)
    T_tab_min = 1; T_tab_max = 0
    if (allocated(kf_tab)) deallocate(kf_tab)
    if (allocated(kb_tab)) deallocate(kb_tab)
    if (allocated(kinf_lind_tab)) deallocate(kinf_lind_tab)
    if (allocated(k0_lind_tab)) deallocate(k0_lind_tab)
    if (allocated(kc_lind_tab)) deallocate(kc_lind_tab)
    if (allocated(kinf_troe_tab)) deallocate(kinf_troe_tab)
    if (allocated(k0_troe_tab)) deallocate(k0_troe_tab)
    if (allocated(kc_troe_tab)) deallocate(kc_troe_tab)
    if (allocated(Fcent_tab)) deallocate(Fcent_tab)
    have_orders = .false.
    orders_block_warned = .false.
    if (allocated(ord_arrh_tab)) deallocate(ord_arrh_tab)
    if (allocated(ni1_arrh_tab)) deallocate(ni1_arrh_tab)
    if (allocated(ni2_arrh_tab)) deallocate(ni2_arrh_tab)
    if (allocated(epsch_arrh_tab)) deallocate(epsch_arrh_tab)
    if (allocated(ni1_lind_tab)) deallocate(ni1_lind_tab)
    if (allocated(ni2_lind_tab)) deallocate(ni2_lind_tab)
    if (allocated(epsch_lind_tab)) deallocate(epsch_lind_tab)
    if (allocated(ni1_troe_tab)) deallocate(ni1_troe_tab)
    if (allocated(ni2_troe_tab)) deallocate(ni2_troe_tab)
    if (allocated(epsch_troe_tab)) deallocate(epsch_troe_tab)
  end subroutine free_chemistry_data

  !> WARNING (once) of the general procedure for a chemistry-info.txt without the 'Reaction orders' block
  subroutine warn_no_orders_block()
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none
    character(len=*), parameter :: msg = "[WARNING] FLINT: chemistry-info.txt has no 'Reaction orders' block"// &
      " (the file comes from an older table writer): the general procedure takes the stoichiometric reactant"// &
      " coefficients as orders; regenerate the chemistry tables so that they carry the orders of the mechanism"
    if (orders_block_warned) return
    orders_block_warned = .true.
    write(*,'(A)') msg
    write(error_unit,'(A)') msg
  end subroutine warn_no_orders_block

end module FLINT_Lib_Chemistry_data