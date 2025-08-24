module U_Lib_Chemistry_data
  implicit none

  integer                              :: nrc
  integer, dimension(:), allocatable   :: rxn_type     ! 0 -> Arrhenius, 1 -> Troe
  ! Arrhenius
  integer                              :: nrc_arrh
  real(8), dimension(:,:), allocatable :: kf_tab
  real(8), dimension(:,:), allocatable :: kb_tab
  ! Troe
  integer                              :: nrc_troe
  real(8), dimension(:,:), allocatable :: kinf_tab
  real(8), dimension(:,:), allocatable :: k0_tab
  real(8), dimension(:,:), allocatable :: kc_tab
  real(8), dimension(:,:), allocatable :: Fcent_tab
  ! Arrhenius (for the general loop only)
  real(8), dimension(:,:), allocatable :: ni1_arrh_tab
  real(8), dimension(:,:), allocatable :: ni2_arrh_tab
  real(8), dimension(:,:), allocatable :: epsch_arrh_tab
  ! Troe (for the general loop only)
  real(8), dimension(:,:), allocatable :: ni1_troe_tab
  real(8), dimension(:,:), allocatable :: ni2_troe_tab
  real(8), dimension(:,:), allocatable :: epsch_troe_tab

contains

  pure function comp_ch_tabT(ireact,tab,Tint,Tdiff) result(result)
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

  pure function f_kc(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = kc_tab(Tint(1),ireact)
    b = kc_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_kc

  pure function f_kinf(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = kinf_tab(Tint(1),ireact)
    b = kinf_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_kinf

  pure function f_k0(ireact,Tint,Tdiff) result(result)
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
      
    a = k0_tab(Tint(1),ireact)
    b = k0_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff

  end function f_k0

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
    if (allocated(kf_tab)) deallocate(kf_tab)
    if (allocated(kb_tab)) deallocate(kb_tab)
    if (allocated(kinf_tab)) deallocate(kinf_tab)
    if (allocated(k0_tab)) deallocate(k0_tab)
    if (allocated(kc_tab)) deallocate(kc_tab)
    if (allocated(Fcent_tab)) deallocate(Fcent_tab)
    if (allocated(ni1_arrh_tab)) deallocate(ni1_arrh_tab)
    if (allocated(ni2_arrh_tab)) deallocate(ni2_arrh_tab)
    if (allocated(epsch_arrh_tab)) deallocate(epsch_arrh_tab)
    if (allocated(ni1_troe_tab)) deallocate(ni1_troe_tab)
    if (allocated(ni2_troe_tab)) deallocate(ni2_troe_tab)
    if (allocated(epsch_troe_tab)) deallocate(epsch_troe_tab)
  end subroutine free_chemistry_data

end module U_Lib_Chemistry_data