module U_Lib_Chemistry_data
  implicit none

  integer :: nrc
  real(8), dimension(:,:), allocatable :: ni1_tab
  real(8), dimension(:,:), allocatable :: ni2_tab
  real(8), dimension(:,:), allocatable :: kf_tab
  real(8), dimension(:,:), allocatable :: kb_tab
  real(8), dimension(:,:), allocatable :: epsch_tab

contains

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

end module U_Lib_Chemistry_data