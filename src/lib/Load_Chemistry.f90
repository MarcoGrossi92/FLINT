module U_IO_chemistry
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none

contains

  function read_chemistry_file( folder, mech_name ) result(ios)
    use Lib_ORION_Data
    use Lib_Tecplot
    use U_Lib_Chemistry_data
    use U_Lib_Thermodynamic, only: nsc, U_phase_prefix
    implicit none
    character(len=*), intent(in), optional :: folder
    character(len=*), intent(out), optional :: mech_name
    integer :: idum, unitfile
    type(ORION_data)  :: orion
    integer :: i, j, ios
    integer :: Ti1, Ti2
    character(len=16):: chardum

    if (present(folder)) then
      ios = tec_read_points_multivars(orion,2,trim(folder)//'/'//trim(U_phase_prefix)//'chemistry-rate.dat')
    else
      ios = tec_read_points_multivars(orion,2,trim('INPUT/')//trim(U_phase_prefix)//'chemistry-rate.dat')
    endif
    if (ios/=0) then
      write(*,*) '[WARNING] chemistry files not found'
      return
    endif
    Ti1 = nint(orion%block(1)%mesh(1,1,1,1))
    Ti2 = Ti1 + orion%block(1)%Ni - 1
    nrc = size(orion%block)
    allocate(kf_tab(Ti1:Ti2, 1:nrc))
    allocate(kb_tab(Ti1:Ti2, 1:nrc))
    do i = 1, nrc
      kf_tab(Ti1:Ti2,i) = orion%block(i)%vars(1,1:orion%block(1)%Ni,1,1)
      kb_tab(Ti1:Ti2,i) = orion%block(i)%vars(2,1:orion%block(1)%Ni,1,1)
    enddo

    if (present(folder)) then
      open(newunit=unitfile,file=trim(folder)//'/'//trim(U_phase_prefix)//'chemistry-stoich.txt',form='formatted',status='old',action='read')
    else
      open(newunit=unitfile,file='INPUT/'//trim(U_phase_prefix)//'chemistry-stoich.txt',form='formatted',status='old',action='read')
    endif
    read(unitfile,*) mech_name
    allocate(ni1_tab(1:nsc+1,1:nrc))
    allocate(ni2_tab(1:nsc+1,1:nrc))
    allocate(epsch_tab(1:nsc+1,1:nrc))
    do j = 1, nrc; do i = 1, nsc+1
        read(unitfile,*)idum,chardum,ni1_tab(i,j),ni2_tab(i,j),epsch_tab(i,j)
    enddo; enddo 
    close(unitfile)

  end function read_chemistry_file


end module U_IO_chemistry