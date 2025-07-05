module U_IO_chemistry
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none

contains

  function read_chemistry_file( file1, file2, mech_name ) result(ios)
    use Lib_ORION_Data
    use Lib_Tecplot
    use U_Lib_ChemMech, only: kf_tab, kb_tab, ni1_tab, ni2_tab, epsch_tab, nrc
    use U_Lib_Thermodynamic, only: nsc
    implicit none
    character(len=*), intent(in), optional :: file1, file2
    character(len=*), intent(out), optional :: mech_name
    integer :: idum, unitfile
    type(ORION_data)  :: orion
    integer :: i, j, ios
    integer :: Ti1, Ti2
    character(len=16):: chardum

    if (present(file1)) then
      ios = tec_read_points_multivars(orion,2,trim(file1))
    else
      ios = tec_read_points_multivars(orion,2,trim('INPUT/chemistry-rate.dat'))
    endif
    if (ios/=0) return
    Ti1 = nint(orion%block(1)%mesh(1,1,1,1))
    Ti2 = Ti1 + orion%block(1)%Ni - 1
    nrc = size(orion%block)
    allocate(kf_tab(Ti1:Ti2, 1:nrc))
    allocate(kb_tab(Ti1:Ti2, 1:nrc))
    do i = 1, nrc
      kf_tab(Ti1:Ti2,i) = orion%block(i)%vars(1,1:orion%block(1)%Ni,1,1)
      kb_tab(Ti1:Ti2,i) = orion%block(i)%vars(2,1:orion%block(1)%Ni,1,1)
    enddo

    if (present(file2)) then
      open(newunit=unitfile,file=file2,form='formatted',status='old',action='read')
    else
      open(newunit=unitfile,file='INPUT/chemistry-stoich.txt',form='formatted',status='old',action='read')
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