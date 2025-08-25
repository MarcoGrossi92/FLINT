module FLINT_IO_chemistry
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none

contains

  function read_chemistry_file( folder, mech_name ) result(ios)
    use Lib_ORION_Data
    use Lib_Tecplot
    use FLINT_Lib_Chemistry_data
    use FLINT_Lib_Thermodynamic, only: nsc, FLINT_phase_prefix
    implicit none
    character(len=*), intent(in), optional :: folder
    character(len=*), intent(out), optional :: mech_name
    integer :: idum, unitfile
    type(ORION_data)  :: orion
    integer :: i, j, ios, j0, j1
    integer :: Ti1, Ti2
    character(len=16):: chardum

    nrc_arrh = 0
    nrc_troe = 0

    !! Info
    if (present(folder)) then
      open(newunit=unitfile,file=trim(folder)//'/'//trim(FLINT_phase_prefix)//'chemistry-info.txt',form='formatted',status='old',action='read',iostat=ios)
    else
      open(newunit=unitfile,file='INPUT/'//trim(FLINT_phase_prefix)//'chemistry-info.txt',form='formatted',status='old',action='read',iostat=ios)
    endif
    if (ios/=0) then
      write(*,*) '[WARNING] chemistry-info file not found'
      return
    endif
    ! Read mechanism name
    read(unitfile,*) mech_name
    read(unitfile,*)
    read(unitfile,'(A17,I4)') chardum, nrc
    do j = 1, 4; read(unitfile,*); enddo
    ! Read reaction type
    allocate(rxn_type(1:nrc))
    do j = 1, nrc
      read(unitfile,*) idum, chardum
      if (index(trim(chardum),'Troe')>0) then
        nrc_troe = nrc_troe + 1
        rxn_type(j) = 1
      else
        nrc_arrh = nrc_arrh + 1
        rxn_type(j) = 0
      endif
    enddo
    ! Read info for general loop
    allocate(ni1_arrh_tab(1:nsc+1,1:nrc_arrh))
    allocate(ni2_arrh_tab(1:nsc+1,1:nrc_arrh))
    allocate(epsch_arrh_tab(1:nsc+1,1:nrc_arrh))
    if (nrc_troe>0) then
      allocate(ni1_troe_tab(1:nsc+1,1:nrc_troe))
      allocate(ni2_troe_tab(1:nsc+1,1:nrc_troe))
      allocate(epsch_troe_tab(1:nsc+1,1:nrc_troe))
    endif
    ! Read reaction info
    read(unitfile,*)
    read(unitfile,*)
    j0 = 0; j1 = 0
    do j = 1, nrc
      if (rxn_type(j)==0) then
        j0 = j0+1
        do i = 1, nsc+1
          read(unitfile,*)idum,chardum,ni1_arrh_tab(i,j0),ni2_arrh_tab(i,j0),epsch_arrh_tab(i,j0)
        enddo
      elseif (rxn_type(j)==1) then
        j1 = j1+1
        do i = 1, nsc+1
          read(unitfile,*)idum,chardum,ni1_troe_tab(i,j1),ni2_troe_tab(i,j1),epsch_troe_tab(i,j1)
        enddo
      endif
    enddo 
    close(unitfile)

    !! Rate Arrhenius
    if (present(folder)) then
      ios = tec_read_points_multivars(orion,2,trim(folder)//'/'//trim(FLINT_phase_prefix)//'chemistry-Arrhenius.dat')
    else
      ios = tec_read_points_multivars(orion,2,trim('INPUT/')//trim(FLINT_phase_prefix)//'chemistry-Arrhenius.dat')
    endif
    if (ios/=0) then
      write(*,*) '[ERROR] chemistry-Arrhenius file not found'
      stop
    endif
    Ti1 = nint(orion%block(1)%mesh(1,1,1,1))
    Ti2 = Ti1 + orion%block(1)%Ni - 1
    allocate(kf_tab(Ti1:Ti2, 1:nrc_arrh))
    allocate(kb_tab(Ti1:Ti2, 1:nrc_arrh))
    do i = 1, nrc_arrh
      kf_tab(Ti1:Ti2,i) = orion%block(i)%vars(1,1:orion%block(1)%Ni,1,1)
      kb_tab(Ti1:Ti2,i) = orion%block(i)%vars(2,1:orion%block(1)%Ni,1,1)
    enddo

    !! Rate Troe
    if (nrc_troe/=0) then
      deallocate(orion%block)
      if (present(folder)) then
        ios = tec_read_points_multivars(orion,4,trim(folder)//'/'//trim(FLINT_phase_prefix)//'chemistry-Troe.dat')
      else
        ios = tec_read_points_multivars(orion,4,trim('INPUT/')//trim(FLINT_phase_prefix)//'chemistry-Troe.dat')
      endif
      if (ios/=0) then
        write(*,*) '[ERROR] chemistry-Troe file not found'
        stop
      endif
      allocate(kinf_tab(Ti1:Ti2, 1:nrc_troe))
      allocate(k0_tab(Ti1:Ti2, 1:nrc_troe))
      allocate(kc_tab(Ti1:Ti2, 1:nrc_troe))
      allocate(Fcent_tab(Ti1:Ti2, 1:nrc_troe))
      do i = 1, nrc_troe
        kinf_tab(Ti1:Ti2,i)  = orion%block(i)%vars(1,1:orion%block(1)%Ni,1,1)
        k0_tab(Ti1:Ti2,i)    = orion%block(i)%vars(2,1:orion%block(1)%Ni,1,1)
        kc_tab(Ti1:Ti2,i)    = orion%block(i)%vars(3,1:orion%block(1)%Ni,1,1)
        Fcent_tab(Ti1:Ti2,i) = orion%block(i)%vars(4,1:orion%block(1)%Ni,1,1)
      enddo
    endif

  end function read_chemistry_file


end module FLINT_IO_chemistry