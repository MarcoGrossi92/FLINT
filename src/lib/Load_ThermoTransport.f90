! read_idealgas_thermo ios map:
! 0 -> no error
! 1 -> phase.txt not found
! 2 -> phase.txt found but error reading
! 3 -> thermo file not found
! 4 -> thermo file found but error reading
! 5 -> composition file not found (non-fatal, treated as 0 by caller)
! 6 -> composition file found but error reading
!
! read_idealgas_composition ios map:
! 0 -> no error
! 1 -> file not found
! 2 -> file found but error reading
!
! read_idealgas_transport ios map:
! 0 -> no error
! 1 -> transport file not found
! 2 -> error reading transport data

module FLINT_Load_ThermoTransport
  use FLINT_Lib_Thermodynamic
  implicit none

contains

  function read_idealgas_thermo(folder) result(ios)
    use strings, only: parse
    use Lib_Tecplot
    use Lib_ORION_data
    implicit none
    character(len=*), intent(in), optional :: folder
    ! Local
    integer           :: ios, i, unitfile, start, dummy1, dummy23
    logical           :: exists_dat, exists_szplt
    character(256)    :: wholestring, args(2)
    character(512)    :: wmfile, thermofile(2)
    type(ORION_data)  :: orion

    if (present(folder)) then
      wmfile = trim(folder)//'/'//trim(FLINT_phase_prefix)//'phase.txt'
      thermofile(1) = trim(folder)//'/'//trim(FLINT_phase_prefix)//'thermo.dat'
      thermofile(2) = trim(folder)//'/'//trim(FLINT_phase_prefix)//'thermo.szplt'
    else
      wmfile = 'INPUT/'//trim(FLINT_phase_prefix)//'phase.txt'
      thermofile(1) = 'INPUT/'//trim(FLINT_phase_prefix)//'thermo.dat'
      thermofile(2) = 'INPUT/'//trim(FLINT_phase_prefix)//'thermo.szplt'
    endif

    ! File 1: phase
    open(newunit=unitFile,file=trim(wmfile),status='old',iostat=ios)
    if (ios/=0) then
      ios = 1
      return
    endif
    ios = 0; ns = -1
    read(unitFile,*,iostat=ios)
    if (ios/=0) then; close(unitFile); ios = 2; return; endif
    ios = 0
    do while(ios==0)
      read(unitFile,'(A)',iostat=ios) wholestring
      ns = ns + 1
    enddo
    allocate(wm_tab(1:ns))
    allocate(Ri_tab(1:ns))
    allocate(species_names(1:ns))
    rewind(unitFile)
    read(unitFile,*)
    do i = 1, ns
      read(unitFile,'(A)') wholestring
      call parse(wholestring,' ',args)
      species_names(i) = trim(args(1))
      read(args(2),*) wm_tab(i)
    end do
    close(unitFile)

    Ri_tab = Runiv/wm_tab

    ! File 2: thermo
    inquire(file=trim(thermofile(1)), exist=exists_dat)
    inquire(file=trim(thermofile(2)), exist=exists_szplt)
    if (.not. exists_dat .and. .not. exists_szplt) then
      ios = 3
      return
    endif
    ios = -1
    if (exists_dat)   ios = tec_read_points_multivars(orion,4,trim(thermofile(1)))
    if (ios/=0 .and. exists_szplt) &
      ios = tec_read_structured_multiblock(orion=orion, filename=trim(thermofile(2)))
    if (ios/=0) then
      ios = 4
      return
    endif
    dummy1  = lbound(orion%block(1)%mesh, dim=2)
    dummy23 = lbound(orion%block(1)%mesh, dim=3)
    Tmin = nint(orion%block(1)%mesh(1,dummy1,dummy23,dummy23))
    Tmax = Tmin + ubound(orion%block(1)%mesh, dim=2) - dummy1
    start = Tmin
    if (Tmin==1) Tmin = 0
    allocate(h_tab(Tmin:Tmax, 1:ns))
    allocate(s_tab, cp_tab, dcpi_tab, mold=h_tab)
    dummy23 = lbound(orion%block(1)%vars, dim=3)
    do i = 1, ns
      cp_tab(start:Tmax,i) = orion%block(i)%vars(1,:,dummy23,dummy23)
      h_tab(start:Tmax,i) = orion%block(i)%vars(2,:,dummy23,dummy23)
      s_tab(start:Tmax,i) = orion%block(i)%vars(3,:,dummy23,dummy23)
      dcpi_tab(start:Tmax,i) = orion%block(i)%vars(4,:,dummy23,dummy23)
    enddo
    if (Tmin == 0) then
      cp_tab(0,:) = cp_tab(1,:)
      h_tab(0,:) = cp_tab(1,:)*0.d0
      s_tab(0,:) = s_tab(1,:)
      dcpi_tab(0,:) = dcpi_tab(1,:)
    endif

    ios = read_idealgas_composition(folder)
    if (ios == 1) ios = 5
    if (ios == 2) ios = 6

  end function read_idealgas_thermo


  function read_idealgas_composition(folder) result(ios)
    use strings, only: parse
    implicit none
    character(len=*), intent(in), optional :: folder
    ! Local
    integer           :: ios, i, unitfile
    character(512)    :: file

    if (present(folder)) then
      file = trim(folder)//'/'//trim(FLINT_phase_prefix)//'composition.txt'
    else
      file = 'INPUT/'//trim(FLINT_phase_prefix)//'composition.txt'
    endif

    ! File 1: phase
    open(newunit=unitFile,file=trim(file),status='old',iostat=ios)
    if (ios/=0) then
      ios = 1
      return
    endif
    read(unitfile,*,iostat=ios) ne
    if (ios/=0) then; close(unitFile); ios = 2; return; endif
    allocate(elements_names(1:ne))
    allocate(species_composition(1:ne,1:ns))
    do i = 1, ne
      read(unitFile,'(A)',iostat=ios) elements_names(i)
    enddo
    read(unitFile,*,iostat=ios) ! blank
    read(unitFile,*,iostat=ios) ! composition word
    do i = 1, ns
      read(unitFile,*,iostat=ios) species_composition(:,i)
      if (ios/=0) then; close(unitFile); ios = 2; return; endif
    end do
    close(unitFile)

  end function read_idealgas_composition


  function read_idealgas_transport(folder) result(ios)
    use strings, only: parse
    use Lib_Tecplot
    use Lib_ORION_data
    implicit none
    character(len=*), intent(in), optional :: folder
    ! Local
    integer           :: ios, i, start, dummy1, dummy23
    logical           :: exists, attempted
    character(512)    :: transfile(2)
    type(ORION_data)  :: orion

    ios = 1
    attempted = .false.

    if (present(folder)) then
      transfile(1) = trim(folder)//'/'//trim(FLINT_phase_prefix)//'transport.dat'
      transfile(2) = trim(folder)//'/'//trim(FLINT_phase_prefix)//'transport.szplt'
    else
      transfile(1) = 'INPUT/'//trim(FLINT_phase_prefix)//'transport.dat'
      transfile(2) = 'INPUT/'//trim(FLINT_phase_prefix)//'transport.szplt'
    endif

    inquire(file=trim(transfile(1)),exist=exists)
    if (exists) then
      attempted = .true.
      ios = tec_read_points_multivars(orion,2,trim(transfile(1)))
    endif
    inquire(file=trim(transfile(2)),exist=exists)
    if (exists) then
      attempted = .true.
      ios = tec_read_structured_multiblock(orion=orion, filename=trim(transfile(2)))
    endif
    if (ios/=0) then
      if (attempted) ios = 2
      return
    endif
    dummy1  = lbound(orion%block(1)%mesh, dim=2)
    dummy23 = lbound(orion%block(1)%mesh, dim=3)
    Tmin = nint(orion%block(1)%mesh(1,dummy1,dummy23,dummy23))
    Tmax = Tmin + ubound(orion%block(1)%mesh, dim=2) - dummy1
    start = Tmin
    if (Tmin==1) Tmin = 0
    allocate(mi_tab(Tmin:Tmax, 1:ns))
    allocate(k_tab(Tmin:Tmax, 1:ns))
    dummy23 = lbound(orion%block(1)%vars, dim=3)
    do i = 1, ns
      mi_tab(start:Tmax,i) = orion%block(i)%vars(1,:,dummy23,dummy23)
      k_tab(start:Tmax,i) = orion%block(i)%vars(2,:,dummy23,dummy23)
    enddo
    if (Tmin==0) then
      mi_tab(0,:) = mi_tab(1,:)
      k_tab(0,:) = k_tab(1,:)
    endif

  end function read_idealgas_transport


end module FLINT_Load_ThermoTransport
