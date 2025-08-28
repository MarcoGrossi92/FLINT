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
    integer           :: ios, i, unitfile, start
    character(256)    :: wholestring, args(2)
    character(512)    :: wmfile, thermofile
    type(ORION_data)  :: orion

    if (present(folder)) then
      wmfile = trim(folder)//'/'//trim(FLINT_phase_prefix)//'phase.txt'
      thermofile = trim(folder)//'/'//trim(FLINT_phase_prefix)//'thermo.dat'
    else
      wmfile = 'INPUT/'//trim(FLINT_phase_prefix)//'phase.txt'
      thermofile = 'INPUT/'//trim(FLINT_phase_prefix)//'thermo.dat'
    endif

    open(newunit=unitFile,file=trim(wmfile),status='old',iostat=ios)
    if (ios/=0) then
      write(*,*) '[ERROR] phase.txt type file not found'
      return
    endif
    ios = 0; ns = -1
    read(unitFile,*)
    do while(ios==0)
      read(unitFile,'(A)',iostat=ios) wholestring
      ns = ns + 1
    enddo
    allocate(wm_tab(1:ns))
    allocate(Ri_tab(1:ns))
    rewind(unitFile)
    read(unitFile,*)
    do i = 1, ns
      read(unitFile,'(A)') wholestring
      call parse(wholestring,' ',args)
      read(args(2),*) wm_tab(i)
    end do
    close(unitFile)

    Ri_tab = Runiv/wm_tab

    ios = tec_read_points_multivars(orion,4,trim(thermofile))
    if (ios/=0) then
      write(*,*) '[ERROR] Reading thermo file'
      return
    endif
    Tmin = nint(orion%block(1)%mesh(1,1,1,1))
    Tmax = Tmin + orion%block(1)%Ni - 1
    start = Tmin
    if (Tmin==1) Tmin = 0
    allocate(h_tab(Tmin:Tmax, 1:ns))
    allocate(s_tab, cp_tab, dcpi_tab, mold=h_tab)
    do i = 1, ns
      cp_tab(start:Tmax,i) = orion%block(i)%vars(1,:,1,1)
      h_tab(start:Tmax,i) = orion%block(i)%vars(2,:,1,1)
      s_tab(start:Tmax,i) = orion%block(i)%vars(3,:,1,1)
      dcpi_tab(start:Tmax,i) = orion%block(i)%vars(4,:,1,1)
    enddo
    if (Tmin == 0) then
      cp_tab(0,:) = cp_tab(1,:)
      h_tab(0,:) = cp_tab(1,:)*0.d0
      s_tab(0,:) = s_tab(1,:)
      dcpi_tab(0,:) = dcpi_tab(1,:)
    endif

  end function read_idealgas_thermo


  function read_idealgas_transport(folder) result(ios)
    use strings, only: parse
    use Lib_Tecplot
    use Lib_ORION_data
    implicit none
    character(len=*), intent(in), optional :: folder
    ! Local
    integer           :: ios, i, unitfile, start
    character(256)    :: wholestring, args(2)
    character(512)    :: transfile
    type(ORION_data)  :: orion

    if (present(folder)) then
      transfile = trim(folder)//'/'//trim(FLINT_phase_prefix)//'transport.dat'
    else
      transfile = 'INPUT/'//trim(FLINT_phase_prefix)//'transport.dat'
    endif

    ios = tec_read_points_multivars(orion,2,trim(transfile))
    if (ios/=0) return
    Tmin = nint(orion%block(1)%mesh(1,1,1,1))
    Tmax = Tmin + orion%block(1)%Ni - 1
    start = Tmin
    if (Tmin==1) Tmin = 0
    allocate(mi_tab(Tmin:Tmax, 1:ns))
    allocate(k_tab(Tmin:Tmax, 1:ns))
    do i = 1, ns
      mi_tab(start:Tmax,i) = orion%block(i)%vars(1,:,1,1)
      k_tab(start:Tmax,i) = orion%block(i)%vars(2,:,1,1)
    enddo
    if (Tmin==0) then
      mi_tab(0,:) = mi_tab(1,:)
      k_tab(0,:) = k_tab(1,:)
    endif

  end function read_idealgas_transport


end module FLINT_Load_ThermoTransport
