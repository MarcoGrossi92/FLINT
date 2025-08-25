module FLINT_IO_Table
  use FLINT_Lib_Thermodynamic
  implicit none

contains

  subroutine Read_IdealGas_Properties(folder)
    use strings, only: parse
    use Lib_Tecplot
    use Lib_ORION_data
    implicit none
    character(len=*), intent(in), optional :: folder
    ! Local
    integer           :: ios_trans, ios, i, unitfile, start
    character(256)    :: wholestring, args(2)
    character(512)    :: wmfile, thermofile, transfile
    type(ORION_data)  :: orion

    if (present(folder)) then
      wmfile = trim(folder)//'/'//trim(FLINT_phase_prefix)//'phase.txt'
      thermofile = trim(folder)//'/'//trim(FLINT_phase_prefix)//'thermo.dat'
      transfile = trim(folder)//'/'//trim(FLINT_phase_prefix)//'transport.dat'
    else
      wmfile = 'INPUT/'//trim(FLINT_phase_prefix)//'phase.txt'
      thermofile = 'INPUT/'//trim(FLINT_phase_prefix)//'thermo.dat'
      transfile = 'INPUT/'//trim(FLINT_phase_prefix)//'transport.dat'
    endif

    open(newunit=unitFile,file=trim(wmfile),status='old',iostat=ios)
    if (ios/=0) then
      write(*,*) '[ERROR] phase.txt type file not found'
      stop
    endif
    ios = 0; nsc = -1
    read(unitFile,*)
    do while(ios==0)
      read(unitFile,'(A)',iostat=ios) wholestring
      nsc = nsc + 1
    enddo
    allocate(wm_tab(1:nsc))
    allocate(Ri_tab(1:nsc))
    rewind(unitFile)
    read(unitFile,*)
    do i = 1, nsc
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
    allocate(h_tab(Tmin:Tmax, 1:nsc))
    allocate(s_tab, cp_tab, dcpi_tab, mold=h_tab)
    do i = 1, nsc
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

    deallocate(orion%block)
    ios_trans = tec_read_points_multivars(orion,2,trim(transfile))
    if (ios_trans/=0) return
    Tmin = nint(orion%block(1)%mesh(1,1,1,1))
    Tmax = Tmin + orion%block(1)%Ni - 1
    start = Tmin
    if (Tmin==1) Tmin = 0
    allocate(mi_tab(Tmin:Tmax, 1:nsc))
    allocate(k_tab(Tmin:Tmax, 1:nsc))
    do i = 1, nsc
      mi_tab(start:Tmax,i) = orion%block(i)%vars(1,:,1,1)
      k_tab(start:Tmax,i) = orion%block(i)%vars(2,:,1,1)
    enddo
    if (Tmin==0) then
      mi_tab(0,:) = mi_tab(1,:)
      k_tab(0,:) = k_tab(1,:)
    endif

  end subroutine Read_IdealGas_Properties


end module FLINT_IO_Table
