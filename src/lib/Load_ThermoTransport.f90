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

! read_realfluid_thermo ios map:
! 0 -> no error
! 1 -> phase.txt not found
! 2 -> phase.txt found but error reading
! 3 -> thermo file not found
! 4 -> thermo file found but error reading
! 5 -> too many blocks in thermo file
!
! read_realfluid_transport ios map:
! 0 -> no error
! 1 -> transport file not found
! 2 -> error reading transport data
! 3 -> too many blocks in transport file
! 4 -> error reading transport data: mesh size mismatch with thermo data
!
! ph2pT ios map:
! 0 -> no error
! 1 -> h_tab2D already allocated
! 2 -> table dimensions too small (Nh<1 or Np<0)
! 3 -> T_tab not allocated
! 4 -> found = .false. (interpolation failed); fatal errors override warnings


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

    ! Precompute species-pair constants used by co_fiij (Wilke mixing rule)
    if (allocated(Mi_Mj_pow_m025)) deallocate(Mi_Mj_pow_m025)
    if (allocated(inv_sqrt8_1p))   deallocate(inv_sqrt8_1p)
    allocate(Mi_Mj_pow_m025(1:ns, 1:ns))
    allocate(inv_sqrt8_1p  (1:ns, 1:ns))
    block
      integer :: si, sj
      real(kind=8) :: ratio
      do sj = 1, ns
        do si = 1, ns
          ratio = wm_tab(si) / wm_tab(sj)
          Mi_Mj_pow_m025(si,sj) = ratio**(-0.25d0)
          inv_sqrt8_1p  (si,sj) = 1.d0 / dsqrt(8.d0 * (1.d0 + ratio))
        end do
      end do
    end block

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
    integer           :: Tmin_dummy, Tmax_dummy
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
    Tmin_dummy = nint(orion%block(1)%mesh(1,dummy1,dummy23,dummy23))
    Tmax_dummy = Tmin_dummy + ubound(orion%block(1)%mesh, dim=2) - dummy1

    if (Tmax_dummy /= Tmax) then
      ios = 3
      return
    endif

    start = Tmin_dummy
    if (Tmin_dummy==1) Tmin_dummy = 0
    allocate(mi_tab(Tmin_dummy:Tmax, 1:ns))
    allocate(k_tab(Tmin_dummy:Tmax, 1:ns))
    dummy23 = lbound(orion%block(1)%vars, dim=3)
    do i = 1, ns
      mi_tab(start:Tmax,i) = orion%block(i)%vars(1,:,dummy23,dummy23)
      k_tab(start:Tmax,i) = orion%block(i)%vars(2,:,dummy23,dummy23)
    enddo
    if (Tmin_dummy==0) then
      mi_tab(0,:) = mi_tab(1,:)
      k_tab(0,:) = k_tab(1,:)
    endif

  end function read_idealgas_transport


  ! Reads binary diffusion coefficients D_ij(T) from <prefix>diffusion.dat.
  ! One ZONE per unique unordered species pair, upper-triangular order
  ! (1,2),(1,3),..,(1,ns),(2,3),..  Each zone holds a single variable (Dij) vs T,
  ! on the same temperature grid as the thermo/transport tables.
  ! Populates dij_tab(Tmin:Tmax, 1:ndij) with ndij = ns*(ns-1)/2.
  function read_idealgas_diffusion(folder) result(ios)
    use Lib_Tecplot
    use Lib_ORION_data
    implicit none
    character(len=*), intent(in), optional :: folder
    ! Local
    integer           :: ios, i, start, dummy1, dummy23
    integer           :: Tmin_dummy, Tmax_dummy
    integer           :: u, ios2, ip
    logical           :: exists, attempted
    character(512)    :: difffile(2), titleline
    type(ORION_data)  :: orion

    ios = 1
    attempted = .false.
    dij_pref = 101325.d0   ! default reference pressure (1 atm) if TITLE has no Pref tag

    ! Single-species mixtures have no pairs: nothing to load.
    ndij = ns*(ns-1)/2
    if (ndij < 1) return

    if (present(folder)) then
      difffile(1) = trim(folder)//'/'//trim(FLINT_phase_prefix)//'diffusion.dat'
      difffile(2) = trim(folder)//'/'//trim(FLINT_phase_prefix)//'diffusion.szplt'
    else
      difffile(1) = 'INPUT/'//trim(FLINT_phase_prefix)//'diffusion.dat'
      difffile(2) = 'INPUT/'//trim(FLINT_phase_prefix)//'diffusion.szplt'
    endif

    inquire(file=trim(difffile(1)),exist=exists)
    if (exists) then
      attempted = .true.
      ! Parse the reference pressure from the TITLE line: "... (Pref=<value> Pa)".
      open(newunit=u, file=trim(difffile(1)), status='old', iostat=ios2)
      if (ios2==0) then
        read(u,'(A)',iostat=ios2) titleline
        if (ios2==0) then
          ip = index(titleline,'Pref=')
          if (ip>0) read(titleline(ip+5:),*,iostat=ios2) dij_pref
        endif
        close(u)
      endif
      ios = tec_read_points_multivars(orion,1,trim(difffile(1)))
    endif
    inquire(file=trim(difffile(2)),exist=exists)
    if (exists) then
      attempted = .true.
      ios = tec_read_structured_multiblock(orion=orion, filename=trim(difffile(2)))
    endif
    if (ios/=0) then
      if (attempted) ios = 2
      return
    endif

    if (size(orion%block) /= ndij) then
      ios = 4
      return
    endif

    dummy1  = lbound(orion%block(1)%mesh, dim=2)
    dummy23 = lbound(orion%block(1)%mesh, dim=3)
    Tmin_dummy = nint(orion%block(1)%mesh(1,dummy1,dummy23,dummy23))
    Tmax_dummy = Tmin_dummy + ubound(orion%block(1)%mesh, dim=2) - dummy1

    if (Tmax_dummy /= Tmax) then
      ios = 3
      return
    endif

    start = Tmin_dummy
    if (Tmin_dummy==1) Tmin_dummy = 0
    allocate(dij_tab(Tmin_dummy:Tmax, 1:ndij))
    dummy23 = lbound(orion%block(1)%vars, dim=3)
    do i = 1, ndij
      dij_tab(start:Tmax,i) = orion%block(i)%vars(1,:,dummy23,dummy23)
    enddo
    if (Tmin_dummy==0) dij_tab(0,:) = dij_tab(1,:)

  end function read_idealgas_diffusion


  function read_realfluid_thermo(folder) result(ios)
    use strings, only: parse
    use Lib_Tecplot
    use Lib_ORION_data
    implicit none
    character(len=*), intent(in), optional :: folder
    ! Local
    integer              :: ios, i, j, unitfile
    character(256)       :: wholestring, args(2), phasefile, thermofile(2)
    type(ORION_data)     :: orion
    logical              :: exists_dat, exists_szplt
    real(8), allocatable :: p_(:), h_(:)

    if (present(folder)) then
      phasefile = trim(folder)//'/'//trim(FLINT_phase_prefix)//'phase.txt'
      thermofile(1) = trim(folder)//'/'//trim(FLINT_phase_prefix)//'thermo.dat'
      thermofile(2) = trim(folder)//'/'//trim(FLINT_phase_prefix)//'thermo.szplt'
    else
      phasefile = 'INPUT/'//trim(FLINT_phase_prefix)//'phase.txt'
      thermofile(1) = 'INPUT/'//trim(FLINT_phase_prefix)//'thermo.dat'
      thermofile(2) = 'INPUT/'//trim(FLINT_phase_prefix)//'thermo.szplt'
    endif

    !! File 1: phase — read name of the phase
    print*, 'Reading real fluid thermo from phase file: ', trim(phasefile)
    open(newunit=unitFile,file=trim(phasefile),status='old',iostat=ios)
    if (ios/=0) then
      ios = 1
      return
    endif
    ios = 0
    read(unitfile,*)!skip first line
    read(unitfile,'(A)',iostat=ios) wholestring
    if (ios/=0) then; close(unitfile); ios = 2; return; endif
    call parse(wholestring,' ',args)
    real_species_names(1) = trim(adjustl(args(1)))
    close(unitfile)

    !! File 2: thermo — use ORION to get mesh (p, h node coords)
    inquire(file=trim(thermofile(1)), exist=exists_dat)
    inquire(file=trim(thermofile(2)), exist=exists_szplt)
    if (.not. exists_dat .and. .not. exists_szplt) then
      ios = 3
      return
    endif
    ios = -1
    if (exists_dat)   ios = tec_read_structured_multiblock(orion=orion, filename=trim(thermofile(1)))
    if (ios/=0 .and. exists_szplt) ios = tec_read_structured_multiblock(orion=orion, filename=trim(thermofile(2)))
    if (ios/=0) then
      ios = 4
      return
    endif

    if (size(orion%block) /= 1) then
      ios = 5
      return
    endif

    ! p: I-axis → mesh(2, i=0..Ni, j=0, 0)
    ! h: J-axis → mesh(1, i=0, j=0..Nj, 0)
    allocate(p_(0:orion%block(1)%Ni))
    allocate(h_(0:orion%block(1)%Nj))
    p_ = orion%block(1)%mesh(1, :, 0, 0)
    h_ = orion%block(1)%mesh(2, 0, :, 0)
    pmin = p_(0)
    hmin = h_(0)
    deltap = p_(1) - p_(0)
    deltah = h_(1) - h_(0)
    allocate(rho_tab(0:orion%block(1)%Ni, 0:orion%block(1)%Nj))
    allocate(T_tab, dT_tab, rh_tab, hT_tab, s_tab2D, rp_tab, sound_tab, mold=rho_tab)
    ! rho_tab   = density
    ! T_tab     = temperature
    ! dT_tab    = derivative of density with respect to temperature @ constant pressure
    ! ht_tab    = derivative of enthalpy with respect to temperature @ constant pressure = cp
    ! s_tab2D   = entropy
    ! rp_tab    = derivative of density with respect to pressure @ constant enthalpy
    ! sound_tab = speed of sound
    do j = 0, orion%block(1)%Nj
      do i = 0, orion%block(1)%Ni
        rho_tab(i, j)   = orion%block(1)%vars(1, i, j, 0)
        T_tab(i, j)     = orion%block(1)%vars(2, i, j, 0)
        dT_tab(i, j)    = orion%block(1)%vars(3, i, j, 0)
        rh_tab(i, j)    = orion%block(1)%vars(4, i, j, 0)
        hT_tab(i, j)    = orion%block(1)%vars(5, i, j, 0)
        s_tab2D(i, j)   = orion%block(1)%vars(6, i, j, 0)
        rp_tab(i, j)    = orion%block(1)%vars(7, i, j, 0)
        sound_tab(i, j) = orion%block(1)%vars(8, i, j, 0)
      enddo
    enddo

  end function read_realfluid_thermo


  function read_realfluid_transport(folder) result(ios)
    use Lib_Tecplot
    use Lib_ORION_data
    implicit none
    character(len=*), intent(in), optional :: folder
    ! Local
    integer           :: ios, i, j
    type(ORION_data)  :: orion
    logical           :: exists_dat, exists_szplt
    character(512)    :: transfile(2)

    !! transport — use ORION to get mesh (p, h node coords)
    if (present(folder)) then
      transfile(1) = trim(folder)//'/'//trim(FLINT_phase_prefix)//'transport.dat'
      transfile(2) = trim(folder)//'/'//trim(FLINT_phase_prefix)//'transport.szplt'
    else
      transfile(1) = 'INPUT/'//trim(FLINT_phase_prefix)//'transport.dat'
      transfile(2) = 'INPUT/'//trim(FLINT_phase_prefix)//'transport.szplt'
    endif

    inquire(file=trim(transfile(1)), exist=exists_dat)
    inquire(file=trim(transfile(2)), exist=exists_szplt)
    if (.not. exists_dat .and. .not. exists_szplt) then
      ios = 1
      return
    endif
    ios = -1
    if (exists_dat)                  ios = tec_read_structured_multiblock(orion=orion, filename=trim(transfile(1)))
    if (ios/=0 .and. exists_szplt)   ios = tec_read_structured_multiblock(orion=orion, filename=trim(transfile(2)))
    if (ios/=0) then
      ios = 2
      return
    endif

    if (size(orion%block) /= 1) then
      ios = 3
      return
    endif

    if (orion%block(1)%Ni+1 /= size(rho_tab,1) .or. orion%block(1)%Nj+1 /= size(rho_tab,2)) then
      ios = 4
      return
    endif

    ! p: I-axis → mesh(2, i=0..Ni, j=0, 0)
    ! h: J-axis → mesh(1, i=0, j=0..Nj, 0)
    allocate(mi_tab2D(0:orion%block(1)%Ni, 0:orion%block(1)%Nj))
    allocate(k_tab2D(0:orion%block(1)%Ni, 0:orion%block(1)%Nj))
    ! mi_tab2D = viscosità dinamica
    ! k_tab2D  = conducibilità termica
    do j = 0, orion%block(1)%Nj
      do i = 0, orion%block(1)%Ni
        mi_tab2D(i, j)  = orion%block(1)%vars(1, i, j, 0)
        k_tab2D(i, j)   = orion%block(1)%vars(2, i, j, 0)
      enddo
    enddo

  end function read_realfluid_transport


  function ph2pT( ) result(ios)

    use Lib_Tecplot
    use Lib_ORION_data

    implicit none

    integer :: ios
    integer :: i, j, k
    integer :: Np, Nh
    real(8) :: T_target, frac, dT_row
    logical :: found

    ios = -1 ! preallocation

    ! Pre-flight checks
    if (.not. allocated(T_tab)) then
      ios = 3
      return
    end if
    if (allocated(h_tab2D) .or. allocated(Tmin2) .or. allocated(deltaT)) then
      ios = 1
      return
    end if

    Np = ubound(T_tab, 1)
    Nh = ubound(T_tab, 2)

    if (Nh < 1 .or. Np < 0) then
      ios = 2
      return
    end if

    allocate(h_tab2D, mold=T_tab)
    allocate(Tmin2(0:Np))
    allocate(deltaT(0:Np))

    ! Building inverse table
    do i = 0, Np
      Tmin2(i)   = T_tab(i, 0)
      deltaT(i) = (T_tab(i, Nh) - T_tab(i, 0)) / dble(Nh)

      do k = 0, Nh
        T_target = Tmin2(i) + k * deltaT(i)
        found = .false.

        do j = 0, Nh - 1
          if (T_target >= T_tab(i, j) .and. T_target <= T_tab(i, j+1)) then
            dT_row = T_tab(i, j+1) - T_tab(i, j)
            frac   = (T_target - T_tab(i, j)) / dT_row
            h_tab2D(i, k) = (hmin + j * deltah) + frac * deltah
            found = .true.
            exit
          end if
        end do

        ! Clamp @ boundary (k=0 e k=Nh are limit cases)
        if (.not. found) then
          if (k == 0) then
            h_tab2D(i, k) = hmin
            found = .true.
          else if (k == Nh) then
            h_tab2D(i, k) = hmin + Nh * deltah
            found = .true.
          end if
        end if

        if (.not. found) then
          deallocate(h_tab2D, Tmin2, deltaT)
          ios = 4
          return
        end if
      end do
    end do

    ios = 0
  endfunction ph2pT

end module FLINT_Load_ThermoTransport
