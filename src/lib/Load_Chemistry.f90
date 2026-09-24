! ios map
! 0 -> no error
! 1 -> info file not found
! 2 -> error reading info file
! 3 -> table file not found
! 4 -> error reading table file (also: fewer zones than reactions of that type)
! 5 -> table value not admissible (a negative falloff rate coefficient, or a non-finite F_cent)
! 6 -> chemistry-Arrhenius.dat does not cover the thermo temperature range, a falloff table is not
!      on the grid of chemistry-Arrhenius.dat, or a rate table is not on a 1 K step

module FLINT_Load_Chemistry
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none

contains

  function read_chemistry( folder, mech_name ) result(ios)
    use Lib_ORION_Data
    use Lib_Tecplot
    use FLINT_Lib_Chemistry_data
    use FLINT_Lib_Thermodynamic, only: ns, FLINT_phase_prefix, species_names, Tmin, Tmax, cp_tab
    use iso_fortran_env, only: error_unit
    implicit none
    character(len=*), intent(in), optional :: folder
    character(len=*), intent(out), optional :: mech_name
    integer :: idum, unitfile
    type(ORION_data)  :: orion
    integer :: i, j, ios, j0, j1, j2
    integer :: Ti1, Ti2, dummy1, dummy23, Tt1, Tt2, Tf
    character(len=32):: chardum
    character(len=256) :: line
    integer :: nord, isp
    real(8) :: order

    nrc_arrh = 0
    nrc_troe = 0
    nrc_lindemann = 0

    !! Info
    if (present(folder)) then
      open(newunit=unitfile,file=trim(folder)//'/'//trim(FLINT_phase_prefix)//'chemistry-info.txt',form='formatted',status='old',action='read',iostat=ios)
    else
      open(newunit=unitfile,file='INPUT/'//trim(FLINT_phase_prefix)//'chemistry-info.txt',form='formatted',status='old',action='read',iostat=ios)
    endif
    if (ios/=0) then
      ios = 1
      return
    endif
    ! Read mechanism name: the whole first line, leading/trailing blanks removed.
    ! (A list-directed read cut the name at the first blank, comma or slash:
    ! 'Aramco 2.0' became 'Aramco'. Every name hooked in Assign_Mechanism is a
    ! single token, so existing INPUT folders select the same routine as before.)
    if (present(mech_name)) then
      read(unitfile,'(A)',iostat=ios) mech_name
      if (ios/=0) then; ios = 2; close(unitfile); return; endif
      ! TAB and CR were blanks for the list-directed read this replaces (hand-edited or
      ! CRLF files): turn them into blanks before trimming, so 'WD<TAB>' still hooks WD.
      do i = 1, len_trim(mech_name)
        if (mech_name(i:i) == achar(9) .or. mech_name(i:i) == achar(13)) mech_name(i:i) = ' '
      enddo
      mech_name = trim(adjustl(mech_name))
    else
      read(unitfile,*,iostat=ios)
      if (ios/=0) then; ios = 2; close(unitfile); return; endif
    endif
    read(unitfile,*,iostat=ios)
    read(unitfile,'(A17,I4)',iostat=ios) chardum, nrc
    if (ios/=0) then; ios = 2; close(unitfile); return; endif
    do j = 1, 4; read(unitfile,*,iostat=ios); enddo
    ! Read reaction type
    allocate(rxn_type(1:nrc))
    do j = 1, nrc
      read(unitfile,*,iostat=ios) idum, chardum
      if (ios/=0) then; ios = 2; close(unitfile); return; endif
      if (index(trim(chardum),'Troe')>0) then
        nrc_troe = nrc_troe + 1
        rxn_type(j) = 1
      elseif (index(trim(chardum),'Lindemann')>0) then
        nrc_lindemann = nrc_lindemann + 1
        rxn_type(j) = 2
      else
        nrc_arrh = nrc_arrh + 1
        rxn_type(j) = 0
      endif
    enddo
    ! Read info for general loop
    allocate(ni1_arrh_tab(1:ns+1,1:nrc_arrh))
    allocate(ni2_arrh_tab(1:ns+1,1:nrc_arrh))
    allocate(epsch_arrh_tab(1:ns+1,1:nrc_arrh))
    if (nrc_troe>0) then
      allocate(ni1_troe_tab(1:ns+1,1:nrc_troe))
      allocate(ni2_troe_tab(1:ns+1,1:nrc_troe))
      allocate(epsch_troe_tab(1:ns+1,1:nrc_troe))
    endif
    if (nrc_lindemann>0) then
      allocate(ni1_lind_tab(1:ns+1,1:nrc_lindemann))
      allocate(ni2_lind_tab(1:ns+1,1:nrc_lindemann))
      allocate(epsch_lind_tab(1:ns+1,1:nrc_lindemann))
    endif
    ! Read reaction info
    read(unitfile,*,iostat=ios)
    read(unitfile,*,iostat=ios)
    j0 = 0; j1 = 0; j2 = 0
    do j = 1, nrc
      if (rxn_type(j)==0) then
        j0 = j0+1
        do i = 1, ns+1
          read(unitfile,*,iostat=ios)idum,chardum,ni1_arrh_tab(i,j0),ni2_arrh_tab(i,j0),epsch_arrh_tab(i,j0)
        enddo
      elseif (rxn_type(j)==1) then
        j1 = j1+1
        do i = 1, ns+1
          read(unitfile,*,iostat=ios)idum,chardum,ni1_troe_tab(i,j1),ni2_troe_tab(i,j1),epsch_troe_tab(i,j1)
        enddo
      elseif (rxn_type(j)==2) then
        j2 = j2+1
        do i = 1, ns+1
          read(unitfile,*,iostat=ios)idum,chardum,ni1_lind_tab(i,j2),ni2_lind_tab(i,j2),epsch_lind_tab(i,j2)
        enddo
      endif
    enddo 
    ! Trailing block: explicit forward reaction orders; a table writer always ends the file with
    ! it (n = 0 when the mechanism has none). An older INPUT folder has no block: the general
    ! procedure then warns once (warn_no_orders_block) and uses the reactant coefficients.
    !   (blank lines)
    !   Reaction orders
    !   <n>
    !   <ir> <species name> <order>     n rows; ir = index in the 'Reaction type' list above
    ! Species not listed keep their stoichiometric reactant coefficient as order (real, as in
    ! Cantera's mass-action law). Orders on falloff reactions are not supported (ios = 2).
    have_orders = .false.
    allocate(ord_arrh_tab(1:ns, 1:nrc_arrh))
    ord_arrh_tab = ni1_arrh_tab(1:ns, 1:nrc_arrh)
    do
      read(unitfile,'(A)',iostat=ios) line
      if (ios /= 0) exit                                   ! end of file: no block
      if (len_trim(line) == 0) cycle
      if (trim(adjustl(line)) /= 'Reaction orders') exit  ! anything else: not a block
      read(unitfile,*,iostat=ios) nord
      if (ios/=0) then; ios = 2; close(unitfile); return; endif
      do i = 1, nord
        read(unitfile,*,iostat=ios) idum, chardum, order
        if (ios/=0) then; ios = 2; close(unitfile); return; endif
        if (idum < 1 .or. idum > nrc) then
          write(*,'(A,I0,A)') '[ERROR] FLINT read_chemistry: Reaction orders: reaction ', idum, ' does not exist'
          write(error_unit,'(A,I0,A)') '[ERROR] FLINT read_chemistry: Reaction orders: reaction ', idum, ' does not exist'
          ios = 2; close(unitfile); return
        endif
        if (rxn_type(idum) /= 0) then
          write(*,'(A,I0,A)') '[ERROR] FLINT read_chemistry: Reaction orders: reaction ', idum, &
            ' is a falloff reaction (orders are supported for Arrhenius-type reactions only)'
          write(error_unit,'(A,I0,A)') '[ERROR] FLINT read_chemistry: Reaction orders: reaction ', idum, &
            ' is a falloff reaction (orders are supported for Arrhenius-type reactions only)'
          ios = 2; close(unitfile); return
        endif
        isp = 0
        do j = 1, ns
          if (trim(species_names(j)) == trim(chardum)) then; isp = j; exit; endif
        enddo
        if (isp == 0) then
          write(*,'(A)') '[ERROR] FLINT read_chemistry: Reaction orders: species '//trim(chardum)//' is not in phase.txt'
          write(error_unit,'(A)') '[ERROR] FLINT read_chemistry: Reaction orders: species '//trim(chardum)//' is not in phase.txt'
          ios = 2; close(unitfile); return
        endif
        ord_arrh_tab(isp, count(rxn_type(1:idum) == 0)) = order
      enddo
      have_orders = .true.
      exit
    enddo
    ios = 0
    close(unitfile)

    !! Rate Arrhenius
    if (present(folder)) then
      ios = tec_read_points_multivars(orion,2,trim(folder)//'/'//trim(FLINT_phase_prefix)//'chemistry-Arrhenius.dat')
      if (ios/=0) ios = tec_read_structured_multiblock(orion=orion, filename=trim(folder)//'/'//trim(FLINT_phase_prefix)//'chemistry-Arrhenius.szplt')
    else
      ios = tec_read_points_multivars(orion,2,trim('INPUT/')//trim(FLINT_phase_prefix)//'chemistry-Arrhenius.dat')
      if (ios/=0) ios = tec_read_structured_multiblock(orion=orion, filename=trim('INPUT/')//trim(FLINT_phase_prefix)//'chemistry-Arrhenius.szplt')
    endif
    if (ios/=0) then
      ios = 3
      return
    endif
    ! A reaction type without tables of its own (e.g. falloff-SRI) is counted as Arrhenius above:
    ! the file then has fewer zones than Arrhenius-type reactions and the copy below indexed past
    ! its last zone (segmentation fault in RELEASE, bounds error under -check all).
    if (size(orion%block) < nrc_arrh) then
      write(line,'(A,I0,A,I0,A)') '[ERROR] FLINT read_chemistry: chemistry-Arrhenius.dat has ', size(orion%block), &
        ' zones for ', nrc_arrh, ' Arrhenius-type reactions (a reaction type without tables, e.g. falloff-SRI?)'
      write(*,'(A)') trim(line)
      write(error_unit,'(A)') trim(line)
      ios = 4
      return
    endif
    dummy1  = lbound(orion%block(1)%mesh, dim=2)
    dummy23 = lbound(orion%block(1)%mesh, dim=3)
    Ti1 = nint(orion%block(1)%mesh(1,dummy1,dummy23,dummy23))
    Ti2 = Ti1 + ubound(orion%block(1)%mesh, dim=2) - dummy1
    ! 1 K step: row T is the rate at T kelvin only if the last row is the first + rows - 1
    Tf = nint(orion%block(1)%mesh(1,ubound(orion%block(1)%mesh, dim=2),dummy23,dummy23))
    if (Tf /= Ti2) then
      write(line,'(A,I0,A,I0,A,I0,A)') '[ERROR] FLINT read_chemistry: chemistry-Arrhenius.dat is not on a 1 K step (', &
        Ti2 - Ti1 + 1, ' rows from ', Ti1, ' to ', Tf, ' K): row T must be the rate at T kelvin'
      write(*,'(A)') trim(line)
      write(error_unit,'(A)') trim(line)
      ios = 6
      return
    endif
    T_tab_min = Ti1; T_tab_max = Ti2   ! row T of every rate table = rate at T kelvin
    ! Table range contract: row T of a rate table is the rate at T kelvin (the
    ! tables are allocated on their own first and last row), and the source terms
    ! need a rate at every temperature of the thermo tables. A rate table that
    ! covers the thermo range is accepted, also when it extends beyond it (e.g.
    ! the 1..15000 K tables of a database folder with thermo tables regenerated
    ! on a narrower range: the rows are read at T kelvin and rhs_native/jac_native
    ! guard both ranges); a rate table that starts above or ends below the thermo
    ! tables is refused: the kinetics would read rows outside it.
    if (allocated(cp_tab)) then
      Tf = merge(1, Tmin, Tmin == 0)
      if (Ti1 > Tf .or. Ti2 < Tmax) then
        write(*,'(A,I0,A,I0,A,I0,A,I0,A)') '[ERROR] FLINT read_chemistry: chemistry-Arrhenius.dat covers ', Ti1, '..', Ti2, &
          ' K, the thermo tables ', Tf, '..', Tmax, ' K: the rate tables must cover the thermo temperature range'
        write(error_unit,'(A,I0,A,I0,A,I0,A,I0,A)') '[ERROR] FLINT read_chemistry: chemistry-Arrhenius.dat covers ', Ti1, '..', Ti2, &
          ' K, the thermo tables ', Tf, '..', Tmax, ' K: the rate tables must cover the thermo temperature range'
        ios = 6
        return
      endif
    endif
    allocate(kf_tab(Ti1:Ti2, 1:nrc_arrh))
    allocate(kb_tab(Ti1:Ti2, 1:nrc_arrh))
    dummy23 = lbound(orion%block(1)%vars, dim=3)
    do i = 1, nrc_arrh
      kf_tab(Ti1:Ti2,i) = orion%block(i)%vars(1,:,dummy23,dummy23)
      kb_tab(Ti1:Ti2,i) = orion%block(i)%vars(2,:,dummy23,dummy23)
    enddo

    !! Rate Troe
    if (nrc_troe/=0) then
      deallocate(orion%block)
      if (present(folder)) then
        ios = tec_read_points_multivars(orion,4,trim(folder)//'/'//trim(FLINT_phase_prefix)//'chemistry-Troe.dat')
        if (ios/=0) ios = tec_read_structured_multiblock(orion=orion, filename=trim(folder)//'/'//trim(FLINT_phase_prefix)//'chemistry-Troe.szplt')
      else
        ios = tec_read_points_multivars(orion,4,trim('INPUT/')//trim(FLINT_phase_prefix)//'chemistry-Troe.dat')
        if (ios/=0) ios = tec_read_structured_multiblock(orion=orion, filename=trim('INPUT/')//trim(FLINT_phase_prefix)//'chemistry-Troe.szplt')
      endif
      if (ios/=0) then
        ios = 3
        return
      endif
      if (size(orion%block) < nrc_troe) then
        write(line,'(A,I0,A,I0,A)') '[ERROR] FLINT read_chemistry: chemistry-Troe.dat has ', size(orion%block), &
          ' zones for ', nrc_troe, ' falloff-Troe reactions'
        write(*,'(A)') trim(line)
        write(error_unit,'(A)') trim(line)
        ios = 4
        return
      endif
      ! Table range contract: the falloff tables share the grid of the Arrhenius table
      dummy1  = lbound(orion%block(1)%mesh, dim=2)
      dummy23 = lbound(orion%block(1)%mesh, dim=3)
      Tt1 = nint(orion%block(1)%mesh(1,dummy1,dummy23,dummy23))
      Tt2 = Tt1 + ubound(orion%block(1)%mesh, dim=2) - dummy1
      Tf = nint(orion%block(1)%mesh(1,ubound(orion%block(1)%mesh, dim=2),dummy23,dummy23))
      if (Tf /= Tt2) then
        write(line,'(A,I0,A,I0,A,I0,A)') '[ERROR] FLINT read_chemistry: chemistry-Troe.dat is not on a 1 K step (', &
          Tt2 - Tt1 + 1, ' rows from ', Tt1, ' to ', Tf, ' K): row T must be the rate at T kelvin'
        write(*,'(A)') trim(line)
        write(error_unit,'(A)') trim(line)
        ios = 6
        return
      endif
      if (Tt1 /= Ti1 .or. Tt2 /= Ti2) then
        write(*,'(A,I0,A,I0,A,I0,A,I0,A)') '[ERROR] FLINT read_chemistry: chemistry-Troe.dat covers ', Tt1, '..', Tt2, &
          ' K, chemistry-Arrhenius.dat ', Ti1, '..', Ti2, ' K: every rate table must share one temperature grid'
        write(error_unit,'(A,I0,A,I0,A,I0,A,I0,A)') '[ERROR] FLINT read_chemistry: chemistry-Troe.dat covers ', Tt1, '..', Tt2, &
          ' K, chemistry-Arrhenius.dat ', Ti1, '..', Ti2, ' K: every rate table must share one temperature grid'
        ios = 6
        return
      endif
      dummy23 = lbound(orion%block(1)%vars, dim=3)
      allocate(Fcent_tab(Ti1:Ti2, 1:nrc_troe))
      allocate(k0_troe_tab, kinf_troe_tab, kc_troe_tab, mold=Fcent_tab)
      do i = 1, nrc_troe
        kinf_troe_tab(Ti1:Ti2,i) = orion%block(i)%vars(1,:,dummy23,dummy23)
        k0_troe_tab(Ti1:Ti2,i)   = orion%block(i)%vars(2,:,dummy23,dummy23)
        kc_troe_tab(Ti1:Ti2,i)   = orion%block(i)%vars(3,:,dummy23,dummy23)
        Fcent_tab(Ti1:Ti2,i) = orion%block(i)%vars(4,:,dummy23,dummy23)
      enddo
      ! Negative limiting rate coefficients and a non-finite F_cent are not admissible;
      ! F_cent <= 0 is (published parameter sets reach it at high T: see f_F)
      do i = 1, nrc_troe
        do j = Ti1, Ti2
          if (Fcent_tab(j,i) /= Fcent_tab(j,i) .or. abs(Fcent_tab(j,i)) > huge(1d0) .or. &
              kinf_troe_tab(j,i) < 0d0 .or. k0_troe_tab(j,i) < 0d0) then
            write(*,'(A,I0,A,I0,A)') '[ERROR] FLINT read_chemistry: chemistry-Troe.dat falloff-Troe reaction ', i, &
              ' at T = ', j, ' K: k_inf/k_0 < 0 or F_cent not finite: table not admissible'
            write(error_unit,'(A,I0,A,I0,A)') '[ERROR] FLINT read_chemistry: chemistry-Troe.dat falloff-Troe reaction ', i, &
              ' at T = ', j, ' K: k_inf/k_0 < 0 or F_cent not finite: table not admissible'
            ios = 5
            return
          endif
        enddo
      enddo
    endif
  
    !! Rate Lindemann
    if (nrc_lindemann/=0) then
      deallocate(orion%block)
      if (present(folder)) then
        ios = tec_read_points_multivars(orion,3,trim(folder)//'/'//trim(FLINT_phase_prefix)//'chemistry-Lindemann.dat')
        if (ios/=0) ios = tec_read_structured_multiblock(orion=orion, filename=trim(folder)//'/'//trim(FLINT_phase_prefix)//'chemistry-Lindemann.szplt')
      else
        ios = tec_read_points_multivars(orion,3,trim('INPUT/')//trim(FLINT_phase_prefix)//'chemistry-Lindemann.dat')
        if (ios/=0) ios = tec_read_structured_multiblock(orion=orion, filename=trim('INPUT/')//trim(FLINT_phase_prefix)//'chemistry-Lindemann.szplt')
      endif
      if (ios/=0) then
        ios = 3
        return
      endif
      if (size(orion%block) < nrc_lindemann) then
        write(line,'(A,I0,A,I0,A)') '[ERROR] FLINT read_chemistry: chemistry-Lindemann.dat has ', size(orion%block), &
          ' zones for ', nrc_lindemann, ' falloff-Lindemann reactions'
        write(*,'(A)') trim(line)
        write(error_unit,'(A)') trim(line)
        ios = 4
        return
      endif
      ! Table range contract: the falloff tables share the grid of the Arrhenius table
      dummy1  = lbound(orion%block(1)%mesh, dim=2)
      dummy23 = lbound(orion%block(1)%mesh, dim=3)
      Tt1 = nint(orion%block(1)%mesh(1,dummy1,dummy23,dummy23))
      Tt2 = Tt1 + ubound(orion%block(1)%mesh, dim=2) - dummy1
      Tf = nint(orion%block(1)%mesh(1,ubound(orion%block(1)%mesh, dim=2),dummy23,dummy23))
      if (Tf /= Tt2) then
        write(line,'(A,I0,A,I0,A,I0,A)') '[ERROR] FLINT read_chemistry: chemistry-Lindemann.dat is not on a 1 K step (', &
          Tt2 - Tt1 + 1, ' rows from ', Tt1, ' to ', Tf, ' K): row T must be the rate at T kelvin'
        write(*,'(A)') trim(line)
        write(error_unit,'(A)') trim(line)
        ios = 6
        return
      endif
      if (Tt1 /= Ti1 .or. Tt2 /= Ti2) then
        write(*,'(A,I0,A,I0,A,I0,A,I0,A)') '[ERROR] FLINT read_chemistry: chemistry-Lindemann.dat covers ', Tt1, '..', Tt2, &
          ' K, chemistry-Arrhenius.dat ', Ti1, '..', Ti2, ' K: every rate table must share one temperature grid'
        write(error_unit,'(A,I0,A,I0,A,I0,A,I0,A)') '[ERROR] FLINT read_chemistry: chemistry-Lindemann.dat covers ', Tt1, '..', Tt2, &
          ' K, chemistry-Arrhenius.dat ', Ti1, '..', Ti2, ' K: every rate table must share one temperature grid'
        ios = 6
        return
      endif
      dummy23 = lbound(orion%block(1)%vars, dim=3)
      allocate(kinf_lind_tab(Ti1:Ti2, 1:nrc_lindemann))
      allocate(k0_lind_tab, kc_lind_tab, mold=kinf_lind_tab)
      do i = 1, nrc_lindemann
        kinf_lind_tab(Ti1:Ti2,i) = orion%block(i)%vars(1,:,dummy23,dummy23)
        k0_lind_tab(Ti1:Ti2,i)   = orion%block(i)%vars(2,:,dummy23,dummy23)
        kc_lind_tab(Ti1:Ti2,i)   = orion%block(i)%vars(3,:,dummy23,dummy23)
      enddo
      do i = 1, nrc_lindemann
        do j = Ti1, Ti2
          if (kinf_lind_tab(j,i) < 0d0 .or. k0_lind_tab(j,i) < 0d0) then
            write(*,'(A,I0,A,I0,A)') '[ERROR] FLINT read_chemistry: chemistry-Lindemann.dat falloff-Lindemann reaction ', i, &
              ' at T = ', j, ' K: k_inf/k_0 < 0: table not admissible'
            write(error_unit,'(A,I0,A,I0,A)') '[ERROR] FLINT read_chemistry: chemistry-Lindemann.dat falloff-Lindemann reaction ', i, &
              ' at T = ', j, ' K: k_inf/k_0 < 0: table not admissible'
            ios = 5
            return
          endif
        enddo
      enddo
    endif

  end function read_chemistry


end module FLINT_Load_Chemistry