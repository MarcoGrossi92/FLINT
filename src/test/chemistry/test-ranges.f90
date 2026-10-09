! test-ranges: the temperature-grid contract of the tables (fixtures made by test/chemistry/ranges/make_ranges.py
! from database/WD, database/CORIA and database/Gerlinger; thermo tables 1400..1600 K):
!  - rate tables on the thermo grid: accepted; starting above it (1450..1600 K): refused (ios = 6;
!    row T of a rate table is the rate at T kelvin, the source terms need both tables at every T);
!  - rate tables wider than the thermo grid (database/WD, 1..15000 K): accepted, and rows 1400..1600
!    are bit-identical to the tables on the thermo grid (rates-equal);
!  - a falloff table on a grid other than the Arrhenius one: refused (ios = 6; before this check the
!    rows were copied into arrays of another shape);
!  - a transport table that starts above the thermo grid (1420..1600 K): refused (ios = 3; the
!    species-contiguous copy would read rows below the table's first row); one that starts below it
!    (1380..1600 K): accepted, with the rows, the copies and the mixture viscosity and conductivity
!    of the table on the thermo grid; the same two cases for a binary-diffusion table;
!  - the same for the LAST row (rates-short, rates-troe-short), for a falloff-Lindemann table and for a
!    binary-diffusion table (ios = 3);
!  - in a child process (argument child-grid) the five refusals are printed on the error unit too,
!    with the reason for the two tables that start above the thermo grid, and the two tables that
!    start below it print none;
!  - a table with fewer zones than reactions of its type (ios = 4) and tables on a 2 K step (ios = 6 / 4);
!  - a negative k_inf / k_0 in a falloff table (ios = 5), and a NaN F_cent (ios = 5, also under the
!    -ffast-math of RELEASE, where the former test Fcent /= Fcent was folded to .false.);
!  - every row of every zone, not only the first and the last row of zone 1: a rate table with zones
!    2.. on another grid (rates-zone2-shift) and one with an interior row missing and another one twice
!    (rates-gap) are refused (ios = 6; before, the rows were copied under other temperatures); the same
!    gap in the last zone of a transport and of a binary-diffusion table (ios = 3) and in zone 1 of the
!    thermo table (thermo-gap, ios = 4).
! Exit code 1 on failure.
program test
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  use FLINT_Load_chemistry
  use FLINT_Lib_Chemistry_data
  implicit none
  integer :: err, nfail
  real(8), allocatable :: kf_eq(:,:), kb_eq(:,:)
  real(8), allocatable :: mi_eq(:,:), k_eq(:,:), miT_eq(:,:), kT_eq(:,:), dij_eq(:,:), rhoi(:), Dm_eq(:,:), Dm(:,:)
  real(8) :: mil_eq(4), kl_eq(4), mil(4), kl(4)
  character(32) :: mech_name
  character(len=512) :: self, arg

  call get_command_argument(0, self); call get_command_argument(1, arg)
  if (arg == 'child-grid') then
    err = read_idealgas_thermo('ranges/thermo-1400/')
    err = read_chemistry(folder='ranges/rates-narrow', mech_name=mech_name); call free_chemistry_data()
    err = read_chemistry(folder='ranges/rates-troe-mismatch', mech_name=mech_name); call free_chemistry_data()
    err = read_chemistry(folder='ranges/rates-lind-mismatch', mech_name=mech_name); call free_chemistry_data()
    err = read_idealgas_transport('ranges/transport-shifted/')
    err = read_idealgas_diffusion('ranges/diffusion-shifted/')
    err = read_idealgas_transport('ranges/transport-lower/')
    err = read_idealgas_diffusion('ranges/diffusion-lower/')
    stop
  endif
  nfail = 0
  err = read_idealgas_thermo('ranges/thermo-1400/')
  call verdict('thermo tables 1400..1600 K loaded', err == 0 .and. merge(1, Tmin, Tmin == 0) == 1400 .and. Tmax == 1600)

  err = read_chemistry(folder='ranges/rates-equal', mech_name=mech_name)
  call verdict('rate tables on the thermo grid (1400..1600 K): accepted', err == 0 .and. T_tab_min == 1400 .and. T_tab_max == 1600)
  kf_eq = kf_tab; kb_eq = kb_tab
  call free_chemistry_data()
  ! a rate table that covers the thermo grid and extends beyond it is read at T kelvin
  err = read_chemistry(folder='../../database/WD', mech_name=mech_name)
  call verdict('rate tables wider than the thermo grid (1..15000 K over 1400..1600 K): accepted', &
    err == 0 .and. T_tab_min == 1 .and. T_tab_max == 15000)
  if (err == 0) then
    call verdict('wider rate tables: rows 1400..1600 K bit-identical to the tables on the thermo grid', &
      all(kf_tab(1400:1600,:) == kf_eq) .and. all(kb_tab(1400:1600,:) == kb_eq))
  endif
  call free_chemistry_data()
  err = read_chemistry(folder='ranges/rates-narrow', mech_name=mech_name)
  call verdict('rate tables on another grid than the thermo one (1450..1600 K): refused with ios = 6', err == 6)
  call free_chemistry_data()
  err = read_chemistry(folder='ranges/rates-troe-equal', mech_name=mech_name)
  call verdict('falloff-Troe table on the Arrhenius grid: accepted', err == 0 .and. nrc_troe == 1)
  call free_chemistry_data()
  err = read_chemistry(folder='ranges/rates-troe-mismatch', mech_name=mech_name)
  call verdict('falloff-Troe table on another grid (1450..1600 K): refused with ios = 6', err == 6)
  call free_chemistry_data()

  err = read_idealgas_transport('ranges/transport-shifted/')
  call verdict('transport table starting at 1420 K with thermo from 1400 K: refused with ios = 3', err == 3)
  err = read_idealgas_transport('ranges/transport-gap/')
  call verdict('transport table with row 1500 K missing and row 1499 K twice in its last zone: refused with ios = 3', err == 3)
  err = read_idealgas_transport('ranges/transport-equal/')
  call verdict('transport table on the thermo grid: accepted', err == 0)
  ! a transport table that starts below the thermo grid covers it: accepted, row T is the value at
  ! T kelvin, and the transport routines give the values of the table on the thermo grid
  allocate(rhoi(ns)); rhoi = 1d0/ns
  if (allocated(mi_tab)) then
    mi_eq = mi_tab(1400:1600,:); k_eq = k_tab(1400:1600,:); miT_eq = mi_tabT; kT_eq = k_tabT
    call transport_values(mil_eq, kl_eq)
    deallocate(mi_tab, k_tab)
  endif
  err = read_idealgas_transport('ranges/transport-lower/')
  call verdict('transport table starting at 1380 K with thermo from 1400 K: accepted', err == 0)
  if (err == 0 .and. allocated(mi_eq)) then
    call transport_values(mil, kl)
    call verdict('transport table from 1380 K: rows 1400..1600 K, their copies and the mixture viscosity and ' // &
      'conductivity bit-identical to the table on the thermo grid', lbound(mi_tab, 1) == 1380 .and. &
      all(mi_tab(1400:1600,:) == mi_eq) .and. all(k_tab(1400:1600,:) == k_eq) .and. all(mi_tabT == miT_eq) .and. &
      all(k_tabT == kT_eq) .and. all(mil == mil_eq) .and. all(kl == kl_eq))
  endif

  ! the last row too, the falloff-Lindemann and the binary-diffusion tables
  err = read_chemistry(folder='ranges/rates-short', mech_name=mech_name)
  call verdict('rate tables ending before the thermo ones (1400..1550 K): refused with ios = 6', err == 6)
  call free_chemistry_data()
  err = read_chemistry(folder='ranges/rates-troe-short', mech_name=mech_name)
  call verdict('falloff-Troe table ending before the Arrhenius one (1400..1550 K): refused with ios = 6', err == 6)
  call free_chemistry_data()
  err = read_chemistry(folder='ranges/rates-lind-equal', mech_name=mech_name)
  call verdict('falloff-Lindemann table on the Arrhenius grid: accepted', err == 0 .and. nrc_lindemann == 1)
  call free_chemistry_data()
  err = read_chemistry(folder='ranges/rates-lind-mismatch', mech_name=mech_name)
  call verdict('falloff-Lindemann table on another grid (1450..1600 K): refused with ios = 6', err == 6)
  call free_chemistry_data()
  err = read_idealgas_diffusion('ranges/diffusion-shifted/')
  call verdict('diffusion table starting at 1420 K with thermo from 1400 K: refused with ios = 3', err == 3)
  err = read_idealgas_diffusion('ranges/diffusion-gap/')
  call verdict('diffusion table with row 1500 K missing and row 1499 K twice in its last pair: refused with ios = 3', err == 3)
  err = read_idealgas_diffusion('ranges/diffusion-equal/')
  call verdict('diffusion table on the thermo grid: accepted', err == 0)
  ! the same for a binary-diffusion table that starts below the thermo grid (other coefficients below 1400 K)
  allocate(Dm_eq(ns,4), Dm(ns,4))
  if (allocated(dij_tab)) then
    dij_eq = dij_tab(1400:1600,:)
    call diffusion_values(Dm_eq)
    deallocate(dij_tab)
  endif
  err = read_idealgas_diffusion('ranges/diffusion-lower/')
  call verdict('diffusion table starting at 1380 K with thermo from 1400 K: accepted', err == 0)
  if (err == 0 .and. allocated(dij_eq)) then
    call diffusion_values(Dm)
    call verdict('diffusion table from 1380 K: rows 1400..1600 K and the mixture-averaged diffusion ' // &
      'coefficients bit-identical to the table on the thermo grid', lbound(dij_tab, 1) == 1380 .and. &
      all(dij_tab(1400:1600,:) == dij_eq) .and. all(Dm == Dm_eq))
  endif
  ! a table with fewer zones than reactions of its type (e.g. a falloff-SRI reaction counted as Arrhenius)
  err = read_chemistry(folder='ranges/rates-missing-zone', mech_name=mech_name)
  call verdict('Arrhenius table with 2 zones for 3 reactions: refused with ios = 4', err == 4)
  call free_chemistry_data()
  ! a negative limiting rate coefficient in a falloff table
  err = read_chemistry(folder='ranges/rates-troe-negk', mech_name=mech_name)
  call verdict('falloff-Troe table with k_inf < 0 at 1500 K: refused with ios = 5', err == 5)
  call free_chemistry_data()
  err = read_chemistry(folder='ranges/rates-troe-nan', mech_name=mech_name)
  call verdict('falloff-Troe table with F_cent = NaN at 1500 K: refused with ios = 5', err == 5)
  call free_chemistry_data()
  err = read_chemistry(folder='ranges/rates-lind-negk', mech_name=mech_name)
  call verdict('falloff-Lindemann table with k_0 < 0 at 1500 K: refused with ios = 5', err == 5)
  call free_chemistry_data()
  ! rows on a 2 K step: the first and the computed last row match the thermo grid, the temperatures do not
  err = read_chemistry(folder='ranges/rates-step2', mech_name=mech_name)
  call verdict('rate table on a 2 K step (1400..1800 K, 201 rows): refused with ios = 6', err == 6)
  call free_chemistry_data()
  ! every zone and every row: zone 1 alone on the 1 K grid, or its first and last row only, is not enough
  err = read_chemistry(folder='ranges/rates-zone2-shift', mech_name=mech_name)
  call verdict('rate table with zone 1 on 1400..1600 K and zones 2.. on 1401..1601 K: refused with ios = 6', err == 6)
  call free_chemistry_data()
  err = read_chemistry(folder='ranges/rates-gap', mech_name=mech_name)
  call verdict('rate table without row 1500 K and with row 1499 K twice (first row, last row and row count ' // &
    'of the 1 K grid): refused with ios = 6', err == 6)
  call free_chemistry_data()
  call execute_command_line(trim(self)//' child-grid > ranges/child-grid.out 2> ranges/child-grid.err', exitstat=err)
  call execute_command_line('test "$(command grep -c -F ''[ERROR] FLINT read_'' ranges/child-grid.err)" = 5', exitstat=err)
  call verdict('child process: the five grid refusals are on the error unit too', err == 0)
  call execute_command_line('test "$(command grep -c -F ''must start at or below the first thermo temperature'' ' // &
    'ranges/child-grid.err)" = 2', exitstat=err)
  call verdict('child process: the transport and diffusion tables that start above the thermo grid are refused ' // &
    'with that reason', err == 0)

  ! thermo tables on a 2 K step and with an interior gap (last: the phase arrays are reloaded)
  call free_thermo()
  err = read_idealgas_thermo('ranges/thermo-step2/')
  call verdict('thermo table on a 2 K step (1400..1800 K): refused with ios = 4', err == 4)
  call free_thermo()
  err = read_idealgas_thermo('ranges/thermo-gap/')
  call verdict('thermo table without row 1500 K and with row 1499 K twice (first row, last row and row count ' // &
    'of the 1 K grid): refused with ios = 4', err == 4)

  if (nfail > 0) then
    write(*,'(A,I0,A)') ' Verdict -> fail (', nfail, ' checks)'
    stop 1
  endif
  write(*,'(A)') ' Verdict -> pass'
contains
  subroutine verdict(what, good)
    character(*), intent(in) :: what
    logical, intent(in) :: good
    if (good) then
      write(*,'(A)') ' [ok]   '//what
    else
      write(*,'(A)') ' [FAIL] '//what
      nfail = nfail + 1
    endif
  end subroutine verdict
  ! the phase arrays and thermo tables, before reloading phase.txt and thermo.dat
  subroutine free_thermo()
    if (allocated(species_names)) deallocate(species_names)
    if (allocated(wm_tab)) deallocate(wm_tab)
    if (allocated(Ri_tab)) deallocate(Ri_tab)
    if (allocated(h_tab)) deallocate(h_tab)
    if (allocated(s_tab)) deallocate(s_tab)
    if (allocated(cp_tab)) deallocate(cp_tab)
    if (allocated(dcpi_tab)) deallocate(dcpi_tab)
  end subroutine free_thermo
  ! mixture viscosity and conductivity (Wilke) at four temperatures of the thermo grid, equal densities
  subroutine transport_values(mil, kl)
    real(8), intent(out) :: mil(4), kl(4)
    real(8), parameter :: Ts(4) = [1400d0, 1400.5d0, 1455.25d0, 1599.75d0]
    integer :: j
    do j = 1, 4
      call co_k_mi_lam_Wilke(rhoi, sum(rhoi), Ts(j), mil(j), kl(j))
    enddo
  end subroutine transport_values
  ! mixture-averaged diffusion coefficients at the same four temperatures, at 1 atm
  subroutine diffusion_values(D)
    real(8), intent(out) :: D(:,:)
    integer, parameter :: Ti(4) = [1400, 1400, 1455, 1599]
    real(8), parameter :: Td(4) = [0d0, 0.5d0, 0.25d0, 0.75d0]
    integer :: j, Tint(2)
    do j = 1, 4
      Tint = [Ti(j), Ti(j)+1]
      call co_DS_expr(rhoi, sum(rhoi), Tint, Td(j), 101325d0, D(:,j))
    enddo
  end subroutine diffusion_values
end program test
