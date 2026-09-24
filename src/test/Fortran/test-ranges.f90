! test-ranges: the temperature-grid contract of the tables (fixtures made by test/ranges/make_ranges.py
! from database/WD, database/CORIA and database/Gerlinger; thermo tables 1400..1600 K):
!  - rate tables on the thermo grid: accepted; starting above it (1450..1600 K): refused (ios = 6;
!    row T of a rate table is the rate at T kelvin, the source terms need both tables at every T);
!  - rate tables wider than the thermo grid (database/WD, 1..15000 K): accepted, and rows 1400..1600
!    are bit-identical to the tables on the thermo grid (rates-equal);
!  - a falloff table on a grid other than the Arrhenius one: refused (ios = 6; before this check the
!    rows were copied into arrays of another shape);
!  - a transport table whose first row differs from the thermo one: refused (ios = 3; before this
!    check the species-contiguous copy read rows below the table's first row).
!  - the same for the LAST row (rates-short, rates-troe-short), for a falloff-Lindemann table and for a
!    binary-diffusion table (ios = 3);
!  - in a child process (argument child-grid) the five refusals are printed on the error unit too;
!  - a table with fewer zones than reactions of its type (ios = 4) and tables on a 2 K step (ios = 6 / 4);
!  - a negative k_inf / k_0 in a falloff table (ios = 5).
! Exit code 1 on failure.
program test
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  use FLINT_Load_chemistry
  use FLINT_Lib_Chemistry_data
  implicit none
  integer :: err, nfail
  real(8), allocatable :: kf_eq(:,:), kb_eq(:,:)
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
  err = read_chemistry(folder='../database/WD', mech_name=mech_name)
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
  err = read_idealgas_transport('ranges/transport-equal/')
  call verdict('transport table on the thermo grid: accepted', err == 0)

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
  err = read_idealgas_diffusion('ranges/diffusion-equal/')
  call verdict('diffusion table on the thermo grid: accepted', err == 0)
  ! a table with fewer zones than reactions of its type (e.g. a falloff-SRI reaction counted as Arrhenius)
  err = read_chemistry(folder='ranges/rates-missing-zone', mech_name=mech_name)
  call verdict('Arrhenius table with 2 zones for 3 reactions: refused with ios = 4', err == 4)
  call free_chemistry_data()
  ! a negative limiting rate coefficient in a falloff table
  err = read_chemistry(folder='ranges/rates-troe-negk', mech_name=mech_name)
  call verdict('falloff-Troe table with k_inf < 0 at 1500 K: refused with ios = 5', err == 5)
  call free_chemistry_data()
  err = read_chemistry(folder='ranges/rates-lind-negk', mech_name=mech_name)
  call verdict('falloff-Lindemann table with k_0 < 0 at 1500 K: refused with ios = 5', err == 5)
  call free_chemistry_data()
  ! rows on a 2 K step: the first and the computed last row match the thermo grid, the temperatures do not
  err = read_chemistry(folder='ranges/rates-step2', mech_name=mech_name)
  call verdict('rate table on a 2 K step (1400..1800 K, 201 rows): refused with ios = 6', err == 6)
  call free_chemistry_data()
  call execute_command_line(trim(self)//' child-grid > ranges/child-grid.out 2> ranges/child-grid.err', exitstat=err)
  call execute_command_line('test "$(command grep -c -F ''[ERROR] FLINT read_'' ranges/child-grid.err)" = 5', exitstat=err)
  call verdict('child process: the five grid refusals are on the error unit too', err == 0)

  ! thermo tables on a 2 K step (last: the phase arrays are reloaded)
  if (allocated(species_names)) deallocate(species_names)
  if (allocated(wm_tab)) deallocate(wm_tab)
  if (allocated(Ri_tab)) deallocate(Ri_tab)
  if (allocated(h_tab)) deallocate(h_tab)
  if (allocated(s_tab)) deallocate(s_tab)
  if (allocated(cp_tab)) deallocate(cp_tab)
  if (allocated(dcpi_tab)) deallocate(dcpi_tab)
  err = read_idealgas_thermo('ranges/thermo-step2/')
  call verdict('thermo table on a 2 K step (1400..1800 K): refused with ios = 4', err == 4)

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
end program test
