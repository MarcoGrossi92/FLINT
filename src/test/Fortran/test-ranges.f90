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
