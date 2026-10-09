! test-rhs-range: rhs_native and jac_native bail out (F = -1, zero Jacobian) when the temperature is
! outside the range of the RATE tables (T_tab_min..T_tab_max), as they already did outside the thermo
! tables. read_chemistry refuses rate tables on another grid than the thermo one, so for tables loaded
! through it the two guards coincide; the rate guard is a defence for tables set by another path (a
! driver's in-memory tables, a loader of another format): below the first row of the rate tables the
! source term used to be an out-of-bounds read (garbage in RELEASE, runtime error under -fcheck=all).
! Fixture: test/chemistry/tables/WD-100K (thermo and rate tables 100..400 K, one grid); the rate range is then
! narrowed in memory to 200..300 K. Exit code 1 on failure; compiles against FLINT <= 2223136.
program test
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  use FLINT_Load_chemistry
  use FLINT_Lib_Chemistry_data
  use FLINT_Lib_Chemistry_wdot
  use FLINT_Lib_Chemistry_rhs
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  character(32) :: mech_name
  integer :: err, nz, nfail
  real(8), allocatable :: Z(:), F(:), DFY(:,:)
  real(8) :: rpar(1)
  integer :: ipar(1)

  nfail = 0; rpar = 0.d0; ipar = 0
  err = read_idealgas_thermo('tables/WD-100K/')
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_idealgas_thermo tables/WD-100K: ios=', err; stop 1; endif
  err = read_chemistry(folder='tables/WD-100K', mech_name=mech_name)
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_chemistry tables/WD-100K: ios=', err; stop 1; endif
  call Assign_Mechanism(mech_name)
  nz = ns + 1
  allocate(Z(nz), F(nz), DFY(nz, nz))
  Z(1:ns) = 1d-20; Z(1) = 0.2d0; Z(2) = 0.8d0
  write(*,'(A,I0,A,I0,A,I0,A,I0,A)') ' thermo tables ', merge(1, Tmin, Tmin == 0), '..', Tmax, ' K; rate tables ', &
    lbound(kf_tab, dim=1), '..', ubound(kf_tab, dim=1), ' K'
  call verdict('tables/WD-100K loaded on one grid: rate range 100..400 K', T_tab_min == 100 .and. T_tab_max == 400)
  call check_T(50.d0,  .true.,  '50 K (below the thermo and the rate tables)')
  call check_T(250.d0, .false., '250 K (inside both)')
  call check_T(ieee_value(1d0, ieee_quiet_nan), .true., 'NaN (bit test: survives -ffast-math and FPE traps)')
  ! the rate range narrowed in memory (tables set by another path): 200..300 K inside the thermo tables
  T_tab_min = 200; T_tab_max = 300
  call check_T(150.d0,  .true.,  '150 K (inside the thermo tables, below the rate range)')
  call check_T(350.d0,  .true.,  '350 K (inside the thermo tables, above the rate range)')
  call check_T(299.5d0, .false., '299.5 K (inside both)')
  call free_chemistry_data()
  if (nfail > 0) then
    write(*,'(A,I0,A)') ' Verdict -> fail (', nfail, ' checks)'
    stop 1
  endif
  write(*,'(A)') ' Verdict -> pass'
contains
  subroutine check_T(T, outside, what)
    real(8), intent(in) :: T
    logical, intent(in) :: outside
    character(*), intent(in) :: what
    Z(nz) = T
    F = 7.d0
    call rhs_native(nz, 0.d0, Z, F)
    if (outside) then
      call verdict(what//': rhs_native returns F = -1', all(F == -1.d0))
      DFY = 7.d0
      call jac_native(nz, 0.d0, Z, DFY, nz, rpar, ipar)
      call verdict(what//': jac_native returns a zero Jacobian', all(DFY == 0.d0))
    else
      call verdict(what//': rhs_native returns finite rates, not the bail-out value', &
        all(F == F) .and. .not. all(F == -1.d0) .and. F(1) < 0.d0)
    endif
  end subroutine check_T
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
