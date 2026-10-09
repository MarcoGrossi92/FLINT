! test-inert: species appended after the slots of a compiled mechanism routine are inert in EVERY
! code path, on the fixture test/chemistry/inert/WD-plus2 (database/WD + N2 + AR, made by
! test/chemistry/inert/make_WD-plus2.py). The omegadot/dwdr/dwdT dummies of the routines are INTENT(OUT):
! the standard leaves them undefined on entry, so a routine that assigns only its own slots leaves
! the appended ones undefined whatever the caller zeroed before the call. Every check fills the
! output with a sentinel value first; a routine that defines its whole block leaves no sentinel.
! Meant to run on optimised builds too (the solvers' -O3 flags). Exit code 1 on failure.
program test
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  use FLINT_Load_chemistry
  use FLINT_Lib_Chemistry_data
  use FLINT_Lib_Chemistry_wdot
  use FLINT_Lib_Chemistry_rhs
  use globH2_mod
  implicit none
  integer, parameter :: ns_r = 5
  real(8), parameter :: sentinel = 7.0d0, T = 1500.5d0
  real(8), allocatable :: roi(:), w(:), Z(:), F(:), DFY(:,:), dwdr(:,:), dwdT(:)
  real(8) :: rpar(1)
  integer :: ipar(1), err, nz, nfail
  character(32) :: mech_name

  nfail = 0
  err = read_idealgas_thermo('inert/WD-plus2/')
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_idealgas_thermo inert/WD-plus2: ios=', err; stop 1; endif
  err = read_chemistry(folder='inert/WD-plus2', mech_name=mech_name)
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_chemistry inert/WD-plus2: ios=', err; stop 1; endif
  call Assign_Mechanism(mech_name)
  call verdict('WD-plus2: 7 species loaded = 5 routine slots + 2 appended, WD hooked', &
    ns == ns_r + 2 .and. associated(chemistry_source))
  allocate(roi(ns), w(ns), Z(ns+1), F(ns+1), DFY(ns+1,ns+1), dwdr(ns,ns), dwdT(ns))
  roi = [0.05d0, 0.30d0, 0.10d0, 0.08d0, 0.12d0, 0.40d0, 0.02d0]

  ! 1) the routine called directly with a sentinel-filled omegadot
  w = sentinel
  call chemistry_source(roi, T, w)
  call verdict('WD routine, direct call: the two appended slots are exactly 0', all(w(ns_r+1:ns) == 0d0))
  call verdict('WD routine, direct call: no sentinel left in any slot', all(w /= sentinel))
  call verdict('WD routine, direct call: own slots finite and reacting', all(w(1:ns_r) == w(1:ns_r)) .and. any(w(1:ns_r) /= 0d0))

  ! 2) the rhs used by the ODE integrator
  nz = ns + 1; Z(1:ns) = roi; Z(nz) = T
  F = sentinel
  call rhs_native(nz, 0d0, Z, F)
  call verdict('rhs_native: F of the two appended species is exactly 0', all(F(ns_r+1:ns) == 0d0))
  call verdict('rhs_native: no sentinel left (species and temperature rows)', all(F /= sentinel))

  ! 3) the analytical Jacobian path on the same 7-species phase: Frolov_nopressure has 3 slots
  chemistry_source => Frolov_nopressure
  chemistry_jacobian => Frolov_nopressure_jac
  w = sentinel
  call chemistry_source(roi, T, w)
  call verdict('Frolov_nopressure, direct call: slots 4..7 are exactly 0', all(w(4:ns) == 0d0))
  dwdr = sentinel; dwdT = sentinel
  call chemistry_jacobian(roi, T, dwdr, dwdT)
  call verdict('Frolov_nopressure_jac, direct call: rows 4..7 of dwdr and dwdT are exactly 0', &
    all(dwdr(4:ns,:) == 0d0) .and. all(dwdT(4:ns) == 0d0))
  call verdict('Frolov_nopressure_jac, direct call: no sentinel left', all(dwdr /= sentinel) .and. all(dwdT /= sentinel))
  DFY = sentinel
  call jac_native(nz, 0d0, Z, DFY, nz, rpar, ipar)
  call verdict('jac_native: rows 4..7 (species after the routine slots) are exactly 0', all(DFY(4:ns,:) == 0d0))
  call verdict('jac_native: no sentinel left', all(DFY /= sentinel))

  call free_chemistry_data()
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
