! test-stoich: the general procedure reproduces Cantera's net production rates for reactions with FRACTIONAL
! stoichiometric coefficients and no explicit orders: the mass-action law with the real coefficients (forward
! [H2] [O2]^0.5 for H2 + 0.5 O2, reverse [H2] [O2]^0.5 for the fractional products of H2O <=> H2 + 0.5 O2,
! kb = kf/Kc of the same reaction), for Arrhenius, three-body, Troe and Lindemann reactions, and its net rates
! vanish at Cantera's equilibrium composition. The integer-rounded law used before by the Arrhenius loop
! (nint: [O2]^1 for 0.5 O2, [O2]^0 for 0.25 O2, [H2O]^2 for 1.5 H2O) fails both. stoich-int (2 H2 + O2 <=>
! 2 H2O) is the integer control. Fixtures test/stoich/<name> made by test/stoich/make_stoich.py (constructed
! yaml mechanisms, tables of a table writer from the thermo of the yaml, Cantera references at table nodes, lean, rich,
! near-equilibrium and equilibrium states). Also checks that general has no analytical Jacobian (the
! integrator differentiates the rates numerically). Needs no Cantera at run time. Exit code 1 on failure.
program test
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_chemistry
  use FLINT_Lib_Chemistry_data
  use FLINT_Lib_Chemistry_wdot
  implicit none
  integer, parameter :: ncase = 6
  character(len=16), parameter :: cases(ncase) = [character(len=16) :: 'stoich-frac', 'stoich-prod', 'stoich-3b', &
    'stoich-troe', 'stoich-lind', 'stoich-int']
  real(8), parameter :: tol = 1d-10
  character(32) :: mech_name
  character(len=512) :: line
  integer :: c, k, u, err, nstate, flag, nfail, kworst
  real(8) :: T, pres, rerr, emax, eqmax
  real(8), allocatable :: roi(:), w(:), wct(:), gross(:)

  nfail = 0
  do c = 1, ncase
    open(newunit=u, file='stoich/'//trim(cases(c))//'/reference.txt', status='old', action='read')
    read(u,'(A)') line
    read(u,*) ns, nstate
    if (allocated(wm_tab)) deallocate(wm_tab)
    if (allocated(Ri_tab)) deallocate(Ri_tab)
    if (allocated(species_names)) deallocate(species_names)
    allocate(wm_tab(ns), Ri_tab(ns), species_names(ns), roi(ns), w(ns), wct(ns), gross(ns))
    read(u,*) species_names
    read(u,*) wm_tab
    Ri_tab = Runiv/wm_tab
    err = read_chemistry(folder='stoich/'//trim(cases(c)), mech_name=mech_name)
    if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_chemistry stoich/'//trim(cases(c))//': ios=', err; stop 1; endif
    call Assign_Mechanism(mech_name)   ! not hooked: general (WARNING expected)
    call verdict(trim(cases(c))//': the general procedure is selected', associated(chemistry_source, general))
    call verdict(trim(cases(c))//': general has no analytical Jacobian', .not. associated(chemistry_jacobian))
    emax = 0d0; eqmax = 0d0; kworst = 0
    do k = 1, nstate
      read(u,*) T, pres, flag
      read(u,*) roi
      read(u,*) wct
      read(u,*) gross
      w = 0d0
      call chemistry_source(roi, T, w)
      rerr = maxval(abs(w - wct))/maxval(gross)
      if (rerr > emax) then; emax = rerr; kworst = k; endif
      if (flag == 1) eqmax = max(eqmax, maxval(abs(w))/maxval(gross))
    enddo
    close(u)
    write(*,'(A,I0,A,ES9.2,A,I0,A,ES9.2)') ' '//trim(cases(c))//': ', nstate, ' states, max|w - Cantera|/max(gross) = ', &
      emax, ' (state ', kworst, '), at Cantera''s equilibrium max|w|/max(gross) = ', eqmax
    call verdict(trim(cases(c))//': general = Cantera net production rates (real stoichiometric coefficients)', emax <= tol)
    call verdict(trim(cases(c))//': general net rates vanish at Cantera''s equilibrium composition', eqmax <= tol)
    call free_chemistry_data()
    deallocate(roi, w, wct, gross)
  enddo
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
