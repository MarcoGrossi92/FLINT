! test-andersen: the WD-Andersen routine (Westbrook-Dryer steps 1-2, step 3 = CO2 dissociation with the
! Andersen orders [CO2] [H2O]^0.5 [O2]^-0.25) reproduces the net production rates of Cantera on the
! tables written by a table writer from WD-Andersen.yaml (fixture test/andersen/WD-Andersen: a complete INPUT
! folder on one
! grid, 1100..1400 K, made by test/andersen/make_WD-Andersen.py with the references embedded) at random
! states on the table nodes (relative 1e-12), including states with a zero concentration of a species
! with a negative order (O2 = 0: the rate of step 3 is zero, the convention of Cantera, and not
! 0**(-0.25) = +Infinity); and the same zero-concentration convention at the other hand-written site
! with a negative order (Coronetti with H2 = 0: finite source, no divide-by-zero exception).
! Needs no Cantera at run time. Exit code 1 on failure.
program test
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_get_flag, ieee_set_flag, ieee_divide_by_zero
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  use FLINT_Load_chemistry
  use FLINT_Lib_Chemistry_data
  use FLINT_Lib_Chemistry_wdot
  use WD_mod, only: Andersen
  use coronetti_mod, only: Coronetti
  implicit none
  character(32) :: mech_name
  character(len=256) :: line
  character(len=16) :: tag
  integer :: err, u, nsr, nstate, k, it, ir, nfail
  real(8) :: T, rel, tol
  real(8), allocatable :: roi(:), w(:), wref(:)
  logical :: flag, fin
  integer(8) :: seed

  nfail = 0; tol = 1d-12
  ! 1) the fixture, loaded as a solver does: thermo, chemistry, then Assign_Mechanism (contract check)
  err = read_idealgas_thermo('andersen/WD-Andersen/')
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_idealgas_thermo andersen/WD-Andersen: ios=', err; stop 1; endif
  err = read_chemistry(folder='andersen/WD-Andersen', mech_name=mech_name)
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_chemistry andersen/WD-Andersen: ios=', err; stop 1; endif
  call verdict('mechanism name read from the fixture tables: WD-Andersen', trim(mech_name) == 'WD-Andersen')
  call Assign_Mechanism(mech_name)
  call verdict('Assign_Mechanism hooks the Andersen routine (contract check passed on the fixture layout)', &
    associated(chemistry_source, Andersen))
  allocate(roi(ns), w(ns), wref(ns))
  open(newunit=u, file='andersen/WD-Andersen/reference.txt', status='old', action='read')
  read(u,'(A)') line
  read(u,*) nsr, nstate
  if (nsr /= ns) then; write(*,'(A)') '[FAIL] reference species count'; stop 1; endif
  do k = 1, nstate
    read(u,*) T, tag; read(u,*) roi; read(u,*) wref
    w = 0d0
    call chemistry_source(roi, T, w)
    fin = all(ieee_is_finite(w))
    if (maxval(abs(wref)) > 0d0) then
      rel = maxval(abs(w - wref))/maxval(abs(wref))
      write(*,'(A,F7.1,A,A,A,ES10.3)') ' T = ', T, ' K ', trim(tag), ': max |w - Cantera| / max |Cantera| = ', rel
      call verdict('Andersen = Cantera net production rates, '//trim(tag), fin .and. rel <= tol)
    else
      write(*,'(A,F7.1,A,A,A,ES10.3)') ' T = ', T, ' K ', trim(tag), ': Cantera rates are 0, max |w| = ', maxval(abs(w))
      call verdict('Andersen = 0 where Cantera gives 0, '//trim(tag)//' (no WD step proceeds without O2; step 3 not +Inf)', &
        fin .and. maxval(abs(w)) == 0d0)
    endif
  enddo
  close(u)

  ! 2) the other hand-written site with a negative order: Coronetti (its 9 slots O2, C4H6, H2O, CO, CO2, H2, O, H, OH)
  !    with H2 = 0 on synthetic 1 K tables: the reverse term of H2 + 1/2 O2 <-> H2O is 0, 0**(-0.75) not evaluated
  call free_chemistry_data()
  deallocate(wm_tab, Ri_tab, species_names, h_tab, cp_tab, dcpi_tab, s_tab)
  ns = 9
  allocate(wm_tab(ns), Ri_tab(ns))
  wm_tab = [31.998d0, 54.092d0, 18.015d0, 28.010d0, 44.009d0, 2.016d0, 15.999d0, 1.008d0, 17.007d0]
  Ri_tab = 8314.46d0/wm_tab
  allocate(kf_tab(1:3000, 1:10), kb_tab(1:3000, 1:10)); nrc_arrh = 10; T_tab_min = 1; T_tab_max = 3000
  seed = 20260925_8
  do ir = 1, 10
    do it = 1, 3000
      kf_tab(it, ir) = 10d0**(6d0*urand(seed) - 3d0)
      kb_tab(it, ir) = 10d0**(6d0*urand(seed) - 3d0)
    enddo
  enddo
  deallocate(roi, w); allocate(roi(ns), w(ns))
  roi = [2d-2, 1d-2, 3d-3, 4d-3, 5d-4, 0d0, 1d-4, 1d-5, 2d-4]
  call ieee_set_flag(ieee_divide_by_zero, .false.)
  w = 0d0
  call Coronetti(roi, 1500.5d0, w)
  call ieee_get_flag(ieee_divide_by_zero, flag)
  call verdict('Coronetti with H2 = 0: finite source (reverse term of H2 + 1/2 O2 <-> H2O is 0, not 0**(-0.75))', &
    all(ieee_is_finite(w)))
  call verdict('Coronetti with H2 = 0: no divide-by-zero exception (a trap under -fpe0 / -ffpe-trap=zero)', .not. flag)

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
  function urand(st) result(x)
    integer(8), intent(inout) :: st
    real(8) :: x
    st = mod(1103515245_8*st + 12345_8, 2147483648_8)
    x = dble(st)/2147483648d0
  end function urand
end program test
