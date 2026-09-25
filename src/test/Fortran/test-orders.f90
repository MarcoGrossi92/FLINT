! test-orders: the general procedure with the optional 'Reaction orders' block of chemistry-info.txt
! reproduces the FORWARD production rates of Cantera for a mechanism with explicit yaml orders (the reverse rate
! constants of the table come from the thermo database chosen by the table writer, not from the yaml: the driver
! zeroes kb_tab, so the comparison isolates the reaction orders, which act on the forward rate)
! (JLR-frassoldati, reaction 1: CH4^0.5 O2^1.3), and without the block keeps the integer-rounded
! stoichiometric law (old INPUT folders unchanged). Fixture test/orders/JLR-frassoldati made by
! test/orders/make_JLR-frassoldati.py from the tables of a table writer (1200..1800 K) with Cantera references
! at three states. Needs no Cantera at run time. Exit code 1 on failure.
program test
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_chemistry
  use FLINT_Lib_Chemistry_data
  use FLINT_Lib_Chemistry_wdot
  implicit none
  character(32) :: mech_name
  character(len=512) :: line
  integer :: err, u, nsr, nstate, k, i
  real(8) :: T, tol
  real(8), allocatable :: roi(:), w(:), w_orders(:), w_nint(:), roi0(:)
  integer :: nfail
  real(8), parameter :: wm_ref(9) = [31.998d0, 16.043d0, 18.015d0, 28.010d0, 44.009d0, 2.016d0, 1.008d0, 15.999d0, 17.007d0]
  character(len=s_str_len), parameter :: names_ref(9) = [character(len=s_str_len) :: 'O2', 'CH4', 'H2O', 'CO', 'CO2', 'H2', 'H', 'O', 'OH']

  nfail = 0; tol = 1d-8
  ! species of phase.txt (no thermo tables are needed: general is called directly)
  ns = 9
  allocate(wm_tab(ns), Ri_tab(ns), species_names(ns))
  wm_tab = wm_ref; Ri_tab = Runiv/wm_tab; species_names = names_ref
  allocate(roi(ns), w(ns), w_orders(ns), w_nint(ns), roi0(ns))

  ! helper
  call verdict('pow_order(0, -0.75) = 0 (Cantera: zero rate at zero concentration)', pow_order(0d0, -0.75d0) == 0d0)
  call verdict('pow_order(2, 2) = 4 (integer path)', pow_order(2d0, 2d0) == 4d0)
  call verdict('pow_order(2, 0.5) = sqrt(2) (real path)', pow_order(2d0, 0.5d0) == sqrt(2d0))
  call verdict('pow_order(-1, 0.5) = 0 (negative concentration clipped)', pow_order(-1d0, 0.5d0) == 0d0)

  ! 1) with the block: Cantera with the yaml orders
  err = read_chemistry(folder='orders/JLR-frassoldati', mech_name=mech_name)
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_chemistry orders/JLR-frassoldati: ios=', err; stop 1; endif
  call verdict('block read: have_orders', have_orders)
  call verdict('block read: orders of reaction 1 = CH4 0.5, O2 1.3, others stoichiometric', &
    ord_arrh_tab(2,1) == 0.5d0 .and. ord_arrh_tab(1,1) == 1.3d0 .and. ord_arrh_tab(1,2) == ni1_arrh_tab(1,2))
  call Assign_Mechanism(mech_name)   ! JLR-Frassoldati is not hooked: general (WARNING expected)
  kb_tab = 0d0   ! forward part only: the table's kb comes from the writer's thermo database, not from the yaml
  open(newunit=u, file='orders/JLR-frassoldati/reference.txt', status='old', action='read')
  read(u,'(A)') line
  read(u,*) nsr, nstate
  if (nsr /= ns) then; write(*,'(A)') '[FAIL] reference species count'; stop 1; endif
  do k = 1, nstate
    read(u,*) T; read(u,*) roi0; read(u,*) w_orders; read(u,*) w_nint
    roi = roi0; w = 0d0
    call general(roi, T, w)
    write(*,'(A,F7.1,A,ES10.3,A,ES10.3)') ' T = ', T, ' K: max |w - Cantera(orders)| / max|w| = ', &
      maxval(abs(w - w_orders))/maxval(abs(w_orders)), ', vs nint law = ', maxval(abs(w - w_nint))/maxval(abs(w_orders))
    call verdict('with the block: general = Cantera forward rates with the yaml orders', maxval(abs(w - w_orders)) <= tol*maxval(abs(w_orders)))
    call verdict('with the block: general differs from the integer-rounded law', maxval(abs(w - w_nint)) > 1d-3*maxval(abs(w_orders)))
  enddo
  close(u)
  call free_chemistry_data()

  ! 2) without the block (old INPUT): the integer-rounded stoichiometric law, unchanged
  call execute_command_line('mkdir -p orders/noblock && cp orders/JLR-frassoldati/chemistry-Arrhenius.dat orders/noblock/ && ' // &
    'cp orders/JLR-frassoldati/chemistry-info-noblock.txt orders/noblock/chemistry-info.txt')
  err = read_chemistry(folder='orders/noblock', mech_name=mech_name)
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_chemistry orders/noblock: ios=', err; stop 1; endif
  call verdict('no block: have_orders is false', .not. have_orders)
  kb_tab = 0d0
  open(newunit=u, file='orders/JLR-frassoldati/reference.txt', status='old', action='read')
  read(u,'(A)') line; read(u,*) nsr, nstate
  do k = 1, nstate
    read(u,*) T; read(u,*) roi0; read(u,*) w_orders; read(u,*) w_nint
    roi = roi0; w = 0d0
    call general(roi, T, w)
    call verdict('no block: general = integer-rounded stoichiometric law (Cantera forward rate constants)', &
      maxval(abs(w - w_nint)) <= tol*maxval(abs(w_nint)))
  enddo
  close(u)
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
