!> test-runiv: the universal gas constant of FLINT is the exact SI value R = N_A k_B = 8314.46261815324
!> J/(kmol K) (SI 2019 / CODATA 2018 exact value; Cantera 3.0.1 gas_constant) in both places that define
!> it (Runiv of FLINT_Lib_Thermodynamic, Rr of FLINT_CEA_data), the ideal-gas loader derives Ri_tab from
!> it, and the quantities that carry it agree with Cantera: the pressure p = rho R_mix T of a Cantera state
!> of database/WD, and the compiled Frolov routine (global-H2.f90: rate proportional to (p/p_atm)^-1.15 with
!> p = sum(roi*Ri_tab)*T) on a Cantera H2/O2/N2 state at 2000 K and 10 atm, against the same law evaluated
!> at the Cantera pressure. With the former Runiv = 8314.51 every FLINT pressure was 5.7e-6 high and the
!> Frolov rate 6.6e-6 low. Run from test/thermo: ../../bin/test/test-runiv. Exit code 1 on any [FAIL].
program test_runiv
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  use FLINT_CEA_data, only: Rr
  use globH2_mod, only: Frolov
  implicit none
  real(8), parameter :: R_SI = 8314.46261815324d0
  ! Cantera 3.0.1 on database/WD/WD.yaml: TPX = 1500 K, 101325 Pa, CH4:1 O2:2 (species order of phase.txt)
  real(8), parameter :: T_ref = 1500d0, p_ref = 101325d0, rho_ref = 2.16756219385544596d-1
  real(8), parameter :: Y_ref(5) = [2.00439785604517778d-1, 7.99560214395482194d-1, 0d0, 0d0, 0d0]
  ! Cantera 3.0.1, species in the slots of the Frolov routine (1 O2, 2 H2O, 3 H2, 4 N2):
  ! TPX = 2000 K, 10 atm, H2:2 O2:1 N2:3.76
  real(8), parameter :: T_fr = 2000d0, p_fr = 1013250d0, rho_fr = 1.27420816282737959d0
  real(8), parameter :: Y_fr(4) = [2.26354006971007327d-1, 0d0, 2.85223875275673958d-2, 7.45123605501425201d-1]
  real(8), parameter :: W_fr(4) = [31.998d0, 18.015d0, 2.016d0, 28.014d0]
  real(8), allocatable :: rhoi(:)
  real(8) :: p, f, roi4(4), w4(4), law, dev
  integer :: err, nok, nfail

  nok = 0; nfail = 0
  call verdict('Runiv = 8314.46261815324 J/(kmol K) (exact SI value, Cantera 3.0.1)', Runiv == R_SI)
  call verdict('CEA Rr = Runiv', Rr == Runiv)
  err = read_idealgas_thermo('../../database/WD/')
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_idealgas_thermo database/WD: ios=', err; stop 1; endif
  call verdict('database/WD loaded with 5 species', ns == 5)
  call verdict('Ri_tab*wm_tab = Runiv for every species (1e-14)', maxval(abs(Ri_tab*wm_tab/Runiv - 1d0)) < 1d-14)
  allocate(rhoi(ns)); rhoi = rho_ref*Y_ref
  p = sum(rhoi)*f_Rtot(rhoi)*T_ref
  write(*,'(A,ES23.15,A,ES10.3)') 'p(rho, T, Y) = ', p, ' Pa, relative deviation from Cantera ', abs(p/p_ref - 1d0)
  call verdict('p(rho, T, Y) = Cantera p within 1e-9 (5.7e-6 with 8314.51)', abs(p/p_ref - 1d0) < 1d-9)
  f = (p/101325d0)**(-1.15d0)
  write(*,'(A,ES23.15,A,ES10.3)') '(p/p_atm)^-1.15 = ', f, ', deviation from 1 ', abs(f - 1d0)
  call verdict('(p/p_atm)^-1.15 = 1 within 1e-9 at the Cantera state (6.6e-6 with 8314.51)', abs(f - 1d0) < 1d-9)

  ! the compiled Frolov routine on its own species layout
  deallocate(wm_tab, Ri_tab); ns = 4
  allocate(wm_tab(ns), Ri_tab(ns)); wm_tab = W_fr; Ri_tab = Runiv/wm_tab
  roi4 = rho_fr*Y_fr
  call Frolov(roi4, T_fr, w4)
  law = -0.5d0*8d11*(p_fr/101325d0)**(-1.15d0)*(roi4(3)/W_fr(3))**2*(roi4(1)/W_fr(1))*exp(-1d4/T_fr)
  dev = max(abs(w4(1)/(W_fr(1)*law) - 1d0), abs(w4(2)/(-2d0*W_fr(2)*law) - 1d0), abs(w4(3)/(2d0*W_fr(3)*law) - 1d0))
  write(*,'(A,ES23.15,A,ES10.3)') 'Frolov omegadot(O2) = ', w4(1), ' kg/(m3 s), deviation from the law ', dev
  call verdict('Frolov routine = the p^-1.15 law at the Cantera p within 1e-9 (6.6e-6 with 8314.51)', &
               dev < 1d-9 .and. w4(4) == 0d0)
  write(*,'(A,I0,A,I0,A)') 'test-runiv: ', nok, ' ok / ', nfail, ' FAIL'
  if (nfail > 0) stop 1
contains
  subroutine verdict(what, ok)
    character(len=*), intent(in) :: what
    logical, intent(in) :: ok
    if (ok) then
      write(*,'(A)') '[ok]   '//what; nok = nok + 1
    else
      write(*,'(A)') '[FAIL] '//what; nfail = nfail + 1
    endif
  end subroutine verdict
end program test_runiv
