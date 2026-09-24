! test-tables: the rate tables are indexed by temperature (row T = rate at T kelvin) whatever the
! first temperature of the table.
!  1. database/WD (first row 1 K): f_kf/f_kb return row T at every T (regression guard);
!  2. the same rows restricted to 100..400 K (test/tables/WD-100K, made by test/tables/make_WD-100K.py):
!     row T is still the rate at T kelvin (f_kf/f_kb == the full 1 K table at the same T);
!  3. every hand-written routine (WD, Andersen, OSK, JLR, Frassoldati, CKJLR10sp, singh, Singh_WC32,
!     singhC3H6, Coronetti, Nassini_4, Frolov_nopressure, Frolov) with synthetic in-memory tables:
!     omegadot with the tables starting at 50, 100, 300 and 799 K is BIT-IDENTICAL to omegadot with
!     the same rows in tables starting at 1 K;
!  4. positive control: the accessor of FLINT <= 2223136 (assumed-shape dummy tab(:,:), copied below
!     as old_comp_ch_tabT) returns the rate of row T + Tmin - 1 when the table starts at Tmin /= 1;
!  5. the public comp_ch_tabT of the library (dummy tab(T_tab_min:,:)) equals f_kf/f_kb at every row of
!     the 1 K table (where it also equals the old accessor), of test/tables/WD-100K and of the synthetic
!     tables starting at 100 K (where the old accessor differs in every sample).
! Needs no Cantera. Exit code 1 on failure.
program test
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  use FLINT_Load_chemistry
  use FLINT_Lib_Chemistry_data
  use WD_mod
  use JLRs_mod
  use singh_mod
  use coronetti_mod
  use globH2_mod
  implicit none
  integer, parameter :: T0 = 100, T1 = 400
  integer, parameter :: nsyn = 12, nrc_syn = 12, Tmax_syn = 3000, nroutine = 13, nsamp = 3, ntmin = 4
  integer, parameter :: Tmins(ntmin) = [50, 100, 300, 799]
  real(8), parameter :: Tsamp(nsamp) = [800.0d0, 1234.5d0, 2998.37d0]
  real(8), parameter :: Tctrl(nsamp) = [800.0d0, 1234.5d0, 2500.37d0]
  character(len=18), parameter :: rname(nroutine) = [character(len=18) :: 'WD', 'Andersen', 'OSK', 'JLR', &
    'Frassoldati', 'CKJLR10sp', 'singh', 'Singh_WC32', 'singhC3H6', 'Coronetti', 'Nassini_4', &
    'Frolov_nopressure', 'Frolov']
  real(8), allocatable :: kf_full(:,:), kb_full(:,:), kf_ref(:,:), kb_ref(:,:)
  real(8) :: Tdiff(2), a, b, c, dmax, td
  real(8) :: roi0(nsyn), roi(nsyn), w(nsyn), w_ref(nsyn, nsamp)
  integer :: err, ir, T, k, m, Tint(2), nfail1, nfail2, nfail3, nctrl, nctrl_tot, lb, ub, nbad
  integer :: nctrl1, nctrl2, nctrl2_tot, nrest1, nrest2, nrest4
  character(32) :: mech_name

  Tdiff = [0.d0, 0.37d0]
  err = read_idealgas_thermo('../database/WD/')
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_idealgas_thermo ../database/WD: ios=', err; stop 1; endif
  err = read_chemistry(folder='../database/WD/', mech_name=mech_name)
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_chemistry ../database/WD: ios=', err; stop 1; endif
  lb = lbound(kf_tab, dim=1); ub = ubound(kf_tab, dim=1)
  write(*,'(A,I0,A,I0,A,I0)') ' WD full table: rows ', lb, '..', ub, ' K, nrc_arrh = ', nrc_arrh

  ! 1) table starting at 1 K: f_kf/f_kb must return row T (the rate at T kelvin) at every T;
  !    the assumed-shape accessor of FLINT <= 2223136 (old_comp_ch_tabT below) is right here too
  nfail1 = 0; nctrl1 = 0; nrest1 = 0
  do ir = 1, nrc_arrh
    do T = lb, ub-1
      Tint = [T, T+1]
      do k = 1, 2
        a = kf_tab(T,ir) + (kf_tab(T+1,ir) - kf_tab(T,ir))*Tdiff(k); b = f_kf(ir, Tint, Tdiff(k))
        if (a /= b) nfail1 = nfail1 + 1
        if (old_comp_ch_tabT(ir, kf_tab, Tint, Tdiff(k)) /= b) nctrl1 = nctrl1 + 1
        if (comp_ch_tabT(ir, kf_tab, Tint, Tdiff(k)) /= b) nrest1 = nrest1 + 1
        a = kb_tab(T,ir) + (kb_tab(T+1,ir) - kb_tab(T,ir))*Tdiff(k); b = f_kb(ir, Tint, Tdiff(k))
        if (a /= b) nfail1 = nfail1 + 1
        if (old_comp_ch_tabT(ir, kb_tab, Tint, Tdiff(k)) /= b) nctrl1 = nctrl1 + 1
        if (comp_ch_tabT(ir, kb_tab, Tint, Tdiff(k)) /= b) nrest1 = nrest1 + 1
      enddo
    enddo
  enddo
  write(*,'(A,I0,A,I0,A,I0)') ' 1 K table   : f_kf|f_kb /= row T mismatches = ', nfail1, &
    '; old assumed-shape accessor /= f_kf|f_kb: ', nctrl1, '; comp_ch_tabT /= f_kf|f_kb: ', nrest1
  allocate(kf_full(T0:T1, nrc_arrh), kb_full(T0:T1, nrc_arrh))
  kf_full = kf_tab(T0:T1, :); kb_full = kb_tab(T0:T1, :)
  call free_chemistry_data()

  ! 2) the same rows in a table starting at 100 K: row T must still be the rate at T kelvin
  err = read_chemistry(folder='tables/WD-100K', mech_name=mech_name)
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_chemistry tables/WD-100K: ios=', err; stop 1; endif
  lb = lbound(kf_tab, dim=1); ub = ubound(kf_tab, dim=1)
  write(*,'(A,I0,A,I0)') ' WD-100K table: rows ', lb, '..', ub, ' K'
  if (lb /= T0 .or. ub /= T1) then; write(*,'(A)') '[FAIL] unexpected table bounds'; stop 1; endif
  nfail2 = 0; dmax = 0.d0; nctrl2 = 0; nctrl2_tot = 0; nrest2 = 0
  do ir = 1, nrc_arrh
    do T = T0, T1-1
      Tint = [T, T+1]
      do k = 1, 2
        b = f_kf(ir, Tint, Tdiff(k))
        c = kf_full(T,ir) + (kf_full(T+1,ir) - kf_full(T,ir))*Tdiff(k)
        if (b /= c) nfail2 = nfail2 + 1
        ! the old accessor renumbers the 301 rows from 1: row T is T + 99 K, and T > 300 is out of
        ! bounds (a bounds-checked build aborts there; an optimised one reads past the table)
        if (T + 1 <= T1 - T0 + 1) then
          a = old_comp_ch_tabT(ir, kf_tab, Tint, Tdiff(k)); nctrl2_tot = nctrl2_tot + 1
          if (a /= b) then
            nctrl2 = nctrl2 + 1
            if (b /= 0.d0) dmax = max(dmax, abs(a-b)/abs(b))
          endif
        endif
        if (comp_ch_tabT(ir, kf_tab, Tint, Tdiff(k)) /= f_kf(ir, Tint, Tdiff(k))) nrest2 = nrest2 + 1
        b = f_kb(ir, Tint, Tdiff(k))
        c = kb_full(T,ir) + (kb_full(T+1,ir) - kb_full(T,ir))*Tdiff(k)
        if (b /= c) nfail2 = nfail2 + 1
        if (comp_ch_tabT(ir, kb_tab, Tint, Tdiff(k)) /= b) nrest2 = nrest2 + 1
        ! (no kb control: the three WD reactions are irreversible, kb = 0 on every row)
      enddo
    enddo
  enddo
  write(*,'(A,I0,A,I0,A,I0,A,ES10.2)') ' 100 K table : f_kf|f_kb /= full 1 K table mismatches = ', nfail2, &
    '; old assumed-shape accessor (kf, rows 100..300) /= f_kf: ', nctrl2, ' of ', nctrl2_tot, ', max rel. error ', dmax
  write(*,'(A,I0,A)') ' 100 K table : comp_ch_tabT /= f_kf|f_kb at every row: ', nrest2, ' (expected 0)'
  call free_chemistry_data()

  ! 3) hand-written routines with synthetic tables: bit identity for tables starting at Tmin /= 1
  ns = nsyn
  if (allocated(wm_tab)) deallocate(wm_tab)
  if (allocated(Ri_tab)) deallocate(Ri_tab)
  allocate(wm_tab(ns), Ri_tab(ns))
  wm_tab = [2.016d0, 31.998d0, 18.015d0, 28.010d0, 44.009d0, 16.043d0, 17.007d0, 1.008d0, 15.999d0, &
            168.3d0, 42.08d0, 28.014d0]
  Ri_tab = 8314.46d0/wm_tab
  allocate(kf_ref(1:Tmax_syn, nrc_syn), kb_ref(1:Tmax_syn, nrc_syn))
  do ir = 1, nrc_syn
    do T = 1, Tmax_syn
      kf_ref(T,ir) = 1d3*(1d0 + 0.1d0*ir)*sqrt(dble(T))*exp(-1500d0*ir/dble(T))
      kb_ref(T,ir) = kf_ref(T,ir)*(0.05d0 + 0.01d0*ir)*exp(-200d0*ir/dble(T))
    enddo
  enddo
  roi0 = [0.05d0, 0.30d0, 0.10d0, 0.08d0, 0.12d0, 0.20d0, 0.01d0, 0.002d0, 0.004d0, 0.05d0, 0.03d0, 0.40d0]
  nfail3 = 0
  do ir = 1, nroutine
    call set_tables(1)
    do k = 1, nsamp
      roi = roi0; w = 0.d0
      call call_routine(ir, roi, Tsamp(k), w)
      w_ref(:,k) = w
    enddo
    nbad = 0
    do m = 1, ntmin
      call set_tables(Tmins(m))
      do k = 1, nsamp
        roi = roi0; w = 0.d0
        call call_routine(ir, roi, Tsamp(k), w)
        if (any(w /= w_ref(:,k)) .or. any(w /= w)) nbad = nbad + 1
      enddo
    enddo
    write(*,'(A,A18,A,I0,A,I0,A,ES10.3)') ' routine ', rname(ir), ': ', nbad, ' of ', ntmin*nsamp, &
      ' (Tmin, T) states differ from the 1 K tables; max |omegadot| at 1 K = ', maxval(abs(w_ref))
    nfail3 = nfail3 + nbad
  enddo

  ! 4) positive control: the assumed-shape accessor of FLINT <= 2223136 reads row T + Tmin - 1
  call set_tables(100)
  nctrl = 0; nctrl_tot = 0; nrest4 = 0
  do ir = 1, nrc_syn
    do k = 1, nsamp
      Tint = [int(Tctrl(k)), int(Tctrl(k)) + 1]; td = Tctrl(k) - int(Tctrl(k))
      nctrl_tot = nctrl_tot + 1
      if (old_comp_ch_tabT(ir, kf_tab, Tint, td) /= f_kf(ir, Tint, td)) nctrl = nctrl + 1
      if (comp_ch_tabT(ir, kf_tab, Tint, td) /= f_kf(ir, Tint, td)) nrest4 = nrest4 + 1
      if (comp_ch_tabT(ir, kb_tab, Tint, td) /= f_kb(ir, Tint, td)) nrest4 = nrest4 + 1
    enddo
  enddo
  write(*,'(A,I0,A,I0,A)') ' positive control (tables from 100 K): the old assumed-shape accessor differs from f_kf in ', &
    nctrl, ' of ', nctrl_tot, ' samples (expected: all)'
  write(*,'(A,I0,A,I0,A)') ' comp_ch_tabT (tables from 100 K): differs from f_kf|f_kb in ', nrest4, ' of ', 2*nctrl_tot, &
    ' samples (expected 0)'
  call free_chemistry_data()
  if (nfail1 + nfail2 + nfail3 > 0 .or. nctrl /= nctrl_tot .or. nctrl1 /= 0 .or. nctrl2 /= nctrl2_tot &
      .or. nrest1 + nrest2 + nrest4 /= 0) then
    write(*,'(A)') ' Verdict -> fail (row T of a rate table must be the rate at T kelvin)'
    stop 1
  endif
  write(*,'(A)') ' Verdict -> pass'

contains

  subroutine set_tables(Tlo)
    integer, intent(in) :: Tlo
    if (allocated(kf_tab)) deallocate(kf_tab)
    if (allocated(kb_tab)) deallocate(kb_tab)
    allocate(kf_tab(Tlo:Tmax_syn, nrc_syn), kb_tab(Tlo:Tmax_syn, nrc_syn))
    kf_tab = kf_ref(Tlo:Tmax_syn, :); kb_tab = kb_ref(Tlo:Tmax_syn, :)
    nrc_arrh = nrc_syn
    T_tab_min = Tlo; T_tab_max = Tmax_syn      ! as read_chemistry sets them
  end subroutine set_tables

  subroutine call_routine(ir, roi, temp, w)
    integer, intent(in) :: ir
    real(8), intent(inout) :: roi(nsyn), w(nsyn)
    real(8), intent(in) :: temp
    select case (ir)
    case (1);  call WD(roi, temp, w)
    case (2);  call Andersen(roi, temp, w)
    case (3);  call OSK(roi, temp, w)
    case (4);  call JLR(roi, temp, w)
    case (5);  call Frassoldati(roi, temp, w)
    case (6);  call CKJLR10sp(roi, temp, w)
    case (7);  call singh(roi, temp, w)
    case (8);  call Singh_WC32(roi, temp, w)
    case (9);  call singhC3H6(roi, temp, w)
    case (10); call Coronetti(roi, temp, w)
    case (11); call Nassini_4(roi, temp, w)
    case (12); call Frolov_nopressure(roi, temp, w)
    case (13); call Frolov(roi, temp, w)
    end select
  end subroutine call_routine

  ! The accessor of FLINT <= 2223136 (Lib_Chemistry_data.f90:36-48): an assumed-shape dummy is
  ! renumbered from 1, so row Tint(1) is the rate at Tint(1) + Tmin - 1 kelvin.
  pure function old_comp_ch_tabT(ireact, tab, Tint, Tdiff) result(res)
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: tab(:,:), Tdiff
    real(8) :: res, a, b
    a = tab(Tint(1), ireact); b = tab(Tint(2), ireact)
    res = a + (b-a)*Tdiff
  end function old_comp_ch_tabT
end program test
