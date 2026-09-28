! test-falloff: the Troe/Lindemann rates where the tables vanish and the k_c convention of the
! falloff tables. In-memory tables: no fixture, no Cantera. Exit code 1
! on failure; compiles against FLINT <= 2223136 too (there it fails: NaN/Infinity in RELEASE, SIGFPE
! under the -ffpe-trap/-fpe0 DEBUG flags).
!  - Pr = 0 (no third body, [M] = 0) and k_inf = 0 (tabulated rate underflowed at low T) give a
!       zero, finite rate instead of NaN (log10(0), 0/0).
!  - k_c <= 0 in the table marks an irreversible reaction: rate(2) = 0; k_c > 0 keeps
!       rate(2) = rate(1)/k_c bit for bit.
program test
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use FLINT_Lib_Chemistry_data
  use FLINT_Lib_Chemistry_falloff, only: f_k_troe, f_k_lindemann, f_F
  implicit none
  integer, parameter :: Tlo = 1, Thi = 3000
  integer :: T, nfail, Tint(2)
  real(8) :: rate(2), Tdiff, kc

  nfail = 0
  nrc_troe = 2; nrc_lindemann = 2
  allocate(kinf_troe_tab(Tlo:Thi,2), k0_troe_tab(Tlo:Thi,2), kc_troe_tab(Tlo:Thi,2), Fcent_tab(Tlo:Thi,2))
  allocate(kinf_lind_tab(Tlo:Thi,2), k0_lind_tab(Tlo:Thi,2), kc_lind_tab(Tlo:Thi,2))
  do T = Tlo, Thi
    kinf_troe_tab(T,:) = 1d6*exp(-3000d0/dble(T))     ! underflows to 0 at T = 1 (exp(-3000))
    k0_troe_tab(T,:)   = 1d4*exp(-100d0/dble(T))      ! > 0 everywhere
    kc_troe_tab(T,1)   = 1d-3*exp(500d0/dble(T))      ! reaction 1: reversible
    kc_troe_tab(T,2)   = 0d0                          ! reaction 2: irreversible (k_c <= 0)
    Fcent_tab(T,:)     = 0.6d0
  enddo
  kinf_lind_tab = kinf_troe_tab; k0_lind_tab = k0_troe_tab; kc_lind_tab = kc_troe_tab

  ! reference state: T = 1500.25 K, [M] = 0.04 kmol/m3
  Tint = [1500, 1501]; Tdiff = 0.25d0
  rate = f_k_troe(1, Tint, Tdiff, 0.04d0); kc = f_kc_troe(1, Tint, Tdiff)
  call verdict('Troe reversible: rate(1) > 0 and finite', rate(1) > 0d0 .and. ieee_is_finite(rate(1)))
  call verdict('Troe reversible: rate(2) == rate(1)/k_c (unchanged path)', rate(2) == rate(1)/kc)
  rate = f_k_troe(2, Tint, Tdiff, 0.04d0)
  call verdict('Troe k_c = 0: rate(1) > 0 and rate(2) == 0', rate(1) > 0d0 .and. rate(2) == 0d0)
  rate = f_k_troe(1, Tint, Tdiff, 0d0)
  call verdict('Troe [M] = 0 (Pr = 0): rate == 0 and finite', all(rate == 0d0) .and. all(ieee_is_finite(rate)))
  ! the guard itself (the rate 0 * F is folded to 0 by -ffast-math whatever F): F(Pr = 0) is exactly 1
  call verdict('Troe Pr = 0: the broadening factor is exactly 1 (no log10(0))', f_F(0d0, 0.6d0) == 1d0)
  call verdict('Troe Pr = 1: the broadening factor is finite and in (0, 1]', &
    ieee_is_finite(f_F(1d0, 0.6d0)) .and. f_F(1d0, 0.6d0) > 0d0 .and. f_F(1d0, 0.6d0) <= 1d0)
  call verdict('table: k_inf underflowed to 0 at 1 K, k_0 > 0', kinf_troe_tab(1,1) == 0d0 .and. k0_troe_tab(1,1) > 0d0)
  Tint = [1, 2]; Tdiff = 0d0
  rate = f_k_troe(1, Tint, Tdiff, 0.04d0)
  call verdict('Troe k_inf = 0: rate == 0 and finite', all(rate == 0d0) .and. all(ieee_is_finite(rate)))
  ! F_cent <= 0 (AramcoMech 2.0 C2H4 + H (+M) above 4871 K): Cantera's convention log10(max(F_cent, 1e-300))
  Fcent_tab(1500:1501,2) = -2.8365d-3
  Tint = [1500, 1501]; Tdiff = 0.25d0
  rate = f_k_troe(2, Tint, Tdiff, 0.04d0)
  call verdict('Troe F_cent < 0: rate finite and equal to the Cantera convention', all(ieee_is_finite(rate)) .and. &
    abs(rate(1) - troe_ref(2, 0.04d0, 1d-300)) <= 1d-12*troe_ref(2, 0.04d0, 1d-300) .and. rate(1) > 0d0)
  Fcent_tab(1500:1501,2) = 0d0
  rate = f_k_troe(2, Tint, Tdiff, 0.04d0)
  call verdict('Troe F_cent = 0: rate finite (Cantera convention)', all(ieee_is_finite(rate)))
  Fcent_tab(1500:1501,2) = 0.6d0
  rate = f_k_troe(2, Tint, Tdiff, 0.04d0)
  call verdict('Troe F_cent = 0.6: rate equal to the reference formula', abs(rate(1) - troe_ref(2, 0.04d0, 0.6d0)) <= &
    1d-12*troe_ref(2, 0.04d0, 0.6d0))

  Tint = [1500, 1501]; Tdiff = 0.25d0
  rate = f_k_lindemann(1, Tint, Tdiff, 0.04d0); kc = f_kc_lind(1, Tint, Tdiff)
  call verdict('Lindemann reversible: rate(1) > 0 and finite', rate(1) > 0d0 .and. ieee_is_finite(rate(1)))
  call verdict('Lindemann reversible: rate(2) == rate(1)/k_c (unchanged path)', rate(2) == rate(1)/kc)
  rate = f_k_lindemann(2, Tint, Tdiff, 0.04d0)
  call verdict('Lindemann k_c = 0: rate(1) > 0 and rate(2) == 0', rate(1) > 0d0 .and. rate(2) == 0d0)
  rate = f_k_lindemann(1, Tint, Tdiff, 0d0)
  call verdict('Lindemann [M] = 0: rate == 0 and finite', all(rate == 0d0) .and. all(ieee_is_finite(rate)))
  Tint = [1, 2]; Tdiff = 0d0
  rate = f_k_lindemann(1, Tint, Tdiff, 0.04d0)
  call verdict('Lindemann k_inf = 0: rate == 0 and finite', all(rate == 0d0) .and. all(ieee_is_finite(rate)))

  call free_chemistry_data()
  if (nfail > 0) then
    write(*,'(A,I0,A)') ' Verdict -> fail (', nfail, ' checks)'
    stop 1
  endif
  write(*,'(A)') ' Verdict -> pass'
contains
  ! Troe rate at Tint/Tdiff of the tables (interpolated as FLINT does), F_cent given
  function troe_ref(ir, cm, fc) result(k)
    integer, intent(in) :: ir
    real(8), intent(in) :: cm, fc
    real(8) :: k, ki, k0, pr, lf, c, n, f1
    ki = kinf_troe_tab(Tint(1),ir) + (kinf_troe_tab(Tint(2),ir) - kinf_troe_tab(Tint(1),ir))*Tdiff
    k0 = k0_troe_tab(Tint(1),ir) + (k0_troe_tab(Tint(2),ir) - k0_troe_tab(Tint(1),ir))*Tdiff
    pr = k0*cm/ki; lf = log10(max(fc, 1d-300)); c = -0.4d0 - 0.67d0*lf; n = 0.75d0 - 1.27d0*lf
    f1 = (log10(pr) + c)/(n - 0.14d0*(log10(pr) + c))
    k = ki*(pr/(1d0 + pr))*10d0**(lf/(1d0 + f1*f1))
  end function troe_ref
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
