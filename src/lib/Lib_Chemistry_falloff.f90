! Source: https://cantera.org/3.1/reference/kinetics/rate-constants.html#sec-falloff-rate
module FLINT_Lib_Chemistry_falloff
  use FLINT_Lib_Chemistry_data
  implicit none
  private
  public :: f_k_troe
  public :: f_k_lindemann
  public :: f_reverse
  public :: f_F

contains

  ! Computes the Lindemann rates for a reaction given temperature, third-body concentration, and reaction index
  pure function f_k_lindemann(ireact, Tint, Tdiff, coM) result(rate)
    real(8), intent(in) :: coM
    integer, intent(in) :: ireact
    integer, intent(in) :: Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: rate(2)
    real(8) :: kinf, k0, kc, Fcent, Pr

    ! Compute the Troe rate using the precomputed tables
    kinf = f_kinf_lind(ireact, Tint, Tdiff)
    k0 = f_k0_lind(ireact, Tint, Tdiff)
    kc = f_kc_lind(ireact, Tint, Tdiff)

    ! Vanished high-pressure limit (underflow of the tabulated rate at low T):
    ! k = kinf*Pr/(1+Pr) <= kinf = 0, and Pr = k0*coM/kinf would be 0/0.
    if (kinf <= 0d0) then
      rate = 0d0
      return
    endif

    ! Reduced pressure
    Pr = k0 * coM / kinf

    ! Forward rate
    rate(1) = kinf * ( Pr / (1d0+Pr) )
    ! Backward rate: kc <= 0 marks an irreversible reaction (no reverse step)
    rate(2) = f_reverse(rate(1), kc)

  end function f_k_lindemann

  ! Computes the Troe rates for a reaction given temperature, third-body concentration, and reaction index
  pure function f_k_troe(ireact, Tint, Tdiff, coM) result(rate)
    real(8), intent(in) :: coM
    integer, intent(in) :: ireact
    integer, intent(in) :: Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: rate(2)
    real(8) :: kinf, k0, kc, Fcent, Pr

    ! Compute the Troe rate using the precomputed tables
    kinf = f_kinf_troe(ireact, Tint, Tdiff)
    k0 = f_k0_troe(ireact, Tint, Tdiff)
    kc = f_kc_troe(ireact, Tint, Tdiff)
    Fcent = f_Fcent(ireact, Tint, Tdiff)

    ! Vanished high-pressure limit (underflow of the tabulated rate at low T):
    ! k = kinf*Pr/(1+Pr)*F <= kinf = 0, and Pr = k0*coM/kinf would be 0/0.
    if (kinf <= 0d0) then
      rate = 0d0
      return
    endif

    ! Reduced pressure
    Pr = k0 * coM / kinf

    ! Forward rate
    rate(1) = kinf * ( Pr / (1d0+Pr) ) * f_F(Pr, Fcent)
    ! Backward rate: kc <= 0 marks an irreversible reaction (no reverse step)
    rate(2) = f_reverse(rate(1), kc)

  end function f_k_troe

  ! Reverse rate coefficient from the forward one and the equilibrium constant
  ! (concentration units) tabulated in the k_c column of chemistry-Troe.dat /
  ! chemistry-Lindemann.dat. Contract with the table writer (format of these files): a value
  ! kc <= 0 marks an irreversible reaction and gives no reverse step; a table
  ! written before this contract carries the equilibrium constant also for
  ! irreversible falloff reactions and keeps its (spurious) reverse rate.
  pure function f_reverse(kfwd, kc) result(krev)
    real(8), intent(in) :: kfwd, kc
    real(8) :: krev
    if (kc > 0d0) then
      krev = kfwd / kc
    else
      krev = 0d0
    endif
  end function f_reverse

  ! Computes the Troe falloff correction factor F
  pure function f_F(Pr, Fcent) result(F)
    implicit none

    real(8), intent(in)  :: Pr      ! Reduced pressure
    real(8), intent(in)  :: Fcent   ! Known center factor
    real(8)              :: F       ! Final Troe correction

    real(8) :: logFcent, logPr, c, n, f1

    ! Low-pressure limit: Pr = 0 when no third body is present (or k0 underflowed).
    ! The rate kinf*Pr/(1+Pr)*F is then 0 whatever F; the limit of the Troe
    ! factor for Pr -> 0 is Fcent**(1/(1+1/0.14**2)) = Fcent**0.0192 (0.96 for
    ! Fcent = 0.1), so F = 1 is taken there instead of log10(0) = -Infinity
    ! (NaN in F, SIGFPE under -ffpe-trap/-fpe0); the rate stays 0 either way.
    if (Pr <= 0d0) then
      F = 1d0
      return
    endif
    ! The Troe centre factor (1-a)exp(-T/T3)+a exp(-T/T1)+exp(-T2/T) is not positive
    ! for every published parameter set (a < 0 or a > 1, T1 or T3 < 0): e.g.
    ! C2H4 + H (+M) <=> C2H5 (+M) of AramcoMech 2.0/3.0 and FFCM-1 has F_cent <= 0
    ! above 4871 K. Cantera takes log10(max(F_cent, 1e-300)) there (the rate becomes
    ! negligible, it does not stop); the same convention is used here, so that the
    ! tables give Cantera's rates. Bit-identical to log10(Fcent) for Fcent > 1e-300.
    logFcent = log10(max(Fcent, 1d-300))
    logPr = log10(Pr)

    c = -0.4d0 - 0.67d0 * logFcent
    n = 0.75d0 - 1.27d0 * logFcent

    f1 = (logPr + c) / (n - 0.14d0* (logPr + c))

    F = 10d0**(logFcent / (1d0 + f1*f1))

  end function f_F


end module FLINT_Lib_Chemistry_falloff