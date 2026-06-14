!------------------------------------------------------------------------------
! Device-callable copies of the general-ns thermo/transport functions for the
! GPU Phase-2 thermo-cache kernel. The originals (Lib_ThermoTransport) use
! runtime-sized automatics (Xi(ns), fi(ns,ns), cpi(ns)) which would malloc on
! the device per call AND cannot be safely pinned to NSMAX in-place (a CPU run
! with a larger mechanism would overflow). These device copies are fixed-size
! NSMAX and `acc routine seq`; they run ONLY in the device build, only for the
! device-allowed mechanisms (mech_dev guard: ONERA-7 ns=7, Frolov ns=4 -> ns
! <= NSMAX). The CPU functions are left 100% untouched (bit-identical).
!
! Mirrors exactly: co_k_mi_lam_Wilke_expr, co_fiij, f_cp_expr. The leaf helpers
! f_Rtot / f_molecularWeight / f_tabT_expr carry `acc routine seq` and are
! reused from the originals (no general-ns automatic to convert).
!------------------------------------------------------------------------------
module FLINT_Lib_ThermoTransport_dev
  use FLINT_Lib_Thermodynamic
  use FLINT_Lib_Chemistry_data, only: NSMAX
  implicit none
  private
  public :: co_rotot_Rtot_dev, co_k_mi_lam_Wilke_dev, f_cp_dev

contains

  !> rho = sum(roi), R = f_Rtot(roi). (co_rotot_Rtot has no array automatic; a
  !> thin device wrapper for symmetry / a single call site.)
  pure subroutine co_rotot_Rtot_dev ( rhoi, rho, R )
    !$acc routine seq
    implicit none
    real(8), intent(in)  :: rhoi(ns)
    real(8), intent(out) :: rho, R
    rho = sum(rhoi(1:ns))
    R   = f_Rtot(rhoi)
  end subroutine co_rotot_Rtot_dev


  !> Wilke phi_ij matrix. fi is the OUTPUT; its leading dimension is pinned to
  !> NSMAX so sequence association with the caller's fi(NSMAX,NSMAX) is exact
  !> (same leading-extent hazard as the chemistry dwdr). Only (1:ns,1:ns) set.
  pure subroutine co_fiij_dev ( fi, mi )
    !$acc routine seq
    implicit none
    real(8), intent(in)  :: mi(ns)
    real(8), intent(out) :: fi(NSMAX, ns)
    real(8) :: sqrt_mi_num(NSMAX), inv_sqrt_mi_den(NSMAX), ratio, one_plus
    integer :: i, j

    do i = 1, ns
      sqrt_mi_num(i)     = dsqrt(mi(i))
      inv_sqrt_mi_den(i) = 1.d0 / dsqrt(mi(i) + 1d-20)
    enddo
    do j = 1, ns
      do i = 1, ns
        ratio    = sqrt_mi_num(i) * inv_sqrt_mi_den(j) * Mi_Mj_pow_m025(i,j)
        one_plus = 1.d0 + ratio
        fi(i,j)  = one_plus * one_plus * inv_sqrt8_1p(i,j)
      enddo
    enddo
  end subroutine co_fiij_dev


  !> Laminar viscosity + conductivity via Wilke mixing. Fixed-size NSMAX
  !> automatics (1:ns used). Bit-faithful to co_k_mi_lam_Wilke_expr.
  pure subroutine co_k_mi_lam_Wilke_dev ( rhoi, rho, Tint, Tdiff, milam, klam )
    !$acc routine seq
    implicit none
    integer, intent(in)  :: Tint(2)
    real(8), intent(in)  :: rhoi(ns), rho, Tdiff
    real(8), intent(out) :: milam, klam
    real(8) :: lam_den(NSMAX), Xi(NSMAX)
    real(8) :: fi(NSMAX, NSMAX), klam_i(NSMAX), milam_i(NSMAX)
    real(8) :: Wmtot, inv_lam_den
    integer :: s, i, j

    Wmtot = f_molecularWeight(rhoi)
    do s = 1, ns
      Xi(s)      = rhoi(s)*Wmtot/(rho*Wm_tab(s))
      milam_i(s) = f_tabT_expr(s, mi_tab, Tint, Tdiff)
      klam_i(s)  = f_tabT_expr(s, k_tab,  Tint, Tdiff)
    enddo
    call co_fiij_dev(fi, milam_i)
    do i = 1, ns
      lam_den(i) = 0.d0
      do j = 1, ns
        lam_den(i) = lam_den(i) + Xi(j)*fi(i,j)
      enddo
    enddo
    milam = 0.d0
    klam  = 0.d0
    do s = 1, ns
      inv_lam_den = 1.d0/lam_den(s)
      milam = milam + Xi(s)*milam_i(s)*inv_lam_den
      klam  = klam  + Xi(s)*klam_i(s) *inv_lam_den
    enddo
  end subroutine co_k_mi_lam_Wilke_dev


  !> Mixture cp. Fixed-size cpi(NSMAX). Bit-faithful to f_cp_expr.
  pure function f_cp_dev ( rhoi, Tint, Tdiff, rho ) result(res)
    !$acc routine seq
    implicit none
    real(8), intent(in) :: rhoi(ns), Tdiff, rho
    integer, intent(in) :: Tint(2)
    real(8) :: res
    real(8) :: cpi(NSMAX)
    integer :: s
    res = 0.d0
    do s = 1, ns
      cpi(s) = f_tabT_expr(s, cp_tab, Tint, Tdiff)
      res = res + rhoi(s)/rho*cpi(s)
    enddo
  end function f_cp_dev

end module FLINT_Lib_ThermoTransport_dev
