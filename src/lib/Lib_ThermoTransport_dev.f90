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
  public :: co_rotot_Rtot_dev, co_k_mi_lam_Wilke_dev, f_cp_dev, f_ss_dev, f_ss_cp_dev, H0_dev

! MOSE_TT_CT (compile-time dims extension; mirrors MOSE_CHEM_CT in Lib_Radau5_dev):
! when the build pins the mechanism size (MOSE_TT_NS = MOSE_NSC from CMake), the
! species loops below get LITERAL trip counts (full unroll) and the NSMAX-padded
! locals shrink to the exact size (register promotion; Wilke's fi drops 8x8 -> ns^2).
! The species-loop leaves (f_Rtot / f_molecularWeight / f_tabT_expr) are then routed
! to same-module _m copies so nvfortran inlines them into the _dev routines --
! cross-module acc-routine calls are never inlined (the MOSE_CHEM_INLINE finding).
! Executable statements untouched -> bit-exact class. Runtime-guarded in
! flint_acc_upload_thermo (ns == MOSE_TT_NS, error stop on mismatch).
# if defined (MOSE_TT_NS)
#   define NSL MOSE_TT_NS
#   define NSPAD MOSE_TT_NS
#   define TT_F_RTOT f_Rtot_m
#   define TT_F_MOLW f_molecularWeight_m
#   define TT_F_TABT f_tabT_expr_m
# else
#   define NSL ns
#   define NSPAD NSMAX
#   define TT_F_RTOT f_Rtot
#   define TT_F_MOLW f_molecularWeight
#   define TT_F_TABT f_tabT_expr
# endif

contains

  !> Total specific enthalpy, device version. Mirror of H0/H (Lib_ThermoTransport)
  !> built from the already device-safe leaves f_Rtot + f_tabT_expr (h_tab is
  !> `declare create`'d). Kept here in the _dev module so the SHARED CPU thermo
  !> functions stay 100% untouched (annotating them with acc routine seq broke
  !> the FLINT chemistry modules' device globals -> mech_dev declare-create).
  pure function H0_dev(p, rhoi, u) result(entalpy)
    !$acc routine seq
    implicit none
    real(8), intent(in) :: p, rhoi(NSL), u
    real(8) :: entalpy, rho, Rgas, T, Tdiff, h_
    integer :: s, T_i, Tint(2)
    rho  = sum(rhoi(1:NSL))
    Rgas = TT_F_RTOT(rhoi)
    T    = p/(Rgas*rho)
    T_i  = idint(T); Tdiff = T - T_i
    Tint(1) = T_i; Tint(2) = T_i + 1
    entalpy = 0d0
    do s = 1, NSL
      h_ = TT_F_TABT(s, h_tab, Tint, Tdiff)
      entalpy = entalpy + h_*rhoi(s)/rho
    end do
    entalpy = entalpy + 0.5d0*u*u
  end function H0_dev

  !> Frozen speed of sound, device version. Mirror of f_ss (Lib_ThermoTransport)
  !> but uses the fixed-size f_cp_dev (the original f_cp has a general-ns cpi(ns)
  !> automatic that would malloc on the device).
  pure function f_ss_dev(rhoi, p, rho, Rtot) result(res)
    !$acc routine seq
    implicit none
    real(8), intent(in) :: rhoi(NSL), p, rho, Rtot
    real(8) :: res, T, cp, gam, Tdiff
    integer :: T_i, Tint(2)
    T = p/(Rtot*rho)
    T_i = idint(T); Tdiff = T - T_i
    Tint(1) = T_i; Tint(2) = T_i + 1
    cp  = f_cp_dev(rhoi, Tint, Tdiff, rho)
    gam = cp/(cp - Rtot)
    res = dsqrt(gam*Rtot*T)
  end function f_ss_dev

  !> Frozen speed of sound from a PRECOMPUTED mixture cp (the resident thermo-cache
  !> cp_c). Statement-for-statement f_ss_dev with the f_cp_dev call substituted by
  !> the argument: whenever cp carries the bits f_cp_dev would produce at this state
  !> (cp_c does -- same rhoi/rho/Rtot/T operands, see thermo_cache_visc_W_blk), the
  !> result is bit-identical while skipping the full species gather + cp table loop.
  pure function f_ss_cp_dev(cp, p, rho, Rtot) result(res)
    !$acc routine seq
    implicit none
    real(8), intent(in) :: cp, p, rho, Rtot
    real(8) :: res, T, gam
    T = p/(Rtot*rho)
    gam = cp/(cp - Rtot)
    res = dsqrt(gam*Rtot*T)
  end function f_ss_cp_dev

  !> rho = sum(roi), R = f_Rtot(roi). (co_rotot_Rtot has no array automatic; a
  !> thin device wrapper for symmetry / a single call site.)
  pure subroutine co_rotot_Rtot_dev ( rhoi, rho, R )
    !$acc routine seq
    implicit none
    real(8), intent(in)  :: rhoi(NSL)
    real(8), intent(out) :: rho, R
    rho = sum(rhoi(1:NSL))
    R   = TT_F_RTOT(rhoi)
  end subroutine co_rotot_Rtot_dev


  !> Wilke phi_ij matrix. fi is the OUTPUT; its leading dimension is pinned to
  !> NSMAX so sequence association with the caller's fi(NSMAX,NSMAX) is exact
  !> (same leading-extent hazard as the chemistry dwdr). Only (1:ns,1:ns) set.
  pure subroutine co_fiij_dev ( fi, mi )
    !$acc routine seq
    implicit none
    real(8), intent(in)  :: mi(NSL)
    real(8), intent(out) :: fi(NSPAD, NSL)
    real(8) :: sqrt_mi_num(NSPAD), inv_sqrt_mi_den(NSPAD), ratio, one_plus
    integer :: i, j

    do i = 1, NSL
      sqrt_mi_num(i)     = dsqrt(mi(i))
      inv_sqrt_mi_den(i) = 1.d0 / dsqrt(mi(i) + 1d-20)
    enddo
    do j = 1, NSL
      do i = 1, NSL
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
    real(8), intent(in)  :: rhoi(NSL), rho, Tdiff
    real(8), intent(out) :: milam, klam
    real(8) :: lam_den(NSPAD), Xi(NSPAD)
    real(8) :: fi(NSPAD, NSPAD), klam_i(NSPAD), milam_i(NSPAD)
    real(8) :: Wmtot, inv_lam_den
    integer :: s, i, j

    Wmtot = TT_F_MOLW(rhoi)
    do s = 1, NSL
      Xi(s)      = rhoi(s)*Wmtot/(rho*Wm_tab(s))
      milam_i(s) = TT_F_TABT(s, mi_tab, Tint, Tdiff)
      klam_i(s)  = TT_F_TABT(s, k_tab,  Tint, Tdiff)
    enddo
    call co_fiij_dev(fi, milam_i)
    do i = 1, NSL
      lam_den(i) = 0.d0
      do j = 1, NSL
        lam_den(i) = lam_den(i) + Xi(j)*fi(i,j)
      enddo
    enddo
    milam = 0.d0
    klam  = 0.d0
    do s = 1, NSL
      inv_lam_den = 1.d0/lam_den(s)
      milam = milam + Xi(s)*milam_i(s)*inv_lam_den
      klam  = klam  + Xi(s)*klam_i(s) *inv_lam_den
    enddo
  end subroutine co_k_mi_lam_Wilke_dev


  !> Mixture cp. Fixed-size cpi(NSMAX). Bit-faithful to f_cp_expr.
  pure function f_cp_dev ( rhoi, Tint, Tdiff, rho ) result(res)
    !$acc routine seq
    implicit none
    real(8), intent(in) :: rhoi(NSL), Tdiff, rho
    integer, intent(in) :: Tint(2)
    real(8) :: res
    real(8) :: cpi(NSPAD)
    integer :: s
    res = 0.d0
    do s = 1, NSL
      cpi(s) = TT_F_TABT(s, cp_tab, Tint, Tdiff)
      res = res + rhoi(s)/rho*cpi(s)
    enddo
  end function f_cp_dev


# if defined (MOSE_TT_NS)
  !----------------------------------------------------------------------------
  ! MOSE_TT_CT same-module leaf copies -- KEEP IN SYNC with the originals in
  ! Lib_ThermoTransport (FLINT_Lib_Thermodynamic). Bodies VERBATIM; only the
  ! species extent is the literal NSL so the loops unroll and the calls inline
  ! into the _dev routines above (same translation unit).
  !----------------------------------------------------------------------------
  pure function f_tabT_expr_m(sp,tab,Tint,Tdiff) result(result)
    !$acc routine seq
    implicit none
    real(8), intent(in) :: Tdiff, tab(Tmin:Tmax,NSL)
    integer, intent(in) :: sp, Tint(2)
    real(8) :: result
    real(8) :: Vij,Viij

    Vij=tab(Tint(1),sp)       ! int(T)   <- Tint(1)
    Viij=tab(Tint(2),sp)      ! int(T)+1 <- Tint(2)
    result = Vij+(Viij-Vij)*Tdiff

  endfunction f_tabT_expr_m


  pure function f_molecularWeight_m(rhoi) result(result)
    !$acc routine seq
    implicit none
    real(8), intent(in) :: rhoi(NSL)
    real(8) :: result
    integer :: s
    real(8) :: rho

    rho = sum(rhoi)
    result = 0.d0
    do s = 1, NSL
      result = result+(rhoi(s)/rho/Wm_tab(s))
    enddo
    result = 1.d0/result

  endfunction f_molecularWeight_m


  pure function f_Rtot_m(rhoi) result(result)
    !$acc routine seq
    implicit none
    real(8), intent(in)  :: rhoi(NSL)
    real(8) :: result
    integer :: s
    real(8) :: rho

    rho = sum(rhoi)
    result = 0.d0
    do s = 1, NSL
      result = result+Ri_tab(s)*rhoi(s)/rho
    enddo

  endfunction f_Rtot_m
# endif

end module FLINT_Lib_ThermoTransport_dev
