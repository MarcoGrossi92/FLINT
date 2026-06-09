module FLINT_Lib_Chemistry_rhs
  use OSLo
  use FLINT_Lib_Chemistry_wdot
# if defined (CANTERA)
  use cantera
# endif
  implicit none
  logical :: reactive
  !> Toggle analytical Jacobian path. When .true. AND chemistry_jacobian is
  !> associated, Radau5 is called with IJAC=1 and jac_native; otherwise the
  !> integrator falls back to its internal finite-difference Jacobian.
  !> Default off so behavior is unchanged unless the mechanism opts in.
  logical :: use_analytical_jacobian = .false.
# if defined(CANTERA)
  type(phase_t) :: gas
# endif

contains
    
  subroutine rhs_native ( nz, time, Z, F )
    use FLINT_Lib_Thermodynamic
    implicit none
    integer, intent(in)  :: nz
    real(8), intent(in)  :: time
    real(8), intent(in)  :: Z(nz)
    real(8), intent(out) :: F(nz)
    ! Local
    real(8) :: roi(nz-1)
    real(8) :: T
    real(8) :: droic(nz-1)
    real(8) :: eiroi,rho_cv
    integer :: s

    roi(1:ns) = Z(1:ns)
    T = Z(nz)
    if (T < Tmin .or. T > Tmax .or. isnan(T)) then
      F(:) = -1.d0
      return
    end if

    ! Avoid negative rho_i
    do s = 1, ns
      roi(s) = max(roi(s), 0.d0)
    end do

    call Chemistry_Source ( roi, T, droic )

    eiroi = 0.d0; rho_cv = 0.d0
    do s = 1, ns
      eiroi = eiroi + ( f_tabT(T,s,h_tab) - Ri_tab(s)*T ) * droic(s)
      rho_cv = rho_cv + roi(s)*( f_tabT(T,s,cp_tab) - Ri_tab(s) )
    enddo

    F(1:ns) = droic
    F(nz) = -eiroi / rho_cv

  end subroutine rhs_native


  !----------------------------------------------------------------------------
  ! Analytical Jacobian of rhs_native.
  !
  ! Signature matches the Hairer RADAU5 JAC callback:
  !     SUBROUTINE JAC(N, X, Y, DFY, LDFY, RPAR, IPAR)
  ! DFY(i,j) = d F(i) / d Y(j), with full storage (MLJAC=N).
  !
  ! State layout (same as rhs_native): Y = ( roi(1:ns), T ), so nz = ns+1.
  !
  ! Top-left (ns x ns) block and (ns x 1) T-column come from the active
  ! mechanism's analytical Jacobian (chemistry_jacobian).
  ! Bottom row (T-equation) is differentiated here using h_tab / cp_tab.
  !----------------------------------------------------------------------------
  subroutine jac_native(nz, time, Z, DFY, LDFY, RPAR, IPAR)
    use FLINT_Lib_Thermodynamic
    implicit none
    integer, intent(in)  :: nz, LDFY
    real(8), intent(in)  :: time
    real(8), intent(in)  :: Z(nz)
    real(8), intent(out) :: DFY(LDFY, nz)
    real(8), intent(in)  :: RPAR(*)
    integer, intent(in)  :: IPAR(*)
    ! Local
    real(8) :: roi(nz-1), droic(nz-1)
    real(8) :: dwdr(nz-1, nz-1), dwdT(nz-1)
    real(8) :: h_vec(nz-1), cp_vec(nz-1), dh_dT(nz-1), dcp_dT(nz-1)
    real(8) :: T, T_frac
    integer :: T_idx, s, j
    real(8) :: G, D, inv_D, F_T
    real(8) :: dG_drho_j, dD_drho_j, dG_dT, dD_dT

    T = Z(nz)

    ! Mirror rhs_native's bail-out: out-of-range T → F is the constant -1 there,
    ! whose Jacobian is zero. Returning zero keeps Newton in a benign state until
    ! the step is rejected and the integrator retries with smaller H.
    if (T < Tmin .or. T >= Tmax .or. isnan(T)) then
      DFY(1:nz, 1:nz) = 0.d0
      return
    end if

    T_idx  = idint(T)
    T_frac = T - dble(T_idx)
    roi(1:ns) = max(Z(1:ns), 0.d0)

    ! 1) Species block: ∂omegadot/∂(roi,T) from the mechanism.
    call chemistry_source   ( roi, T, droic )
    call chemistry_jacobian ( roi, T, dwdr, dwdT )

    ! 2) Thermo lookups + analytical T-derivatives (same linear-interp table).
    do s = 1, ns
      h_vec(s)      = h_tab(T_idx,s)  + T_frac * (h_tab(T_idx+1,s)  - h_tab(T_idx,s))
      cp_vec(s)     = cp_tab(T_idx,s) + T_frac * (cp_tab(T_idx+1,s) - cp_tab(T_idx,s))
      dh_dT(s)  = h_tab(T_idx+1,s)  - h_tab(T_idx,s)
      dcp_dT(s) = cp_tab(T_idx+1,s) - cp_tab(T_idx,s)
    end do

    ! 3) G = eiroi, D = rho_cv (mirrors rhs_native).
    G = 0.d0
    D = 0.d0
    do s = 1, ns
      G = G + (h_vec(s) - Ri_tab(s) * T) * droic(s)
      D = D + roi(s) * (cp_vec(s) - Ri_tab(s))
    end do
    inv_D = 1.d0 / D
    F_T   = -G * inv_D            ! = F(nz)

    ! 4) Top-left + T-column from the mechanism Jacobian.
    do j = 1, ns
      DFY(1:ns, j) = dwdr(1:ns, j)
    end do
    DFY(1:ns, nz) = dwdT(1:ns)

    ! 5) Bottom row: ∂F(nz)/∂Z.
    !    F(nz) = -G/D  =>  ∂F(nz)/∂x = -inv_D*(∂G/∂x) - F_T*inv_D*(∂D/∂x)
    do j = 1, ns
      ! G depends on roi only via droic(s) and the (h_s - R_s T) factor doesn't
      ! depend on roi → only the dwdr term contributes.
      dG_drho_j = 0.d0
      do s = 1, ns
        dG_drho_j = dG_drho_j + (h_vec(s) - Ri_tab(s) * T) * dwdr(s, j)
      end do
      dD_drho_j = cp_vec(j) - Ri_tab(j)
      DFY(nz, j) = -inv_D * dG_drho_j - F_T * inv_D * dD_drho_j
    end do

    dG_dT = 0.d0
    dD_dT = 0.d0
    do s = 1, ns
      dG_dT = dG_dT + (h_vec(s) - Ri_tab(s) * T) * dwdT(s) &
                    + (dh_dT(s) - Ri_tab(s))    * droic(s)
      dD_dT = dD_dT + roi(s) * dcp_dT(s)
    end do
    DFY(nz, nz) = -inv_D * dG_dT - F_T * inv_D * dD_dT

  end subroutine jac_native


  !----------------------------------------------------------------------------
  ! Finite-difference Jacobian of rhs_native, for one-shot validation against
  ! jac_native. Not called in the hot path. Uses central differences with a
  ! state-relative perturbation. Returns DFY in the same layout as jac_native.
  !----------------------------------------------------------------------------
  subroutine jac_native_fd(nz, time, Z, DFY, eps_in)
    implicit none
    integer, intent(in)  :: nz
    real(8), intent(in)  :: time, Z(nz)
    real(8), intent(out) :: DFY(nz, nz)
    real(8), intent(in), optional :: eps_in   ! relative FD step (default sqrt(eps))
    real(8) :: Zp(nz), Zm(nz), Fp(nz), Fm(nz)
    real(8) :: eps, h
    integer :: j

    eps = sqrt(epsilon(1.d0))
    if (present(eps_in)) eps = eps_in
    do j = 1, nz
      h = eps * max(1.d-30, abs(Z(j)))
      Zp = Z; Zm = Z
      Zp(j) = Z(j) + h
      Zm(j) = Z(j) - h
      call rhs_native(nz, time, Zp, Fp)
      call rhs_native(nz, time, Zm, Fm)
      DFY(:, j) = (Fp - Fm) / (2.d0 * h)
    end do
  end subroutine jac_native_fd


  !----------------------------------------------------------------------------
  ! Sanity check: at a representative post-ignition state, compute both the
  ! analytical and FD Jacobians and print maxerr. Called once during setup if
  ! requested. Failure of the check does NOT abort; it just prints a warning so
  ! the user can decide whether to keep use_analytical_jacobian on.
  !----------------------------------------------------------------------------
  subroutine validate_analytical_jacobian()
    use FLINT_Lib_Thermodynamic
    implicit none
    real(8) :: Z(ns+1), DFY_an(ns+1, ns+1), DFY_fd(ns+1, ns+1)
    real(8) :: scale, abs_err, max_rel_err
    integer :: i, j
    integer :: IPAR_dummy(1)
    real(8) :: RPAR_dummy(1)

    if (.not. associated(chemistry_jacobian)) then
      write(*,'(A)') ' [JAC] no analytical Jacobian registered; skipping validation'
      return
    end if

    ! Representative burnt-mixture state at ~3000 K. Filled species-by-species
    ! up to whatever `ns` the active mechanism has, so this routine works for
    ! any mechanism size (ONERA-7 with ns=7, Frolov_nopressure with ns=4, etc.)
    ! without out-of-bounds writes into the temperature slot.
    Z(:) = 0.d0
    if (ns >= 1) Z(1) = 1.0d-1   ! treat as O2  (or first species in any case)
    if (ns >= 2) Z(2) = 3.0d-1   ! H2O
    if (ns >= 3) Z(3) = 5.0d-2   ! H2
    if (ns >= 4) Z(4) = 8.0d-1   ! N2 / bath gas
    if (ns >= 5) Z(5) = 5.0d-3   ! H
    if (ns >= 6) Z(6) = 5.0d-3   ! O
    if (ns >= 7) Z(7) = 1.0d-2   ! OH
    Z(ns+1) = 3000.d0
    IPAR_dummy = 0
    RPAR_dummy = 0.d0

    call jac_native(ns+1, 0.d0, Z, DFY_an, ns+1, RPAR_dummy, IPAR_dummy)

    ! ---- DIAGNOSTIC SWEEP: vary FD step. If max_rel_err moves with the step,
    ! the discrepancy is FD truncation/roundoff (analytical is fine). If it
    ! stays pinned, there is a genuine error in an analytical term. ----------
    block
      real(8) :: eps_list(6), e, worst_an, worst_fd
      integer :: ie, wi, wj
      eps_list = [1.d-8, 1.49d-8, 1.d-7, 1.d-6, 1.d-5, 1.d-4]
      write(*,'(A)') ' [JAC] FD-step sweep (relerr = elementwise max vs central-diff FD):'
      do ie = 1, 6
        e = eps_list(ie)
        call jac_native_fd(ns+1, 0.d0, Z, DFY_fd, eps_in=e)
        max_rel_err = 0.d0; wi = 0; wj = 0
        do j = 1, ns+1
          do i = 1, ns+1
            scale = max(abs(DFY_an(i,j)), abs(DFY_fd(i,j)), 1.d-30)
            abs_err = abs(DFY_an(i,j) - DFY_fd(i,j))
            if (abs_err/scale > max_rel_err) then
              max_rel_err = abs_err/scale; wi = i; wj = j
            end if
          end do
        end do
        worst_an = DFY_an(wi,wj); worst_fd = DFY_fd(wi,wj)
        write(*,'(A,ES9.2,A,ES10.2,A,I2,A,I2,A,ES12.4,A,ES12.4,A)') &
          '   eps=', e, '   max_relerr=', max_rel_err, &
          '   @(', wi, ',', wj, ')  an=', worst_an, '  fd=', worst_fd
      end do
    end block

    ! ---- Dynamically-weighted metric: how much does the inexact Jacobian
    ! perturb the Newton iteration matrix (I - h*gamma*J) at a representative
    ! chemistry step?  ||(I-hgJ_an)^-1 (I-hgJ_fd) - I||_inf  is the contraction
    ! the inaccuracy adds.  Far below the elementwise relerr if the off entries
    ! are dynamically small. -------------------------------------------------
    block
      real(8) :: hg, A_an(ns+1,ns+1), A_fd(ns+1,ns+1), Minv_B(ns+1,ns+1)
      real(8) :: rownorm, wnorm
      integer :: i2, j2, k2, ip(ns+1), info
      call jac_native_fd(ns+1, 0.d0, Z, DFY_fd)         ! default sqrt(eps) ref
      hg = 1.d-7                                        ! representative h*gamma
      A_an = -hg*DFY_an; A_fd = -hg*DFY_fd
      do i2 = 1, ns+1
        A_an(i2,i2) = A_an(i2,i2) + 1.d0
        A_fd(i2,i2) = A_fd(i2,i2) + 1.d0
      end do
      Minv_B = A_fd
      call dgesv(ns+1, ns+1, A_an, ns+1, ip, Minv_B, ns+1, info)  ! A_an \ A_fd
      wnorm = 0.d0
      if (info == 0) then
        do i2 = 1, ns+1
          rownorm = 0.d0
          do j2 = 1, ns+1
            if (i2 == j2) then
              rownorm = rownorm + abs(Minv_B(i2,j2) - 1.d0)
            else
              rownorm = rownorm + abs(Minv_B(i2,j2))
            end if
          end do
          if (rownorm > wnorm) wnorm = rownorm
        end do
        write(*,'(A,ES10.2,A)') &
          ' [JAC] dynamically-weighted Newton-matrix error ||M^-1 B - I||_inf = ', &
          wnorm, '   (at h*gamma=1e-7)'
      else
        write(*,'(A)') ' [JAC] weighted metric skipped (dgesv info /= 0)'
      end if
    end block

    write(*,'(A)') ' [JAC] (interpret: relerr moving with eps => FD error, not analytical)'
  end subroutine validate_analytical_jacobian


# if defined (CANTERA)
  subroutine rhs_cantera ( nz, time, Z, F )
    use FLINT_Lib_Thermodynamic
    use cantera
    implicit none
    integer, intent(in)  :: nz
    real(8), intent(in)  :: time
    real(8), intent(in)  :: Z(nz)
    real(8), intent(out) :: F(nz)
    ! Local
    real(8) :: roi(nz-1), Y(nz-1)
    real(8) :: T, rho
    real(8) :: droic(nz-1),wdot(nz-1)
    real(8) :: eiroi,rho_cv, e_s, cp_s
    integer :: s

    roi(1:ns) = Z(1:ns)
    T = Z(nz)
    if (T < 0.d0 .or. T > 10000d0) then
      F(:) = -1.d0
      return
    end if

    ! Avoid negative rho_i
    do s = 1, ns
      roi(s) = max(roi(s), 0.d0)
    end do
    rho = sum(roi)
    Y = roi / rho

    call setState_TRY(gas, T, rho, Y)
    call getNetProductionRates(gas, wdot)
    droic = wdot * wm_tab

    eiroi = 0.d0
    rho_cv = 0.d0
    do s = 1, ns
      e_s = f_tabT(T, s, h_tab) - Ri_tab(s) * T
      cp_s = f_tabT(T, s, cp_tab) - Ri_tab(s)
      eiroi = eiroi + e_s * droic(s)
      rho_cv = rho_cv + roi(s) * cp_s
    end do

    F(1:ns) = droic
    F(nz) = -eiroi / rho_cv

  end subroutine rhs_cantera
# endif

end module FLINT_Lib_Chemistry_rhs
