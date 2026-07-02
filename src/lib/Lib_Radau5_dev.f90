!------------------------------------------------------------------------------
! Device-callable RADAU5 for the per-cell chemistry integration (GPU port,
! Phase 1). Faithful specialization of OSlo's hairerRadau5.f (Hairer & Wanner,
! version 1996-07-09 / 2002-01-18, as modified in OSlo: COMMON blocks already
! converted to arguments) for exactly the configuration MOSE's
! wrap_radau5_jac uses:
!
!   ITOL=1 (vector tolerances), IJAC=1 (analytical Jacobian), MLJAC=N (full),
!   IMAS=0 (explicit ODE), IOUT=0 (no dense output), M1=0, no Hessenberg
!   (IJOB=1), STARTN=.FALSE., PRED=.TRUE. (Gustafsson controller), NIT=7,
!   NMAX=100000, UROUND=1.1d-19 (the value OSlo passes via WORK(1)),
!   SAFE=0.9, THET=0.001, QUOT1=1, QUOT2=1.2, FACL=5, FACR=1/8,
!   HMAX=XEND-X, initial H=1e-6.
!
! FCN/JAC procedure arguments are replaced by DIRECT calls to
! rhs_native/jac_native (procedure pointers cannot run on the GPU); those in
! turn call ONERA_7/ONERA_7_jac directly in the _OPENACC build.
! All WRITE/STOP removed: failures return IDID<0 exactly where the original
! did (IDID=-2 too many steps, -3 step too small, -4 repeatedly singular).
! Work arrays are fixed-size per-thread locals (NSMX), no SAVE, no COMMON.
!
! Expressions, constants, evaluation order and control flow are copied from
! the original so that a host build of this routine reproduces the OSlo
! RADAU5 trajectories bit-for-bit (verified by the single-cell spike).
!------------------------------------------------------------------------------
module FLINT_Lib_Radau5_dev
  use FLINT_Lib_Chemistry_rhs, only: rhs_native, jac_native
  use FLINT_Lib_Chemistry_data, only: NZMAX
# if defined (MOSE_CHEM_INLINE)
  ! MOSE_CHEM_INLINE: same-module device copies of the whole per-cell chemistry
  ! chain (rhs/jac + ONERA-7/Frolov mechanisms + kf/kb lookups) so nvfortran can
  ! inline them into radau5_dev -- cross-module acc-routine calls are never
  ! inlined (-Minfo confirmed), costing ~7-10 ABI calls per cell per Newton
  ! sweep. Bodies are VERBATIM copies (suffix _m); KEEP IN SYNC with
  ! Lib_Chemistry_rhs.f90 / Lib_ChemMech/{ONERA-7,Frolov_nopressure}.f90.
  use FLINT_Lib_Thermodynamic, only: ns, Tmin, Tmax, wm_tab, Ri_tab, cp_tab, h_tab
  use FLINT_Lib_Chemistry_data, only: NSCHEM, kf_tab, kb_tab
  use FLINT_Lib_Chemistry_wdot, only: mech_dev
  use, intrinsic :: ieee_arithmetic
# endif
  implicit none
  private
  public :: radau5_dev, flint_acc_upload_tables, flint_acc_upload_thermo

  !> Max system size for the fixed per-thread work arrays (ns+1). Single source
  !> of truth = NZMAX (Lib_Chemistry_data); ONERA-7 -> 8. Keep TIGHT: the four
  !> NSMX x NSMX matrices live in per-thread (CUDA local) memory, and every
  !> process context on a GPU reserves local backing for its max resident
  !> threads -- at NSMX=16 that is ~1.4 GB *per MPI rank*, which OOMed 40 ranks
  !> sharing one 24 GB A30.
  integer, parameter, public :: NSMX = NZMAX

  !> Solver tolerances for the device cell loop (set once at MOSE setup from
  !> the same RT/AT given to setup_odesolver; OSlo keeps its copies private).
  real(8), public :: rtol_dev(NSMX) = 1d-5
  real(8), public :: atol_dev(NSMX) = 1d-5
  !$acc declare create(rtol_dev, atol_dev)
  public :: set_radau5_dev_tols

  !> Hairer ITOL=1 tolerance transform hoisted out of the per-cell integrator
  !> (MOSE_CHEM_CT): rtolt/atolt = transformed rtol/atol, fnewt precomputed.
  !> Computed ON DEVICE in set_radau5_dev_tols so the pow bits are identical to
  !> the in-kernel transform they replace (host libm pow may differ by ULPs).
  real(8), public :: rtolt_dev(NSMX) = 1d-5
  real(8), public :: atolt_dev(NSMX) = 1d-5
  real(8), public :: fnewt_dev = 0d0
  !$acc declare create(rtolt_dev, atolt_dev)

contains

  subroutine set_radau5_dev_tols(n, rt, at)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: rt(n), at(n)
    integer :: i
    real(8) :: expm, quot, tolst, fw
    rtol_dev(1:n) = rt(1:n)
    atol_dev(1:n) = at(1:n)
#   if defined (_OPENACC)
    !$acc update device(rtol_dev, atol_dev)
    ! One-time device-side Hairer ITOL=1 transform (MOSE_CHEM_CT consumes these;
    ! computed on the GPU so the pow/sqrt bits match the in-kernel transform).
    ! uround literal = the 1.1d-19 radau5_dev uses (OSlo WORK(1)).
    !$acc serial copyout(fw) private(i, expm, quot, tolst)
    expm = 2.0d0/3.0d0
    do i = 1, n
      quot         = atol_dev(i)/rtol_dev(i)
      rtolt_dev(i) = 0.1d0*rtol_dev(i)**expm
      atolt_dev(i) = rtolt_dev(i)*quot
    end do
    tolst = rtolt_dev(1)
    fw = max(10*1.1d-19/tolst, min(0.03d0, tolst**0.5d0))
    !$acc end serial
    fnewt_dev = fw
#   else
    ! Host mirror (kept consistent for host builds; radau5_dev's own transform
    ! remains the authority when MOSE_CHEM_NZ is not defined).
    expm = 2.0d0/3.0d0
    do i = 1, n
      quot         = at(i)/rt(i)
      rtolt_dev(i) = 0.1d0*rt(i)**expm
      atolt_dev(i) = rtolt_dev(i)*quot
    end do
    tolst = rtolt_dev(1)
    fnewt_dev = max(10*1.1d-19/tolst, min(0.03d0, tolst**0.5d0))
#   endif
  end subroutine set_radau5_dev_tols

  !----------------------------------------------------------------------------
  ! Upload ONLY the read-only thermo/transport scalars + tables (ns, Tmin, Tmax,
  ! wm_tab, Ri_tab, cp_tab, h_tab + transport tables). These back the device
  ! convective/diffusive/thermo kernels and are needed by EVERY device run,
  ! including FROZEN (no-mechanism) cases. Kept separate from the rate-table
  ! upload because kf_tab/kb_tab are allocatable and stay UNALLOCATED when no
  ! mechanism is read -> updating them then would fault. Host-only; call once
  ! after read_idealgas_thermo (+ read_idealgas_transport for VISC).
  !----------------------------------------------------------------------------
  subroutine flint_acc_upload_thermo()
    use FLINT_Lib_Thermodynamic, only: ns, Tmin, Tmax, wm_tab, Ri_tab, &
                                       cp_tab, dcpi_tab, h_tab, mi_tab, k_tab, &
                                       Mi_Mj_pow_m025, inv_sqrt8_1p
    implicit none
#   if defined (_OPENACC)
    !$acc update device(ns, Tmin, Tmax)
    !$acc update device(wm_tab, Ri_tab, cp_tab, dcpi_tab, h_tab)
    ! Transport tables for the device thermo-cache kernel (VISC cases only;
    ! allocated by read_idealgas_transport before this routine runs).
    if (allocated(mi_tab)) then
      !$acc update device(mi_tab, k_tab, Mi_Mj_pow_m025, inv_sqrt8_1p)
    end if
#   endif
  end subroutine flint_acc_upload_thermo

  !----------------------------------------------------------------------------
  ! Copy the read-only chemistry/thermo tables to the device. Host-only;
  ! call once after read_idealgas_thermo + read_chemistry. The matching
  ! "!$acc declare create" directives live in the table modules. Reactive path
  ! only (kf_tab/kb_tab are the chemistry rate tables); frozen runs must call
  ! flint_acc_upload_thermo instead.
  !----------------------------------------------------------------------------
  subroutine flint_acc_upload_tables()
    use FLINT_Lib_Chemistry_data, only: kf_tab, kb_tab
    implicit none
#   if defined (_OPENACC)
    call flint_acc_upload_thermo()
    !$acc update device(kf_tab, kb_tab)
#   endif
  end subroutine flint_acc_upload_tables


  !----------------------------------------------------------------------------
  ! Integrate Y over [X, XEND]. RTOL/ATOL are NOT modified (the internal
  ! Hairer transform works on copies). Statistics returned per call.
  !----------------------------------------------------------------------------
! MOSE_CHEM_CT: compile-time ODE size. `n` becomes a PARAMETER (= MOSE_CHEM_NZ =
! nsc+1, one binary per mechanism) so every loop trip count below is a literal
! (unroll + register promotion of the per-thread work arrays, which shrink from
! NSMX to n). Executable statements are UNTOUCHED -> bit-exact by construction.
! The launcher guards n_in == MOSE_CHEM_NZ (error stop on mismatch).
# if defined (MOSE_CHEM_NZ)
#   define RD5_LD n
# else
#   define RD5_LD NSMX
# endif
! MOSE_CHEM_INLINE: route the integrator's rhs/jac calls to the same-module
! copies (inlinable); default = the cross-module originals (byte-identical).
# if defined (MOSE_CHEM_INLINE)
#   define RD5_RHS rhs_m
#   define RD5_JAC jac_m
# else
#   define RD5_RHS rhs_native
#   define RD5_JAC jac_native
# endif
  subroutine radau5_dev(n_in, x_in, y, xend, rtol_in, atol_in, &
# if defined (MOSE_CHEM_NZ)
                        fnewt_in, &
# endif
                        idid, nfcn, njac, nstep, naccpt, nrejct)
    !$acc routine seq
    implicit none
    integer, intent(in)    :: n_in
    real(8), intent(in)    :: x_in, xend
    real(8), intent(inout) :: y(n_in)
    real(8), intent(in)    :: rtol_in(n_in), atol_in(n_in)
# if defined (MOSE_CHEM_NZ)
    real(8), intent(in)    :: fnewt_in
# endif
    integer, intent(out)   :: idid
    integer, intent(out)   :: nfcn, njac, nstep, naccpt, nrejct

# if defined (MOSE_CHEM_NZ)
    integer, parameter :: n = MOSE_CHEM_NZ
# else
    integer :: n
# endif

    ! ---- fixed per-thread work arrays (original: partitions of WORK/IWORK)
    real(8) :: z1(RD5_LD), z2(RD5_LD), z3(RD5_LD), y0(RD5_LD), scal(RD5_LD)
    real(8) :: f1(RD5_LD), f2(RD5_LD), f3(RD5_LD)
    real(8) :: fjac(RD5_LD,RD5_LD), e1(RD5_LD,RD5_LD), e2r(RD5_LD,RD5_LD), e2i(RD5_LD,RD5_LD)
    real(8) :: cont(4*RD5_LD)
    integer :: ip1(RD5_LD), ip2(RD5_LD)
    real(8) :: rtol(RD5_LD), atol(RD5_LD)
    real(8) :: rpar_d(1)
    integer :: ipar_d(1)

    ! ---- driver parameters (RADAU5 defaults for the wrap_radau5_jac call)
    real(8) :: uround, expm, quot, safe, thet, tolst, fnewt
    real(8) :: quot1, quot2, hmax, facl, facr
    integer :: nmax, nit
    logical :: startn, pred

    ! ---- RADCOR locals (names kept from the original)
    real(8) :: x, h
    real(8) :: sq6, c1, c2, c1m1, c2m1, c1mc2, dd1, dd2, dd3
    real(8) :: u1, alph, beta, cno
    real(8) :: t11, t12, t13, t21, t22, t23, t31
    real(8) :: ti11, ti12, ti13, ti21, ti22, ti23, ti31, ti32, ti33
    real(8) :: posneg, hmaxn, hold, hopt, faccon, cfac, xold, hhfac
    real(8) :: fac1, alphn, betan, xph
    real(8) :: c3q, c1q, c2q, ak1, ak2, ak3, z1i, z2i, z3i
    real(8) :: a1, a2, a3, dyno, denom, dynold, thq, thqold, theta
    real(8) :: dyth, qnewt, f1i, f2i, f3i, err, fac, hnew
    real(8) :: facgus, hacc, erracc, qt, ak, acont3
    real(8) :: hee1, hee2, hee3, ysafe, delt
    integer :: i, j, ier, n2, n3, newt, nsing
    logical :: reject, first, caljac, last

    ! =========================================================================
    ! RADAU5 driver part (specialized: all input checks statically satisfied)
    ! =========================================================================
# if !defined (MOSE_CHEM_NZ)
    n = n_in
# endif
    nfcn=0; njac=0; nstep=0; naccpt=0; nrejct=0
    uround = 1.1d-19            ! OSlo: WORK(1)=1.1d-19
    nmax   = 100000
    nit    = 7
    startn = .false.
    pred   = .true.
    safe   = 0.9d0
    thet   = 0.001d0
    quot1  = 1.d0
    quot2  = 1.2d0
    facl   = 5.d0
    facr   = 1.d0/8.0d0
    x      = x_in
    hmax   = xend - x

    ! -------- check and change the tolerances (ITOL=1 branch)
# if defined (MOSE_CHEM_NZ)
    ! MOSE_CHEM_CT: tolerances arrive PRE-TRANSFORMED (rtolt_dev/atolt_dev, computed
    ! once on device by set_radau5_dev_tols) + fnewt precomputed -- the per-cell
    ! pow/div transform is skipped. Same bits: the init kernel runs the identical
    ! expressions with the identical device pow.
    do i = 1, n
      rtol(i) = rtol_in(i)
      atol(i) = atol_in(i)
    end do
    fnewt = fnewt_in
# else
    expm = 2.0d0/3.0d0
    do i = 1, n
      quot    = atol_in(i)/rtol_in(i)
      rtol(i) = 0.1d0*rtol_in(i)**expm
      atol(i) = rtol(i)*quot
    end do
    ! --- fnewt (computed from the TRANSFORMED rtol(1), as in the original)
    tolst = rtol(1)
    fnewt = max(10*uround/tolst, min(0.03d0, tolst**0.5d0))
# endif

    ! =========================================================================
    ! RADCOR core (IJOB=1, IMPLCT=F, BANDED=F, M1=0, M2=N, NM1=N, IOUT=0)
    ! =========================================================================
    n2 = 2*n
    n3 = 3*n
    sq6  = dsqrt(6.d0)
    c1   = (4.d0-sq6)/10.d0
    c2   = (4.d0+sq6)/10.d0
    c1m1 = c1-1.d0
    c2m1 = c2-1.d0
    c1mc2= c1-c2
    dd1  = -(13.d0+7.d0*sq6)/3.d0
    dd2  = (-13.d0+7.d0*sq6)/3.d0
    dd3  = -1.d0/3.d0
    u1   = (6.d0+81.d0**(1.d0/3.d0)-9.d0**(1.d0/3.d0))/30.d0
    alph = (12.d0-81.d0**(1.d0/3.d0)+9.d0**(1.d0/3.d0))/60.d0
    beta = (81.d0**(1.d0/3.d0)+9.d0**(1.d0/3.d0))*dsqrt(3.d0)/60.d0
    cno  = alph**2+beta**2
    u1   = 1.0d0/u1
    alph = alph/cno
    beta = beta/cno
    t11  = 9.1232394870892942792d-02
    t12  = -0.14125529502095420843d0
    t13  = -3.0029194105147424492d-02
    t21  = 0.24171793270710701896d0
    t22  = 0.20412935229379993199d0
    t23  = 0.38294211275726193779d0
    t31  = 0.96604818261509293619d0
    ti11 = 4.3255798900631553510d0
    ti12 = 0.33919925181580986954d0
    ti13 = 0.54177053993587487119d0
    ti21 = -4.1787185915519047273d0
    ti22 = -0.32768282076106238708d0
    ti23 = 0.47662355450055045196d0
    ti31 = -0.50287263494578687595d0
    ti32 = 2.5719269498556054292d0
    ti33 = -0.59603920482822492497d0
    posneg = sign(1.d0, xend-x)
    hmaxn  = min(abs(hmax), abs(xend-x))
    h = 1.0d-6                          ! original: H=0 input -> H=1e-6
    h = min(abs(h), hmaxn)
    h = sign(h, posneg)
    hold   = h
    reject = .false.
    first  = .true.
    last   = .false.
    if ((x+h*1.0001d0-xend)*posneg .ge. 0.d0) then
      h = xend-x
      last = .true.
    end if
    hopt   = h
    faccon = 1.d0
    cfac   = safe*(1+2*nit)
    nsing  = 0
    xold   = x
    rpar_d = 0.d0
    ipar_d = 0
    do i = 1, n
      scal(i) = atol(i)+rtol(i)*abs(y(i))
    end do
    hhfac = h
    call RD5_RHS(n, x, y, y0)
    nfcn = nfcn+1

    ! --- basic integration step
 10 continue
    ! *** computation of the Jacobian (IJAC=1: analytically)
    njac = njac+1
    call RD5_JAC(n, x, y, fjac, RD5_LD, rpar_d, ipar_d)
    caljac = .true.
 20 continue
    ! --- compute the matrices E1 and E2 and their decompositions
    fac1  = u1/h
    alphn = alph/h
    betan = beta/h
    ! DECOMR, IJOB=1: B=identity, Jacobian a full matrix
    do j = 1, n
      do i = 1, n
        e1(i,j) = -fjac(i,j)
      end do
      e1(j,j) = e1(j,j)+fac1
    end do
    call dec_dev(n, RD5_LD, e1, ip1, ier)
    if (ier .ne. 0) goto 78
    ! DECOMC, IJOB=1
    do j = 1, n
      do i = 1, n
        e2r(i,j) = -fjac(i,j)
        e2i(i,j) = 0.d0
      end do
      e2r(j,j) = e2r(j,j)+alphn
      e2i(j,j) = betan
    end do
    call decc_dev(n, RD5_LD, e2r, e2i, ip2, ier)
    if (ier .ne. 0) goto 78
 30 continue
    nstep = nstep+1
    if (nstep .gt. nmax) then
      idid = -2                          ! more than NMAX steps needed
      return
    end if
    if (0.1d0*abs(h) .le. abs(x)*uround) then
      idid = -3                          ! step size too small
      return
    end if
    xph = x+h
    ! *** starting values for Newton iteration
    if (first .or. startn) then
      do i = 1, n
        z1(i)=0.d0; z2(i)=0.d0; z3(i)=0.d0
        f1(i)=0.d0; f2(i)=0.d0; f3(i)=0.d0
      end do
    else
      c3q = h/hold
      c1q = c1*c3q
      c2q = c2*c3q
      do i = 1, n
        ak1 = cont(i+n)
        ak2 = cont(i+n2)
        ak3 = cont(i+n3)
        z1i = c1q*(ak1+(c1q-c2m1)*(ak2+(c1q-c1m1)*ak3))
        z2i = c2q*(ak1+(c2q-c2m1)*(ak2+(c2q-c1m1)*ak3))
        z3i = c3q*(ak1+(c3q-c2m1)*(ak2+(c3q-c1m1)*ak3))
        z1(i) = z1i
        z2(i) = z2i
        z3(i) = z3i
        f1(i) = ti11*z1i+ti12*z2i+ti13*z3i
        f2(i) = ti21*z1i+ti22*z2i+ti23*z3i
        f3(i) = ti31*z1i+ti32*z2i+ti33*z3i
      end do
    end if
    ! *** loop for the simplified Newton iteration
    newt = 0
    faccon = max(faccon, uround)**0.8d0
    theta = abs(thet)
 40 continue
    if (newt .ge. nit) goto 78
    ! --- compute the right-hand side
    do i = 1, n
      cont(i) = y(i)+z1(i)
    end do
    call RD5_RHS(n, x+c1*h, cont, z1)
    do i = 1, n
      cont(i) = y(i)+z2(i)
    end do
    call RD5_RHS(n, x+c2*h, cont, z2)
    do i = 1, n
      cont(i) = y(i)+z3(i)
    end do
    call RD5_RHS(n, xph, cont, z3)
    nfcn = nfcn+3
    ! --- solve the linear systems
    do i = 1, n
      a1 = z1(i)
      a2 = z2(i)
      a3 = z3(i)
      z1(i) = ti11*a1+ti12*a2+ti13*a3
      z2(i) = ti21*a1+ti22*a2+ti23*a3
      z3(i) = ti31*a1+ti32*a2+ti33*a3
    end do
    ! SLVRAD, IJOB=1 (IMAS=0)
    do i = 1, n
      a2 = -f2(i)
      a3 = -f3(i)
      z1(i) = z1(i)-f1(i)*fac1
      z2(i) = z2(i)+a2*alphn-a3*betan
      z3(i) = z3(i)+a3*alphn+a2*betan
    end do
    call sol_dev (n, RD5_LD, e1, z1, ip1)
    call solc_dev(n, RD5_LD, e2r, e2i, z2, z3, ip2)
    newt = newt+1
    dyno = 0.d0
    do i = 1, n
      denom = scal(i)
      dyno = dyno+(z1(i)/denom)**2+(z2(i)/denom)**2+(z3(i)/denom)**2
    end do
    dyno = dsqrt(dyno/n3)
    ! --- bad convergence or number of iterations too large
    if (newt.gt.1 .and. newt.lt.nit) then
      thq = dyno/dynold
      if (newt .eq. 2) then
        theta = thq
      else
        theta = sqrt(thq*thqold)
      end if
      thqold = thq
      if (theta .lt. 0.99d0) then
        faccon = theta/(1.0d0-theta)
        dyth = faccon*dyno*theta**(nit-1-newt)/fnewt
        if (dyth .ge. 1.0d0) then
          qnewt = dmax1(1.0d-4, dmin1(20.0d0, dyth))
          hhfac = .8d0*qnewt**(-1.0d0/(4.0d0+nit-1-newt))
          h = hhfac*h
          reject = .true.
          last = .false.
          if (caljac) goto 20
          goto 10
        end if
      else
        goto 78
      end if
    end if
    dynold = max(dyno, uround)
    do i = 1, n
      f1i = f1(i)+z1(i)
      f2i = f2(i)+z2(i)
      f3i = f3(i)+z3(i)
      f1(i) = f1i
      f2(i) = f2i
      f3(i) = f3i
      z1(i) = t11*f1i+t12*f2i+t13*f3i
      z2(i) = t21*f1i+t22*f2i+t23*f3i
      z3(i) = t31*f1i+    f2i
    end do
    if (faccon*dyno .gt. fnewt) goto 40
    ! --- error estimation (ESTRAD, IJOB=1, IMAS=0)
    hee1 = dd1/h
    hee2 = dd2/h
    hee3 = dd3/h
    do i = 1, n
      f2(i) = hee1*z1(i)+hee2*z2(i)+hee3*z3(i)
      cont(i) = f2(i)+y0(i)
    end do
    call sol_dev(n, RD5_LD, e1, cont, ip1)
    err = 0.d0
    do i = 1, n
      err = err+(cont(i)/scal(i))**2
    end do
    err = max(sqrt(err/n), 1.d-10)
    if (err .ge. 1.d0) then
      if (first .or. reject) then
        do i = 1, n
          cont(i) = y(i)+cont(i)
        end do
        call RD5_RHS(n, x, cont, f1)
        nfcn = nfcn+1
        do i = 1, n
          cont(i) = f1(i)+f2(i)
        end do
        call sol_dev(n, RD5_LD, e1, cont, ip1)
        err = 0.d0
        do i = 1, n
          err = err+(cont(i)/scal(i))**2
        end do
        err = max(sqrt(err/n), 1.d-10)
      end if
    end if
    ! --- computation of HNEW; we require .2<=HNEW/H<=8.
    fac  = min(safe, cfac/(newt+2*nit))
    quot = max(facr, min(facl, err**.25d0/fac))
    hnew = h/quot
    ! *** is the error small enough ?
    if (err .lt. 1.d0) then
      ! --- step is accepted
      first = .false.
      naccpt = naccpt+1
      if (pred) then
        ! --- predictive controller of Gustafsson
        if (naccpt .gt. 1) then
          facgus = (hacc/h)*(err**2/erracc)**0.25d0/safe
          facgus = max(facr, min(facl, facgus))
          quot = max(quot, facgus)
          hnew = h/quot
        end if
        hacc = h
        erracc = max(1.0d-2, err)
      end if
      xold = x
      hold = h
      x = xph
      do i = 1, n
        y(i) = y(i)+z3(i)
        z2i = z2(i)
        z1i = z1(i)
        cont(i+n) = (z2i-z3(i))/c2m1
        ak = (z1i-z2i)/c1mc2
        acont3 = z1i/c1
        acont3 = (ak-acont3)/c2
        cont(i+n2) = (ak-cont(i+n))/c1m1
        cont(i+n3) = cont(i+n2)-acont3
      end do
      do i = 1, n
        scal(i) = atol(i)+rtol(i)*abs(y(i))
      end do
      caljac = .false.
      if (last) then
        h = hopt
        idid = 1
        return
      end if
      call RD5_RHS(n, x, y, y0)
      nfcn = nfcn+1
      hnew = posneg*min(abs(hnew), hmaxn)
      hopt = hnew
      hopt = min(h, hnew)
      if (reject) hnew = posneg*min(abs(hnew), abs(h))
      reject = .false.
      if ((x+hnew/quot1-xend)*posneg .ge. 0.d0) then
        h = xend-x
        last = .true.
      else
        qt = hnew/h
        hhfac = h
        if (theta.le.thet .and. qt.ge.quot1 .and. qt.le.quot2) goto 30
        h = hnew
      end if
      hhfac = h
      if (theta .le. thet) goto 20
      goto 10
    else
      ! --- step is rejected
      reject = .true.
      last = .false.
      if (first) then
        h = h*0.1d0
        hhfac = 0.1d0
      else
        hhfac = hnew/h
        h = hnew
      end if
      if (naccpt .ge. 1) nrejct = nrejct+1
      if (caljac) goto 20
      goto 10
    end if
    ! --- unexpected step-rejection / singular matrix
 78 continue
    if (ier .ne. 0) then
      nsing = nsing+1
      if (nsing .ge. 5) then
        idid = -4                        ! matrix is repeatedly singular
        return
      end if
    end if
    h = h*0.5d0
    hhfac = 0.5d0
    reject = .true.
    last = .false.
    if (caljac) goto 20
    goto 10

  end subroutine radau5_dev


  !----------------------------------------------------------------------------
  ! DEC: LU triangularization with partial pivoting (Moler, ACM alg. 423).
  ! Verbatim port of hairerToolbox.f::DEC.
  !----------------------------------------------------------------------------
  subroutine dec_dev(n, ndim, a, ip, ier)
    !$acc routine seq
    implicit none
    integer, intent(in)    :: n, ndim
    real(8), intent(inout) :: a(ndim,n)
    integer, intent(out)   :: ip(n), ier
    integer :: nm1, k, kp1, m, i, j
    real(8) :: t

    ier = 0
    ip(n) = 1
    if (n .eq. 1) goto 70
    nm1 = n-1
    do k = 1, nm1
      kp1 = k+1
      m = k
      do i = kp1, n
        if (dabs(a(i,k)) .gt. dabs(a(m,k))) m = i
      end do
      ip(k) = m
      t = a(m,k)
      if (m .ne. k) then
        ip(n) = -ip(n)
        a(m,k) = a(k,k)
        a(k,k) = t
      end if
      if (t .eq. 0.d0) goto 80
      t = 1.d0/t
      do i = kp1, n
        a(i,k) = -a(i,k)*t
      end do
      do j = kp1, n
        t = a(m,j)
        a(m,j) = a(k,j)
        a(k,j) = t
        if (t .ne. 0.d0) then
          do i = kp1, n
            a(i,j) = a(i,j)+a(i,k)*t
          end do
        end if
      end do
    end do
 70 k = n
    if (a(n,n) .eq. 0.d0) goto 80
    return
 80 ier = k
    ip(n) = 0
  end subroutine dec_dev


  !----------------------------------------------------------------------------
  ! SOL: solve A*x=b from DEC factors. Verbatim port of hairerToolbox.f::SOL.
  !----------------------------------------------------------------------------
  subroutine sol_dev(n, ndim, a, b, ip)
    !$acc routine seq
    implicit none
    integer, intent(in)    :: n, ndim
    real(8), intent(in)    :: a(ndim,n)
    real(8), intent(inout) :: b(n)
    integer, intent(in)    :: ip(n)
    integer :: nm1, k, kp1, m, i, kb, km1
    real(8) :: t

    if (n .eq. 1) goto 50
    nm1 = n-1
    do k = 1, nm1
      kp1 = k+1
      m = ip(k)
      t = b(m)
      b(m) = b(k)
      b(k) = t
      do i = kp1, n
        b(i) = b(i)+a(i,k)*t
      end do
    end do
    do kb = 1, nm1
      km1 = n-kb
      k = km1+1
      b(k) = b(k)/a(k,k)
      t = -b(k)
      do i = 1, km1
        b(i) = b(i)+a(i,k)*t
      end do
    end do
 50 b(1) = b(1)/a(1,1)
  end subroutine sol_dev


  !----------------------------------------------------------------------------
  ! DECC: complex LU (real/imag pair storage). Verbatim port of
  ! hairerToolbox.f::DECC including the TI=0 / TR=0 fast paths.
  !----------------------------------------------------------------------------
  subroutine decc_dev(n, ndim, ar, ai, ip, ier)
    !$acc routine seq
    implicit none
    integer, intent(in)    :: n, ndim
    real(8), intent(inout) :: ar(ndim,n), ai(ndim,n)
    integer, intent(out)   :: ip(n), ier
    integer :: nm1, k, kp1, m, i, j
    real(8) :: tr, ti, den, prodr, prodi

    ier = 0
    ip(n) = 1
    if (n .eq. 1) goto 70
    nm1 = n-1
    do k = 1, nm1
      kp1 = k+1
      m = k
      do i = kp1, n
        if (dabs(ar(i,k))+dabs(ai(i,k)) .gt. dabs(ar(m,k))+dabs(ai(m,k))) m = i
      end do
      ip(k) = m
      tr = ar(m,k)
      ti = ai(m,k)
      if (m .ne. k) then
        ip(n) = -ip(n)
        ar(m,k) = ar(k,k)
        ai(m,k) = ai(k,k)
        ar(k,k) = tr
        ai(k,k) = ti
      end if
      if (dabs(tr)+dabs(ti) .eq. 0.d0) goto 80
      den = tr*tr+ti*ti
      tr = tr/den
      ti = -ti/den
      do i = kp1, n
        prodr = ar(i,k)*tr-ai(i,k)*ti
        prodi = ai(i,k)*tr+ar(i,k)*ti
        ar(i,k) = -prodr
        ai(i,k) = -prodi
      end do
      do j = kp1, n
        tr = ar(m,j)
        ti = ai(m,j)
        ar(m,j) = ar(k,j)
        ai(m,j) = ai(k,j)
        ar(k,j) = tr
        ai(k,j) = ti
        if (dabs(tr)+dabs(ti) .eq. 0.d0) then
          continue
        else if (ti .eq. 0.d0) then
          do i = kp1, n
            prodr = ar(i,k)*tr
            prodi = ai(i,k)*tr
            ar(i,j) = ar(i,j)+prodr
            ai(i,j) = ai(i,j)+prodi
          end do
        else if (tr .eq. 0.d0) then
          do i = kp1, n
            prodr = -ai(i,k)*ti
            prodi = ar(i,k)*ti
            ar(i,j) = ar(i,j)+prodr
            ai(i,j) = ai(i,j)+prodi
          end do
        else
          do i = kp1, n
            prodr = ar(i,k)*tr-ai(i,k)*ti
            prodi = ai(i,k)*tr+ar(i,k)*ti
            ar(i,j) = ar(i,j)+prodr
            ai(i,j) = ai(i,j)+prodi
          end do
        end if
      end do
    end do
 70 k = n
    if (dabs(ar(n,n))+dabs(ai(n,n)) .eq. 0.d0) goto 80
    return
 80 ier = k
    ip(n) = 0
  end subroutine decc_dev


  !----------------------------------------------------------------------------
  ! SOLC: complex solve from DECC factors. Verbatim port of
  ! hairerToolbox.f::SOLC.
  !----------------------------------------------------------------------------
  subroutine solc_dev(n, ndim, ar, ai, br, bi, ip)
    !$acc routine seq
    implicit none
    integer, intent(in)    :: n, ndim
    real(8), intent(in)    :: ar(ndim,n), ai(ndim,n)
    real(8), intent(inout) :: br(n), bi(n)
    integer, intent(in)    :: ip(n)
    integer :: nm1, k, kp1, m, i, kb, km1
    real(8) :: tr, ti, den, prodr, prodi

    if (n .eq. 1) goto 50
    nm1 = n-1
    do k = 1, nm1
      kp1 = k+1
      m = ip(k)
      tr = br(m)
      ti = bi(m)
      br(m) = br(k)
      bi(m) = bi(k)
      br(k) = tr
      bi(k) = ti
      do i = kp1, n
        prodr = ar(i,k)*tr-ai(i,k)*ti
        prodi = ai(i,k)*tr+ar(i,k)*ti
        br(i) = br(i)+prodr
        bi(i) = bi(i)+prodi
      end do
    end do
    do kb = 1, nm1
      km1 = n-kb
      k = km1+1
      den = ar(k,k)*ar(k,k)+ai(k,k)*ai(k,k)
      prodr = br(k)*ar(k,k)+bi(k)*ai(k,k)
      prodi = bi(k)*ar(k,k)-br(k)*ai(k,k)
      br(k) = prodr/den
      bi(k) = prodi/den
      tr = -br(k)
      ti = -bi(k)
      do i = 1, km1
        prodr = ar(i,k)*tr-ai(i,k)*ti
        prodi = ai(i,k)*tr+ar(i,k)*ti
        br(i) = br(i)+prodr
        bi(i) = bi(i)+prodi
      end do
    end do
 50 continue
    den = ar(1,1)*ar(1,1)+ai(1,1)*ai(1,1)
    prodr = br(1)*ar(1,1)+bi(1)*ai(1,1)
    prodi = bi(1)*ar(1,1)-br(1)*ai(1,1)
    br(1) = prodr/den
    bi(1) = prodi/den
  end subroutine solc_dev

# if defined (MOSE_CHEM_INLINE)
  ! ===========================================================================
  ! MOSE_CHEM_INLINE same-module copies (inlinable). Bodies VERBATIM from
  ! Lib_Chemistry_data.f90 (f_kf/f_kb), Lib_ChemMech/Frolov_nopressure.f90,
  ! Lib_ChemMech/ONERA-7.f90 and Lib_Chemistry_rhs.f90 (rhs/jac_native).
  ! KEEP IN SYNC with the originals.
  ! ===========================================================================

  pure function f_kf_m(ireact,Tint,Tdiff) result(result)
    !$acc routine seq
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
    a = kf_tab(Tint(1),ireact)
    b = kf_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff
  end function f_kf_m

  pure function f_kb_m(ireact,Tint,Tdiff) result(result)
    !$acc routine seq
    implicit none
    integer, intent(in) :: ireact, Tint(2)
    real(8), intent(in) :: Tdiff
    real(8) :: a, b
    real(8) :: result
    a = kb_tab(Tint(1),ireact)
    b = kb_tab(Tint(2),ireact)
    result = a+(b-a)*Tdiff
  end function f_kb_m

  subroutine frolov_m(roi, temp, omegadot)
    !$acc routine seq
    implicit none
    real(8), intent(inout) :: roi(ns)
    real(8), intent(in)    :: temp
    real(8), intent(out)   :: omegadot(ns)
    real(8) :: coi(NSCHEM), Tdiff
    integer :: is, T_i, Tint(2)
    real(8) :: kf, kb, net_rate

    do is = 1, ns
      coi(is) = roi(is) / Wm_tab(is)   ! kmol/m^3
    end do

    T_i     = int(temp)
    Tdiff   = temp - T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1

    kf = f_kf_m(1, Tint, Tdiff)
    kb = f_kb_m(1, Tint, Tdiff)

    ! Reaction 1:  2 H2 + O2 <=> 2 H2O
    net_rate = kf * coi(3) * coi(3) * coi(1) - kb * coi(2) * coi(2)

    omegadot(1) = -      Wm_tab(1) * net_rate
    omegadot(2) =  2.d0 * Wm_tab(2) * net_rate
    omegadot(3) = -2.d0 * Wm_tab(3) * net_rate
    omegadot(4) = 0.d0
  end subroutine frolov_m

  subroutine frolov_jac_m(roi, temp, dwdr, dwdT)
    !$acc routine seq
    implicit none
    real(8), intent(in)  :: roi(ns), temp
    real(8), intent(out) :: dwdr(NSCHEM, ns)
    real(8), intent(out) :: dwdT(ns)
    real(8) :: coi(NSCHEM), Tdiff
    integer :: is, T_i, j
    real(8) :: kf, kb, dkf_dT, dkb_dT
    real(8) :: dnet_dc(NSCHEM), dnet_dT_r
    real(8) :: dwdr_c(NSCHEM, NSCHEM)

    do is = 1, ns
      coi(is) = roi(is) / Wm_tab(is)
    end do

    T_i   = int(temp)
    Tdiff = temp - T_i

    kf     = kf_tab(T_i  , 1) + Tdiff * (kf_tab(T_i+1, 1) - kf_tab(T_i, 1))
    kb     = kb_tab(T_i  , 1) + Tdiff * (kb_tab(T_i+1, 1) - kb_tab(T_i, 1))
    dkf_dT = kf_tab(T_i+1, 1) - kf_tab(T_i, 1)
    dkb_dT = kb_tab(T_i+1, 1) - kb_tab(T_i, 1)

    dnet_dc(1:ns) = 0.d0
    dnet_dc(1) =          kf * coi(3) * coi(3)
    dnet_dc(2) = -2.d0 *  kb * coi(2)
    dnet_dc(3) =  2.d0 *  kf * coi(3) * coi(1)
    dnet_dT_r  = dkf_dT * coi(3) * coi(3) * coi(1) &
               - dkb_dT * coi(2) * coi(2)

    dwdr_c(1:ns,1:ns) = 0.d0
    dwdr_c(1, 1:ns) = -      Wm_tab(1) * dnet_dc(1:ns)
    dwdr_c(2, 1:ns) =  2.d0 * Wm_tab(2) * dnet_dc(1:ns)
    dwdr_c(3, 1:ns) = -2.d0 * Wm_tab(3) * dnet_dc(1:ns)

    dwdT(:) = 0.d0
    dwdT(1) = -      Wm_tab(1) * dnet_dT_r
    dwdT(2) =  2.d0 * Wm_tab(2) * dnet_dT_r
    dwdT(3) = -2.d0 * Wm_tab(3) * dnet_dT_r

    do j = 1, ns
      dwdr(1:ns, j) = dwdr_c(1:ns, j) / Wm_tab(j)
    end do
  end subroutine frolov_jac_m

  subroutine onera7_m(roi,temp,omegadot)
    !$acc routine seq
    implicit none
    real(8), intent(inout)  :: roi(ns)
    real(8), intent(in)  :: temp
    real(8), intent(out) :: omegadot(ns)
    real(8) :: coi(NSCHEM), Tdiff
    real(8) :: M
    integer :: is, T_i, Tint(2)
    real(8) :: prodf(1:14), prodb(1:14)

    do is = 1, ns
     coi(is)=roi(is)/Wm_tab(is) ! kmol/m^3
    enddo
    T_i = int(temp)
    Tdiff  = temp-T_i
    Tint(1) = T_i
    Tint(2) = T_i + 1
    ! reac n. 1: H2 + O2 => 2 OH
    prodf(1)=f_kf_m(1,Tint,Tdiff)*(coi(3))*(coi(1))
    prodb(1)=f_kb_m(1,Tint,Tdiff)*(coi(6)*coi(6))
    ! reac n. 2: 2 OH => H2 + O2
    prodf(2)=f_kf_m(2,Tint,Tdiff)*(coi(6)*coi(6))
    prodb(2)=f_kb_m(2,Tint,Tdiff)*(coi(3))*(coi(1))
    ! reac n. 3: H + O2 => O + OH
    prodf(3)=f_kf_m(3,Tint,Tdiff)*(coi(4))*(coi(1))
    prodb(3)=f_kb_m(3,Tint,Tdiff)*(coi(5))*(coi(6))
    ! reac n. 4: O + OH => H + O2
    prodf(4)=f_kf_m(4,Tint,Tdiff)*(coi(5))*(coi(6))
    prodb(4)=f_kb_m(4,Tint,Tdiff)*(coi(4))*(coi(1))
    ! reac n. 5: H2 + OH => H + H2O
    prodf(5)=f_kf_m(5,Tint,Tdiff)*(coi(3))*(coi(6))
    prodb(5)=f_kb_m(5,Tint,Tdiff)*(coi(4))*(coi(2))
    ! reac n. 6: H + H2O => H2 + OH
    prodf(6)=f_kf_m(6,Tint,Tdiff)*(coi(4))*(coi(2))
    prodb(6)=f_kb_m(6,Tint,Tdiff)*(coi(3))*(coi(6))
    ! reac n. 7: H2 + O => H + OH
    prodf(7)=f_kf_m(7,Tint,Tdiff)*(coi(3))*(coi(5))
    prodb(7)=f_kb_m(7,Tint,Tdiff)*(coi(4))*(coi(6))
    ! reac n. 8: H + OH => H2 + O
    prodf(8)=f_kf_m(8,Tint,Tdiff)*(coi(4))*(coi(6))
    prodb(8)=f_kb_m(8,Tint,Tdiff)*(coi(3))*(coi(5))
    ! reac n. 9: 2 OH => H2O + O
    prodf(9)=f_kf_m(9,Tint,Tdiff)*(coi(6)*coi(6))
    prodb(9)=f_kb_m(9,Tint,Tdiff)*(coi(2))*(coi(5))
    ! reac n. 10: H2O + O => 2 OH
    prodf(10)=f_kf_m(10,Tint,Tdiff)*(coi(2))*(coi(5))
    prodb(10)=f_kb_m(10,Tint,Tdiff)*(coi(6)*coi(6))
    ! reac n. 11: H + OH + M => H2O + M
    M=coi(1)+coi(2)*12+coi(3)*2.5+coi(4)+coi(5)+coi(6)+coi(7)
    prodf(11)=f_kf_m(11,Tint,Tdiff)*(coi(4))*(coi(6))*M
    prodb(11)=f_kb_m(11,Tint,Tdiff)*(coi(2))*M
    ! reac n. 12: H2O + M => H + OH + M
    M=coi(1)+coi(2)*12+coi(3)*2.5+coi(4)+coi(5)+coi(6)+coi(7)
    prodf(12)=f_kf_m(12,Tint,Tdiff)*(coi(2))*M
    prodb(12)=f_kb_m(12,Tint,Tdiff)*(coi(4))*(coi(6))*M
    ! reac n. 13: 2 H + M => H2 + M
    M=coi(1)+coi(2)*12+coi(3)*2.5+coi(4)+coi(5)+coi(6)+coi(7)
    prodf(13)=f_kf_m(13,Tint,Tdiff)*(coi(4)*coi(4))*M
    prodb(13)=f_kb_m(13,Tint,Tdiff)*(coi(3))*M
    ! reac n. 14: H2 + M => 2 H + M
    M=coi(1)+coi(2)*12+coi(3)*2.5+coi(4)+coi(5)+coi(6)+coi(7)
    prodf(14)=f_kf_m(14,Tint,Tdiff)*(coi(3))*M
    prodb(14)=f_kb_m(14,Tint,Tdiff)*(coi(4)*coi(4))*M
    ! species source terms
    omegadot(1)=Wm_tab(1)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(4)-prodb(4)))
    omegadot(2)=Wm_tab(2)*(+(1.0-0.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(10)-prodb(10))+(1.0-0.0)*(prodf(11)-prodb(11))+(0.0-1.0)*(prodf(12)-prodb(12)))
    omegadot(3)=Wm_tab(3)*(+(0.0-1.0)*(prodf(1)-prodb(1))+(1.0-0.0)*(prodf(2)-prodb(2))+(0.0-1.0)*(prodf(5)-prodb(5))+(1.0-0.0)*(prodf(6)-prodb(6))+(0.0-1.0)*(prodf(7)-prodb(7))+(1.0-0.0)*(prodf(8)-prodb(8))+(1.0-0.0)*(prodf(13)-prodb(13))+(0.0-1.0)*(prodf(14)-prodb(14)))
    omegadot(4)=Wm_tab(4)*(+(0.0-1.0)*(prodf(3)-prodb(3))+(1.0-0.0)*(prodf(4)-prodb(4))+(1.0-0.0)*(prodf(5)-prodb(5))+(0.0-1.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(7)-prodb(7))+(0.0-1.0)*(prodf(8)-prodb(8))+(0.0-1.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(12)-prodb(12))+(0.0-2.0)*(prodf(13)-prodb(13))+(2.0-0.0)*(prodf(14)-prodb(14)))
    omegadot(5)=Wm_tab(5)*(+(1.0-0.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(4)-prodb(4))+(0.0-1.0)*(prodf(7)-prodb(7))+(1.0-0.0)*(prodf(8)-prodb(8))+(1.0-0.0)*(prodf(9)-prodb(9))+(0.0-1.0)*(prodf(10)-prodb(10)))
    omegadot(6)=Wm_tab(6)*(+(2.0-0.0)*(prodf(1)-prodb(1))+(0.0-2.0)*(prodf(2)-prodb(2))+(1.0-0.0)*(prodf(3)-prodb(3))+(0.0-1.0)*(prodf(4)-prodb(4))+(0.0-1.0)*(prodf(5)-prodb(5))+(1.0-0.0)*(prodf(6)-prodb(6))+(1.0-0.0)*(prodf(7)-prodb(7))+(0.0-1.0)*(prodf(8)-prodb(8))+(0.0-2.0)*(prodf(9)-prodb(9))+(2.0-0.0)*(prodf(10)-prodb(10))+(0.0-1.0)*(prodf(11)-prodb(11))+(1.0-0.0)*(prodf(12)-prodb(12)))
    omegadot(7)=0.d0
  end subroutine onera7_m

  subroutine onera7_jac_m(roi, temp, dwdr, dwdT)
    !$acc routine seq
    implicit none
    real(8), intent(in)  :: roi(ns), temp
    real(8), intent(out) :: dwdr(NSCHEM,ns)
    real(8), intent(out) :: dwdT(ns)
    real(8) :: coi(NSCHEM), Tdiff
    integer :: T_i, j
    real(8) :: kf_r(14), kb_r(14)
    real(8) :: dkf_dT(14), dkb_dT(14)
    real(8) :: M, dM_dc(NSCHEM)
    real(8) :: dRf_dc(NSCHEM), dRb_dc(NSCHEM), dnet_dc(NSCHEM)
    real(8) :: dRf_dT_r, dRb_dT_r, dnet_dT_r
    real(8) :: dwdr_c(NSCHEM,NSCHEM)
    real(8), parameter :: epsM(7) = [1.d0, 12.d0, 2.5d0, 1.d0, 1.d0, 1.d0, 1.d0]

    do j = 1, ns
      coi(j) = roi(j)/Wm_tab(j)
    enddo
    T_i   = int(temp)
    Tdiff = temp - T_i

    do j = 1, 14
      kf_r(j)   =  kf_tab(T_i  ,j) + Tdiff*(kf_tab(T_i+1,j) - kf_tab(T_i,j))
      kb_r(j)   =  kb_tab(T_i  ,j) + Tdiff*(kb_tab(T_i+1,j) - kb_tab(T_i,j))
      dkf_dT(j) =  kf_tab(T_i+1,j) - kf_tab(T_i,j)
      dkb_dT(j) =  kb_tab(T_i+1,j) - kb_tab(T_i,j)
    enddo

    M = coi(1) + 12.d0*coi(2) + 2.5d0*coi(3) + coi(4) + coi(5) + coi(6) + coi(7)
    dM_dc(:) = epsM(:)

    dwdr_c = 0.d0
    dwdT   = 0.d0

    ! ---- reac 1: H2(3) + O2(1) => 2 OH(6)
    dRf_dc = 0.d0; dRb_dc = 0.d0
    dRf_dc(1) = kf_r(1)*coi(3);   dRf_dc(3) = kf_r(1)*coi(1)
    dRb_dc(6) = 2.d0*kb_r(1)*coi(6)
    dRf_dT_r = dkf_dT(1)*coi(3)*coi(1)
    dRb_dT_r = dkb_dT(1)*coi(6)*coi(6)
    dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
    dwdr_c(1,:) = dwdr_c(1,:) -      Wm_tab(1)*dnet_dc
    dwdr_c(3,:) = dwdr_c(3,:) -      Wm_tab(3)*dnet_dc
    dwdr_c(6,:) = dwdr_c(6,:) + 2.d0*Wm_tab(6)*dnet_dc
    dwdT(1) = dwdT(1) -      Wm_tab(1)*dnet_dT_r
    dwdT(3) = dwdT(3) -      Wm_tab(3)*dnet_dT_r
    dwdT(6) = dwdT(6) + 2.d0*Wm_tab(6)*dnet_dT_r

    ! ---- reac 2: 2 OH => H2 + O2
    dRf_dc = 0.d0; dRb_dc = 0.d0
    dRf_dc(6) = 2.d0*kf_r(2)*coi(6)
    dRb_dc(1) = kb_r(2)*coi(3);   dRb_dc(3) = kb_r(2)*coi(1)
    dRf_dT_r = dkf_dT(2)*coi(6)*coi(6)
    dRb_dT_r = dkb_dT(2)*coi(3)*coi(1)
    dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
    dwdr_c(1,:) = dwdr_c(1,:) +      Wm_tab(1)*dnet_dc
    dwdr_c(3,:) = dwdr_c(3,:) +      Wm_tab(3)*dnet_dc
    dwdr_c(6,:) = dwdr_c(6,:) - 2.d0*Wm_tab(6)*dnet_dc
    dwdT(1) = dwdT(1) +      Wm_tab(1)*dnet_dT_r
    dwdT(3) = dwdT(3) +      Wm_tab(3)*dnet_dT_r
    dwdT(6) = dwdT(6) - 2.d0*Wm_tab(6)*dnet_dT_r

    ! ---- reac 3: H(4) + O2(1) => O(5) + OH(6)
    dRf_dc = 0.d0; dRb_dc = 0.d0
    dRf_dc(1) = kf_r(3)*coi(4);   dRf_dc(4) = kf_r(3)*coi(1)
    dRb_dc(5) = kb_r(3)*coi(6);   dRb_dc(6) = kb_r(3)*coi(5)
    dRf_dT_r = dkf_dT(3)*coi(4)*coi(1)
    dRb_dT_r = dkb_dT(3)*coi(5)*coi(6)
    dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
    dwdr_c(1,:) = dwdr_c(1,:) - Wm_tab(1)*dnet_dc
    dwdr_c(4,:) = dwdr_c(4,:) - Wm_tab(4)*dnet_dc
    dwdr_c(5,:) = dwdr_c(5,:) + Wm_tab(5)*dnet_dc
    dwdr_c(6,:) = dwdr_c(6,:) + Wm_tab(6)*dnet_dc
    dwdT(1) = dwdT(1) - Wm_tab(1)*dnet_dT_r
    dwdT(4) = dwdT(4) - Wm_tab(4)*dnet_dT_r
    dwdT(5) = dwdT(5) + Wm_tab(5)*dnet_dT_r
    dwdT(6) = dwdT(6) + Wm_tab(6)*dnet_dT_r

    ! ---- reac 4: O(5) + OH(6) => H(4) + O2(1)
    dRf_dc = 0.d0; dRb_dc = 0.d0
    dRf_dc(5) = kf_r(4)*coi(6);   dRf_dc(6) = kf_r(4)*coi(5)
    dRb_dc(1) = kb_r(4)*coi(4);   dRb_dc(4) = kb_r(4)*coi(1)
    dRf_dT_r = dkf_dT(4)*coi(5)*coi(6)
    dRb_dT_r = dkb_dT(4)*coi(4)*coi(1)
    dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
    dwdr_c(1,:) = dwdr_c(1,:) + Wm_tab(1)*dnet_dc
    dwdr_c(4,:) = dwdr_c(4,:) + Wm_tab(4)*dnet_dc
    dwdr_c(5,:) = dwdr_c(5,:) - Wm_tab(5)*dnet_dc
    dwdr_c(6,:) = dwdr_c(6,:) - Wm_tab(6)*dnet_dc
    dwdT(1) = dwdT(1) + Wm_tab(1)*dnet_dT_r
    dwdT(4) = dwdT(4) + Wm_tab(4)*dnet_dT_r
    dwdT(5) = dwdT(5) - Wm_tab(5)*dnet_dT_r
    dwdT(6) = dwdT(6) - Wm_tab(6)*dnet_dT_r

    ! ---- reac 5: H2(3) + OH(6) => H(4) + H2O(2)
    dRf_dc = 0.d0; dRb_dc = 0.d0
    dRf_dc(3) = kf_r(5)*coi(6);   dRf_dc(6) = kf_r(5)*coi(3)
    dRb_dc(2) = kb_r(5)*coi(4);   dRb_dc(4) = kb_r(5)*coi(2)
    dRf_dT_r = dkf_dT(5)*coi(3)*coi(6)
    dRb_dT_r = dkb_dT(5)*coi(4)*coi(2)
    dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
    dwdr_c(2,:) = dwdr_c(2,:) + Wm_tab(2)*dnet_dc
    dwdr_c(3,:) = dwdr_c(3,:) - Wm_tab(3)*dnet_dc
    dwdr_c(4,:) = dwdr_c(4,:) + Wm_tab(4)*dnet_dc
    dwdr_c(6,:) = dwdr_c(6,:) - Wm_tab(6)*dnet_dc
    dwdT(2) = dwdT(2) + Wm_tab(2)*dnet_dT_r
    dwdT(3) = dwdT(3) - Wm_tab(3)*dnet_dT_r
    dwdT(4) = dwdT(4) + Wm_tab(4)*dnet_dT_r
    dwdT(6) = dwdT(6) - Wm_tab(6)*dnet_dT_r

    ! ---- reac 6: H(4) + H2O(2) => H2(3) + OH(6)
    dRf_dc = 0.d0; dRb_dc = 0.d0
    dRf_dc(2) = kf_r(6)*coi(4);   dRf_dc(4) = kf_r(6)*coi(2)
    dRb_dc(3) = kb_r(6)*coi(6);   dRb_dc(6) = kb_r(6)*coi(3)
    dRf_dT_r = dkf_dT(6)*coi(4)*coi(2)
    dRb_dT_r = dkb_dT(6)*coi(3)*coi(6)
    dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
    dwdr_c(2,:) = dwdr_c(2,:) - Wm_tab(2)*dnet_dc
    dwdr_c(3,:) = dwdr_c(3,:) + Wm_tab(3)*dnet_dc
    dwdr_c(4,:) = dwdr_c(4,:) - Wm_tab(4)*dnet_dc
    dwdr_c(6,:) = dwdr_c(6,:) + Wm_tab(6)*dnet_dc
    dwdT(2) = dwdT(2) - Wm_tab(2)*dnet_dT_r
    dwdT(3) = dwdT(3) + Wm_tab(3)*dnet_dT_r
    dwdT(4) = dwdT(4) - Wm_tab(4)*dnet_dT_r
    dwdT(6) = dwdT(6) + Wm_tab(6)*dnet_dT_r

    ! ---- reac 7: H2(3) + O(5) => H(4) + OH(6)
    dRf_dc = 0.d0; dRb_dc = 0.d0
    dRf_dc(3) = kf_r(7)*coi(5);   dRf_dc(5) = kf_r(7)*coi(3)
    dRb_dc(4) = kb_r(7)*coi(6);   dRb_dc(6) = kb_r(7)*coi(4)
    dRf_dT_r = dkf_dT(7)*coi(3)*coi(5)
    dRb_dT_r = dkb_dT(7)*coi(4)*coi(6)
    dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
    dwdr_c(3,:) = dwdr_c(3,:) - Wm_tab(3)*dnet_dc
    dwdr_c(4,:) = dwdr_c(4,:) + Wm_tab(4)*dnet_dc
    dwdr_c(5,:) = dwdr_c(5,:) - Wm_tab(5)*dnet_dc
    dwdr_c(6,:) = dwdr_c(6,:) + Wm_tab(6)*dnet_dc
    dwdT(3) = dwdT(3) - Wm_tab(3)*dnet_dT_r
    dwdT(4) = dwdT(4) + Wm_tab(4)*dnet_dT_r
    dwdT(5) = dwdT(5) - Wm_tab(5)*dnet_dT_r
    dwdT(6) = dwdT(6) + Wm_tab(6)*dnet_dT_r

    ! ---- reac 8: H(4) + OH(6) => H2(3) + O(5)
    dRf_dc = 0.d0; dRb_dc = 0.d0
    dRf_dc(4) = kf_r(8)*coi(6);   dRf_dc(6) = kf_r(8)*coi(4)
    dRb_dc(3) = kb_r(8)*coi(5);   dRb_dc(5) = kb_r(8)*coi(3)
    dRf_dT_r = dkf_dT(8)*coi(4)*coi(6)
    dRb_dT_r = dkb_dT(8)*coi(3)*coi(5)
    dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
    dwdr_c(3,:) = dwdr_c(3,:) + Wm_tab(3)*dnet_dc
    dwdr_c(4,:) = dwdr_c(4,:) - Wm_tab(4)*dnet_dc
    dwdr_c(5,:) = dwdr_c(5,:) + Wm_tab(5)*dnet_dc
    dwdr_c(6,:) = dwdr_c(6,:) - Wm_tab(6)*dnet_dc
    dwdT(3) = dwdT(3) + Wm_tab(3)*dnet_dT_r
    dwdT(4) = dwdT(4) - Wm_tab(4)*dnet_dT_r
    dwdT(5) = dwdT(5) + Wm_tab(5)*dnet_dT_r
    dwdT(6) = dwdT(6) - Wm_tab(6)*dnet_dT_r

    ! ---- reac 9: 2 OH(6) => H2O(2) + O(5)
    dRf_dc = 0.d0; dRb_dc = 0.d0
    dRf_dc(6) = 2.d0*kf_r(9)*coi(6)
    dRb_dc(2) = kb_r(9)*coi(5);   dRb_dc(5) = kb_r(9)*coi(2)
    dRf_dT_r = dkf_dT(9)*coi(6)*coi(6)
    dRb_dT_r = dkb_dT(9)*coi(2)*coi(5)
    dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
    dwdr_c(2,:) = dwdr_c(2,:) +      Wm_tab(2)*dnet_dc
    dwdr_c(5,:) = dwdr_c(5,:) +      Wm_tab(5)*dnet_dc
    dwdr_c(6,:) = dwdr_c(6,:) - 2.d0*Wm_tab(6)*dnet_dc
    dwdT(2) = dwdT(2) +      Wm_tab(2)*dnet_dT_r
    dwdT(5) = dwdT(5) +      Wm_tab(5)*dnet_dT_r
    dwdT(6) = dwdT(6) - 2.d0*Wm_tab(6)*dnet_dT_r

    ! ---- reac 10: H2O(2) + O(5) => 2 OH(6)
    dRf_dc = 0.d0; dRb_dc = 0.d0
    dRf_dc(2) = kf_r(10)*coi(5);  dRf_dc(5) = kf_r(10)*coi(2)
    dRb_dc(6) = 2.d0*kb_r(10)*coi(6)
    dRf_dT_r = dkf_dT(10)*coi(2)*coi(5)
    dRb_dT_r = dkb_dT(10)*coi(6)*coi(6)
    dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
    dwdr_c(2,:) = dwdr_c(2,:) -      Wm_tab(2)*dnet_dc
    dwdr_c(5,:) = dwdr_c(5,:) -      Wm_tab(5)*dnet_dc
    dwdr_c(6,:) = dwdr_c(6,:) + 2.d0*Wm_tab(6)*dnet_dc
    dwdT(2) = dwdT(2) -      Wm_tab(2)*dnet_dT_r
    dwdT(5) = dwdT(5) -      Wm_tab(5)*dnet_dT_r
    dwdT(6) = dwdT(6) + 2.d0*Wm_tab(6)*dnet_dT_r

    ! ---- reac 11: H(4) + OH(6) + M => H2O(2) + M
    dRf_dc(:) = kf_r(11)*coi(4)*coi(6)*dM_dc(:)
    dRf_dc(4) = dRf_dc(4) + kf_r(11)*coi(6)*M
    dRf_dc(6) = dRf_dc(6) + kf_r(11)*coi(4)*M
    dRb_dc(:) = kb_r(11)*coi(2)*dM_dc(:)
    dRb_dc(2) = dRb_dc(2) + kb_r(11)*M
    dRf_dT_r = dkf_dT(11)*coi(4)*coi(6)*M
    dRb_dT_r = dkb_dT(11)*coi(2)*M
    dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
    dwdr_c(2,:) = dwdr_c(2,:) + Wm_tab(2)*dnet_dc
    dwdr_c(4,:) = dwdr_c(4,:) - Wm_tab(4)*dnet_dc
    dwdr_c(6,:) = dwdr_c(6,:) - Wm_tab(6)*dnet_dc
    dwdT(2) = dwdT(2) + Wm_tab(2)*dnet_dT_r
    dwdT(4) = dwdT(4) - Wm_tab(4)*dnet_dT_r
    dwdT(6) = dwdT(6) - Wm_tab(6)*dnet_dT_r

    ! ---- reac 12: H2O(2) + M => H(4) + OH(6) + M
    dRf_dc(:) = kf_r(12)*coi(2)*dM_dc(:)
    dRf_dc(2) = dRf_dc(2) + kf_r(12)*M
    dRb_dc(:) = kb_r(12)*coi(4)*coi(6)*dM_dc(:)
    dRb_dc(4) = dRb_dc(4) + kb_r(12)*coi(6)*M
    dRb_dc(6) = dRb_dc(6) + kb_r(12)*coi(4)*M
    dRf_dT_r = dkf_dT(12)*coi(2)*M
    dRb_dT_r = dkb_dT(12)*coi(4)*coi(6)*M
    dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
    dwdr_c(2,:) = dwdr_c(2,:) - Wm_tab(2)*dnet_dc
    dwdr_c(4,:) = dwdr_c(4,:) + Wm_tab(4)*dnet_dc
    dwdr_c(6,:) = dwdr_c(6,:) + Wm_tab(6)*dnet_dc
    dwdT(2) = dwdT(2) - Wm_tab(2)*dnet_dT_r
    dwdT(4) = dwdT(4) + Wm_tab(4)*dnet_dT_r
    dwdT(6) = dwdT(6) + Wm_tab(6)*dnet_dT_r

    ! ---- reac 13: 2 H(4) + M => H2(3) + M
    dRf_dc(:) = kf_r(13)*coi(4)*coi(4)*dM_dc(:)
    dRf_dc(4) = dRf_dc(4) + 2.d0*kf_r(13)*coi(4)*M
    dRb_dc(:) = kb_r(13)*coi(3)*dM_dc(:)
    dRb_dc(3) = dRb_dc(3) + kb_r(13)*M
    dRf_dT_r = dkf_dT(13)*coi(4)*coi(4)*M
    dRb_dT_r = dkb_dT(13)*coi(3)*M
    dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
    dwdr_c(3,:) = dwdr_c(3,:) +      Wm_tab(3)*dnet_dc
    dwdr_c(4,:) = dwdr_c(4,:) - 2.d0*Wm_tab(4)*dnet_dc
    dwdT(3) = dwdT(3) +      Wm_tab(3)*dnet_dT_r
    dwdT(4) = dwdT(4) - 2.d0*Wm_tab(4)*dnet_dT_r

    ! ---- reac 14: H2(3) + M => 2 H(4) + M
    dRf_dc(:) = kf_r(14)*coi(3)*dM_dc(:)
    dRf_dc(3) = dRf_dc(3) + kf_r(14)*M
    dRb_dc(:) = kb_r(14)*coi(4)*coi(4)*dM_dc(:)
    dRb_dc(4) = dRb_dc(4) + 2.d0*kb_r(14)*coi(4)*M
    dRf_dT_r = dkf_dT(14)*coi(3)*M
    dRb_dT_r = dkb_dT(14)*coi(4)*coi(4)*M
    dnet_dc = dRf_dc - dRb_dc;  dnet_dT_r = dRf_dT_r - dRb_dT_r
    dwdr_c(3,:) = dwdr_c(3,:) -      Wm_tab(3)*dnet_dc
    dwdr_c(4,:) = dwdr_c(4,:) + 2.d0*Wm_tab(4)*dnet_dc
    dwdT(3) = dwdT(3) -      Wm_tab(3)*dnet_dT_r
    dwdT(4) = dwdT(4) + 2.d0*Wm_tab(4)*dnet_dT_r

    do j = 1, ns
      dwdr(:,j) = dwdr_c(:,j) / Wm_tab(j)
    enddo
  end subroutine onera7_jac_m

  subroutine rhs_m ( nz, time, Z, F )
    !$acc routine seq
    implicit none
    integer, intent(in)  :: nz
    real(8), intent(in)  :: time
    real(8), intent(in)  :: Z(nz)
    real(8), intent(out) :: F(nz)
    real(8) :: roi(NSCHEM)
    real(8) :: T, T_frac
    real(8) :: droic(NSCHEM)
    real(8) :: eiroi,rho_cv
    integer :: s, T_idx
    real(8) :: h_val, cp_val

    T = Z(nz)

    if (T < Tmin .or. T >= Tmax .or. ieee_is_nan(T)) then
       F(:) = -1.0d0
       return
    end if

    T_idx  = idint(T)
    T_frac = T - dble(T_idx)

    roi(1:ns) = max(Z(1:ns), 0.d0)

    if (mech_dev == 2) then
      call frolov_m ( roi, T, droic )
    else
      call onera7_m ( roi, T, droic )
    end if

    eiroi = 0.d0; rho_cv = 0.d0
    do s = 1, ns
        h_val  = h_tab(T_idx, s) + T_frac * (h_tab(T_idx+1, s) - h_tab(T_idx, s))
        cp_val = cp_tab(T_idx, s) + T_frac * (cp_tab(T_idx+1, s) - cp_tab(T_idx, s))

        eiroi  = eiroi + (h_val - Ri_tab(s) * T) * droic(s)
        rho_cv = rho_cv + roi(s) * (cp_val - Ri_tab(s))
    end do

    F(1:ns) = droic(1:ns)
    F(nz) = -eiroi / rho_cv
  end subroutine rhs_m

  subroutine jac_m(nz, time, Z, DFY, LDFY, RPAR, IPAR)
    !$acc routine seq
    implicit none
    integer, intent(in)  :: nz, LDFY
    real(8), intent(in)  :: time
    real(8), intent(in)  :: Z(nz)
    real(8), intent(out) :: DFY(LDFY, nz)
    real(8), intent(in)  :: RPAR(*)
    integer, intent(in)  :: IPAR(*)
    real(8) :: roi(NSCHEM), droic(NSCHEM)
    real(8) :: dwdr(NSCHEM, NSCHEM), dwdT(NSCHEM)
    real(8) :: h_vec(NSCHEM), cp_vec(NSCHEM), dh_dT(NSCHEM), dcp_dT(NSCHEM)
    real(8) :: T, T_frac
    integer :: T_idx, s, j
    real(8) :: G, D, inv_D, F_T
    real(8) :: dG_drho_j, dD_drho_j, dG_dT, dD_dT

    T = Z(nz)

    if (T < Tmin .or. T >= Tmax .or. ieee_is_nan(T)) then
      DFY(1:nz, 1:nz) = 0.d0
      return
    end if

    T_idx  = idint(T)
    T_frac = T - dble(T_idx)
    roi(1:ns) = max(Z(1:ns), 0.d0)

    if (mech_dev == 2) then
      call frolov_m     ( roi, T, droic )
      call frolov_jac_m ( roi, T, dwdr, dwdT )
    else
      call onera7_m     ( roi, T, droic )
      call onera7_jac_m ( roi, T, dwdr, dwdT )
    end if

    do s = 1, ns
      h_vec(s)      = h_tab(T_idx,s)  + T_frac * (h_tab(T_idx+1,s)  - h_tab(T_idx,s))
      cp_vec(s)     = cp_tab(T_idx,s) + T_frac * (cp_tab(T_idx+1,s) - cp_tab(T_idx,s))
      dh_dT(s)  = h_tab(T_idx+1,s)  - h_tab(T_idx,s)
      dcp_dT(s) = cp_tab(T_idx+1,s) - cp_tab(T_idx,s)
    end do

    G = 0.d0
    D = 0.d0
    do s = 1, ns
      G = G + (h_vec(s) - Ri_tab(s) * T) * droic(s)
      D = D + roi(s) * (cp_vec(s) - Ri_tab(s))
    end do
    inv_D = 1.d0 / D
    F_T   = -G * inv_D

    do j = 1, ns
      DFY(1:ns, j) = dwdr(1:ns, j)
    end do
    DFY(1:ns, nz) = dwdT(1:ns)

    do j = 1, ns
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

  end subroutine jac_m
# endif

end module FLINT_Lib_Radau5_dev
