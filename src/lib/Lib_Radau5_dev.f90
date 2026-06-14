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
  implicit none
  private
  public :: radau5_dev, flint_acc_upload_tables

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

contains

  subroutine set_radau5_dev_tols(n, rt, at)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: rt(n), at(n)
    rtol_dev(1:n) = rt(1:n)
    atol_dev(1:n) = at(1:n)
#   if defined (_OPENACC)
    !$acc update device(rtol_dev, atol_dev)
#   endif
  end subroutine set_radau5_dev_tols

  !----------------------------------------------------------------------------
  ! Copy the read-only chemistry/thermo tables to the device. Host-only;
  ! call once after read_idealgas_thermo + read_chemistry. The matching
  ! "!$acc declare create" directives live in the table modules.
  !----------------------------------------------------------------------------
  subroutine flint_acc_upload_tables()
    use FLINT_Lib_Thermodynamic, only: ns, Tmin, Tmax, wm_tab, Ri_tab, &
                                       cp_tab, h_tab, mi_tab, k_tab, &
                                       Mi_Mj_pow_m025, inv_sqrt8_1p
    use FLINT_Lib_Chemistry_data, only: kf_tab, kb_tab
    implicit none
#   if defined (_OPENACC)
    !$acc update device(ns, Tmin, Tmax)
    !$acc update device(wm_tab, Ri_tab, cp_tab, h_tab)
    !$acc update device(kf_tab, kb_tab)
    ! Transport tables for the device thermo-cache kernel (VISC cases only;
    ! allocated by read_idealgas_transport before this routine runs).
    if (allocated(mi_tab)) then
      !$acc update device(mi_tab, k_tab, Mi_Mj_pow_m025, inv_sqrt8_1p)
    end if
#   endif
  end subroutine flint_acc_upload_tables


  !----------------------------------------------------------------------------
  ! Integrate Y over [X, XEND]. RTOL/ATOL are NOT modified (the internal
  ! Hairer transform works on copies). Statistics returned per call.
  !----------------------------------------------------------------------------
  subroutine radau5_dev(n, x_in, y, xend, rtol_in, atol_in, idid, &
                        nfcn, njac, nstep, naccpt, nrejct)
    !$acc routine seq
    implicit none
    integer, intent(in)    :: n
    real(8), intent(in)    :: x_in, xend
    real(8), intent(inout) :: y(n)
    real(8), intent(in)    :: rtol_in(n), atol_in(n)
    integer, intent(out)   :: idid
    integer, intent(out)   :: nfcn, njac, nstep, naccpt, nrejct

    ! ---- fixed per-thread work arrays (original: partitions of WORK/IWORK)
    real(8) :: z1(NSMX), z2(NSMX), z3(NSMX), y0(NSMX), scal(NSMX)
    real(8) :: f1(NSMX), f2(NSMX), f3(NSMX)
    real(8) :: fjac(NSMX,NSMX), e1(NSMX,NSMX), e2r(NSMX,NSMX), e2i(NSMX,NSMX)
    real(8) :: cont(4*NSMX)
    integer :: ip1(NSMX), ip2(NSMX)
    real(8) :: rtol(NSMX), atol(NSMX)
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
    expm = 2.0d0/3.0d0
    do i = 1, n
      quot    = atol_in(i)/rtol_in(i)
      rtol(i) = 0.1d0*rtol_in(i)**expm
      atol(i) = rtol(i)*quot
    end do
    ! --- fnewt (computed from the TRANSFORMED rtol(1), as in the original)
    tolst = rtol(1)
    fnewt = max(10*uround/tolst, min(0.03d0, tolst**0.5d0))

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
    call rhs_native(n, x, y, y0)
    nfcn = nfcn+1

    ! --- basic integration step
 10 continue
    ! *** computation of the Jacobian (IJAC=1: analytically)
    njac = njac+1
    call jac_native(n, x, y, fjac, NSMX, rpar_d, ipar_d)
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
    call dec_dev(n, NSMX, e1, ip1, ier)
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
    call decc_dev(n, NSMX, e2r, e2i, ip2, ier)
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
    call rhs_native(n, x+c1*h, cont, z1)
    do i = 1, n
      cont(i) = y(i)+z2(i)
    end do
    call rhs_native(n, x+c2*h, cont, z2)
    do i = 1, n
      cont(i) = y(i)+z3(i)
    end do
    call rhs_native(n, xph, cont, z3)
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
    call sol_dev (n, NSMX, e1, z1, ip1)
    call solc_dev(n, NSMX, e2r, e2i, z2, z3, ip2)
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
    call sol_dev(n, NSMX, e1, cont, ip1)
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
        call rhs_native(n, x, cont, f1)
        nfcn = nfcn+1
        do i = 1, n
          cont(i) = f1(i)+f2(i)
        end do
        call sol_dev(n, NSMX, e1, cont, ip1)
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
      call rhs_native(n, x, y, y0)
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

end module FLINT_Lib_Radau5_dev
