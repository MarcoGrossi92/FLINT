MODULE FLINT_CEA_data

!***********************************************************************
! FUNDAMENTAL CONSTANTS FROM:  COHEN,E.RICHARD & TAYLOR,BARRY N.,
! THE 1986 CODATA RECOMMENDED VALUES OF THE FUNDAMENTAL PHYSICAL
! CONSTANTS, J.PHYS.CHEM.REF.DATA, VOL.17, NO.4, 1988, PP 1795-1803.
!***********************************************************************

  IMPLICIT NONE

  ! Define module constants using PARAMETER
  INTEGER, PARAMETER :: MAXNGC = 200
  INTEGER, PARAMETER :: MAXNC = 0
  INTEGER, PARAMETER :: NCOL = 8
  INTEGER, PARAMETER :: MAXMAT = 50
  INTEGER, PARAMETER :: MAXTR = 30
  INTEGER, PARAMETER :: MAXR = 24
  INTEGER, PARAMETER :: MAXEL = 20

  REAL(8) :: Deln(MAXNGC), Sln(MAXNGC)

  INTEGER :: Ip, Iplt, It, Nc, Ng, Ngp1, Nlm, Nomit, Nonly, Np, Npr, Npt
  INTEGER :: Ngc, Nspx, Nt
  INTEGER :: Jcond(45), Jx(MAXEL), Ifz(MAXNC)

  REAL(8) :: Bcheck

  INTEGER :: Imat, Iq1, Isv, Jliq, Jsol, Lsave, Msing

  LOGICAL :: Hp, Ions, Massf, Pderiv, Short, Tp, Vol

  REAL(8) :: R, Size, Trace
  REAL(8), PARAMETER :: Rr=8314.51D0,Pi=3.14159265D0
  REAL(8) :: Atwt(MAXEL), X(MAXMAT)
  REAL(8) :: A(MAXEL, MAXNGC)

  CHARACTER(20) :: Elmt(MAXEL)
  CHARACTER(15) :: Prod(0:MAXNGC)

  REAL(8) :: Cpr(NCOL), Dlvpt(NCOL), Dlvtp(NCOL), Hsum(NCOL), &
             Ppp(NCOL), Ssum(NCOL), Totn(NCOL), Ttt(NCOL), Wm(NCOL)

  INTEGER :: Nreac

  REAL(8) :: Cpsum
  REAL(8) :: Temp(2, MAXNC)
  REAL(8) :: Mw(MAXNGC)

  INTEGER :: Iopt, Npp

  INTEGER :: Nm, Nr, Ntape
  INTEGER :: Ind(MAXTR), Jcm(MAXEL)

  ! Saved pristine copies of arrays that CEA_solve mutates during
  ! component switching. Restored at the start of every CEA_solve call.
  REAL(8) :: A_saved(MAXEL, MAXNGC)
  REAL(8) :: Atwt_saved(MAXEL)
  CHARACTER(20) :: Elmt_saved(MAXEL)
  INTEGER :: Jx_saved(MAXEL), Jcm_saved(MAXEL)
  INTEGER :: Nspx_saved, Nlm_saved

  ! Flag for lazy per-thread initialization (OpenMP)
  LOGICAL :: cea_thread_initialized = .false.

  ! --- OpenMP thread safety ---
  ! Every mutable variable is threadprivate so that concurrent
  ! CEA_solve calls from different threads do not race.
  ! Each thread initializes its own copy on first CEA_solve call.
  !$omp threadprivate(Deln, Sln)
  !$omp threadprivate(Ip, Iplt, It, Nc, Ng, Ngp1, Nlm)
  !$omp threadprivate(Nomit, Nonly, Np, Npr, Npt)
  !$omp threadprivate(Ngc, Nspx, Nt)
  !$omp threadprivate(Jcond, Jx, Ifz)
  !$omp threadprivate(Bcheck)
  !$omp threadprivate(Imat, Iq1, Isv, Jliq, Jsol, Lsave, Msing)
  !$omp threadprivate(Hp, Ions, Massf, Pderiv, Short, Tp, Vol)
  !$omp threadprivate(R, Size, Trace)
  !$omp threadprivate(Atwt, X, A)
  !$omp threadprivate(Elmt, Prod)
  !$omp threadprivate(Cpr, Dlvpt, Dlvtp, Hsum)
  !$omp threadprivate(Ppp, Ssum, Totn, Ttt, Wm)
  !$omp threadprivate(Nreac, Cpsum, Temp, Mw)
  !$omp threadprivate(Iopt, Npp)
  !$omp threadprivate(Nm, Nr, Ntape, Ind, Jcm)
  !$omp threadprivate(A_saved, Atwt_saved, Elmt_saved)
  !$omp threadprivate(Jx_saved, Jcm_saved, Nspx_saved, Nlm_saved)
  !$omp threadprivate(cea_thread_initialized)

END MODULE FLINT_CEA_data
