! test-batchF: constant-volume batch reactor with the dedicated routine (explicit), the general procedure
! (general) and, built with Cantera, the Cantera rates (Cantera) for the cases of test/batch/cases.txt.
! Run from test/batch:
!   test-batchF                         asks for the mode, then runs every case
!   test-batchF verification            every case, 1000 steps: <case>/batch-<backend>.dat (batch-verification.py)
!   test-batchF performance             every case, 1 step: the times in comp-batch-<backend>.dat
!   test-batchF check <case>            one case, 1000 steps, against the Cantera reference
!                                       reference/<case>.dat (written once by test-batchCXX --reference):
!                                       final temperature, time of half the temperature rise, mean deviation,
!                                       and general = explicit. Needs no Cantera at run time. Exit code 1 on failure.
! Built with Cantera, the Cantera phase reads element-standard-entropies.yaml (the standard entropies of the
! elements) from the working directory, test/batch.
program test
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use oslo
# if defined (CANTERA)
  use FLINT_cantera_load
# endif
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  use FLINT_Lib_Chemistry_data
  use FLINT_Lib_Chemistry_wdot
  use FLINT_Lib_Chemistry_rhs
  use FLINT_Load_chemistry
  implicit none

  ! acceptance of test-batchF check (FLINT at RT = AT = 1e-7 against Cantera at rtol 1e-10, atol 1e-15;
  ! the worst case of the cases.txt of this commit in brackets)
  real(8), parameter :: tol_Tend = 5d-4   ! |T_end/T_end,ref - 1|                          (8.3e-5, WD)
  real(8), parameter :: tol_ign  = 2d-3   ! |t_half/t_half,ref - 1|, half the rise of T_ref  (9.2e-4, Pelucchi)
  real(8), parameter :: tol_L1   = 1d-3   ! mean |T - T_ref| / |rise of T_ref|              (2.1e-4, Pelucchi)
  real(8), parameter :: tol_gen  = 1d-5   ! max |T_general - T_explicit| / |rise of T_ref|  (1.8e-7, CORIA)
  integer, parameter :: nstep_verification = 1000, maxY = 16

  type :: batch_case
    character(64) :: name = '', yaml = ''
    logical       :: general = .false.
    real(8)       :: tend = 0d0, p = 0d0, T = 0d0
    integer       :: nY = 0
    character(32) :: Ysp(maxY) = ''
    real(8)       :: Yval(maxY) = 0d0
  end type batch_case

  type(batch_case), allocatable :: cases(:)
  real(8), allocatable :: sp_Y(:), Y(:), RT(:), AT(:)
  real(8)              :: timein, timeout, dt
  integer              :: neq, nstep, ncase, c, sim_type, nfail
  integer              :: iopt(3)
  character(32)        :: solver, mech_name
  character(256)       :: mode, only
  logical              :: check, found

# if defined(SUNDIALS)
  solver = 'cvode'
# else
  solver = 'ros4'
# endif
  iopt = 0
  iopt(1) = 1000000

  call get_command_argument(1, mode)
  call get_command_argument(2, only)
  if (mode == '') then
    write(*,*)'what kind of simulation do you want to run?'
    write(*,*)'1) verification'
    write(*,*)'2) performance'
    read(*,*) sim_type
    if (sim_type==1) then
      mode = 'verification'
    elseif (sim_type==2) then
      mode = 'performance'
    else
      stop "choose 1 or 2!"
    endif
  endif
  check = (mode == 'check')
  if (.not.(check .or. mode == 'verification' .or. mode == 'performance') .or. (check .and. only == '')) then
    write(*,'(A)') 'usage: test-batchF [verification | performance | check <case>]'
    stop 2
  endif
  nstep = merge(1, nstep_verification, mode == 'performance')

  call read_cases('cases.txt')

  if (.not. check) then
    open(unit=10, file='comp-batch-general.dat', status='replace', form='formatted')
    open(unit=20, file='comp-batch-explicit.dat', status='replace', form='formatted')
# if defined (CANTERA)
    open(unit=30, file='comp-batch-canteraFor.dat', status='replace', form='formatted')
# endif
  endif

  nfail = 0
  found = .false.
  do c = 1, ncase
    if (check .and. cases(c)%name /= only) cycle
    found = .true.
    call run_case(cases(c))
  enddo

  if (check .and. .not. found) then
    write(*,'(A)') '[FAIL] case '//trim(only)//' is not in cases.txt'
    stop 1
  endif
  if (nfail > 0) then
    write(*,'(A,I0,A)') ' Verdict -> fail (', nfail, ' checks)'
    stop 1
  endif
  if (check) write(*,'(A)') ' Verdict -> pass'

contains

  subroutine run_case(cs)
    type(batch_case), intent(in) :: cs
    character(:), allocatable    :: dir
    real(8), allocatable         :: times(:), T_expl(:), T_gen(:)
    integer                      :: k, i, err

    dir = '../../database/'//trim(cs%name)//'/'
    call execute_command_line('mkdir -p '//trim(cs%name))
    err = read_idealgas_thermo(dir)
    ! ios = 5: no composition.txt, the elements of the species (CEA only)
    if (err /= 0 .and. err /= 5) then; write(*,'(A,I0)') '[FAIL] read_idealgas_thermo '//dir//': ios=', err; stop 1; endif
    err = read_chemistry( folder=dir, mech_name=mech_name )
    if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_chemistry '//dir//': ios=', err; stop 1; endif
# if defined (CANTERA)
    call load_phase(gas, dir//trim(cs%yaml))
# endif

    neq = ns + 1
    allocate(Y(neq), sp_Y(ns), RT(neq), AT(neq))
    allocate(times(nstep), T_expl(nstep), T_gen(nstep))

    sp_Y = 1d-20
    do k = 1, cs%nY
      i = findloc(species_names, cs%Ysp(k), dim=1)
      if (i == 0) then; write(*,'(A)') '[FAIL] '//trim(cs%name)//': no species '//trim(cs%Ysp(k)); stop 1; endif
      sp_Y(i) = cs%Yval(k)
    enddo

    RT = 1d-7
    AT = 1d-7
    call setup_odesolver(N=neq,solver=solver,RT=RT,AT=AT,iopt=iopt)
    dt = cs%tend/nstep

# if defined (CANTERA)
    if (.not. check) call integrate(cs, 'Cantera', 'cantera', .true., 30, times, T_gen)
# endif

    !! Native with coded mechanism
    call Assign_Mechanism(mech_name)
    call integrate(cs, 'explicit', 'explicit', .false., 20, times, T_expl)

    !! Native without coded mechanism
    if (cs%general) then
      call Assign_Mechanism('nemo')
      call integrate(cs, 'general', 'general', .false., 10, times, T_gen)
    endif

    if (check) call compare(cs, times, T_expl, T_gen)

    deallocate(Y, sp_Y, RT, AT, times, T_expl, T_gen)
    deallocate(wm_tab); deallocate(Ri_tab)
    deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
    deallocate(species_names)
    if (allocated(elements_names)) deallocate(elements_names)
    if (allocated(species_composition)) deallocate(species_composition)
    call free_chemistry_data()
  end subroutine run_case

  !> Integrates the case from its initial state, writes <case>/batch-<file>.dat and the time to unit ucomp.
  subroutine integrate(cs, label, file, cantera, ucomp, times, temps)
    type(batch_case), intent(in) :: cs
    character(*), intent(in)     :: label, file
    logical, intent(in)          :: cantera
    integer, intent(in)          :: ucomp
    real(8), intent(out)         :: times(:), temps(:)
    real(8)                      :: time1, time2
    integer                      :: n, u, err

    err = 0

    call initialize(cs)
    open(newunit=u, file=trim(cs%name)//'/batch-'//file//'.dat', status='replace', form='formatted')
    temps = huge(1d0)
    call cpu_time(time1)
    do n = 1, nstep
      timein = timeout; timeout = timeout+dt
      if (cantera) then
# if defined (CANTERA)
        call run_odesolver(neq,timein,timeout,Y,rhs_cantera,no_jacobian,0,err)
# endif
      else
        call run_odesolver(neq,timein,timeout,Y,rhs_native,jac_native,IJAC_chem,err)
      endif
      times(n) = timeout
      temps(n) = y(neq)
      write(u,*) timeout, y(neq)
      if (err < 0) then
        write(*,'(A,I0,A,ES12.5)') '[FAIL] '//trim(cs%name)//' '//label//': solver error ', err, ' at t = ', timeout
        nfail = nfail + 1
        exit
      endif
    enddo
    call cpu_time(time2)
    close(u)

    write(*,*) trim(cs%name)//' '//label//' time =', time2-time1
    if (.not. check) write(ucomp,*) trim(cs%name), time2-time1
  end subroutine integrate

  subroutine initialize(cs)
    type(batch_case), intent(in) :: cs
    real(8) :: R, rho
    R = f_Rtot(sp_Y)
    rho = cs%p/(R*cs%T)
    Y(1:ns) = rho*sp_Y
    Y(neq) = cs%T
    timein  = 0.D0
    timeout = 0.D0
  end subroutine initialize

  !> test-batchF check: the traces against reference/<case>.dat, general against explicit.
  subroutine compare(cs, times, T_expl, T_gen)
    type(batch_case), intent(in) :: cs
    real(8), intent(in)          :: times(:), T_expl(:), T_gen(:)
    real(8), allocatable         :: times_ref(:), T_ref(:)
    character(512)               :: line
    character(:), allocatable    :: file
    real(8)                      :: rise
    integer                      :: u, ios, n

    file = 'reference/'//trim(cs%name)//'.dat'
    allocate(times_ref(nstep), T_ref(nstep))
    open(newunit=u, file=file, status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') '[FAIL] no reference '//file//' (test-batchCXX --reference '//trim(cs%name)//')'
      nfail = nfail + 1
      return
    endif
    n = 0
    do
      read(u,'(A)',iostat=ios) line
      if (ios /= 0) exit
      if (line == '' .or. line(1:1) == '#') cycle
      n = n + 1
      if (n > nstep) exit
      read(line,*) times_ref(n), T_ref(n)
    enddo
    close(u)
    call verdict('reference '//file//': the 1000 steps of test-batchF', n == nstep)
    if (n /= nstep) return
    call verdict('reference '//file//': the same output times', maxval(abs(times/times_ref - 1d0)) < 1d-9)

    rise = T_ref(nstep) - cs%T
    call check_trace('explicit', times_ref, T_ref, cs%T, rise, T_expl)
    if (cs%general) then
      call check_trace('general', times_ref, T_ref, cs%T, rise, T_gen)
      call verdict_value('general = explicit, max |dT| / rise', maxval(abs(T_gen - T_expl))/abs(rise), tol_gen)
    endif
  end subroutine compare

  subroutine check_trace(label, times_ref, T_ref, T0, rise, T)
    character(*), intent(in) :: label
    real(8), intent(in)      :: times_ref(:), T_ref(:), T0, rise, T(:)
    if (.not. all(ieee_is_finite(T)) .or. any(T == huge(1d0))) then
      call verdict(label//': finite temperature at every step', .false.)
      return
    endif
    call verdict_value(label//' vs Cantera: final temperature', abs(T(nstep)/T_ref(nstep) - 1d0), tol_Tend)
    call verdict_value(label//' vs Cantera: time of half the temperature rise', &
      abs(t_half(times_ref, T, T0, rise)/t_half(times_ref, T_ref, T0, rise) - 1d0), tol_ign)
    call verdict_value(label//' vs Cantera: mean |dT| / rise', sum(abs(T - T_ref))/nstep/abs(rise), tol_L1)
  end subroutine check_trace

  !> First time at which T reaches T0 plus half the rise (linear between the output times, from T0 at
  !> time 0); huge if it never does.
  real(8) function t_half(times, T, T0, rise)
    real(8), intent(in) :: times(:), T(:), T0, rise
    real(8) :: Th, time_prev, T_prev
    integer :: i
    Th = T0 + 0.5d0*rise
    time_prev = 0d0
    T_prev = T0
    t_half = huge(1d0)
    do i = 1, size(T)
      if ((T(i) - Th)*sign(1d0, rise) >= 0d0) then
        t_half = time_prev + (Th - T_prev)*(times(i) - time_prev)/(T(i) - T_prev)
        return
      endif
      time_prev = times(i)
      T_prev = T(i)
    enddo
  end function t_half

  subroutine verdict_value(what, value, tol)
    character(*), intent(in) :: what
    real(8), intent(in)      :: value, tol
    character(48)            :: num
    write(num,'(A,ES9.2,A,ES8.1,A)') ' (', value, ' <= ', tol, ')'
    call verdict(what//trim(num), value <= tol)
  end subroutine verdict_value

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

  !> Reads the cases of cases.txt (format in its header).
  subroutine read_cases(file)
    character(*), intent(in) :: file
    character(1024)          :: line
    character(64)            :: tok(6+maxY)
    character(8)             :: general
    integer                  :: u, ios, ntok, k, colon

    allocate(cases(0))
    open(newunit=u, file=file, status='old', action='read', iostat=ios)
    if (ios /= 0) then; write(*,'(A)') '[FAIL] cannot open '//file//' (run from test/batch)'; stop 1; endif
    do
      read(u,'(A)',iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (line == '' .or. line(1:1) == '#') cycle
      call split(line, tok, ntok)
      if (ntok < 7 .or. ntok > 6+maxY) then
        write(*,'(A)') '[FAIL] '//file//': a case is: case yaml general tend p T species:Y ...: '//trim(line); stop 1
      endif
      cases = [cases, batch_case()]
      associate (cs => cases(size(cases)))
        cs%name = tok(1)
        cs%yaml = tok(2)
        general = tok(3)
        cs%general = (general == 'yes')
        read(tok(4),*) cs%tend
        read(tok(5),*) cs%p
        read(tok(6),*) cs%T
        cs%nY = ntok - 6
        do k = 1, cs%nY
          colon = index(tok(6+k), ':')
          cs%Ysp(k) = tok(6+k)(:colon-1)
          read(tok(6+k)(colon+1:),*) cs%Yval(k)
        enddo
      end associate
    enddo
    close(u)
    ncase = size(cases)
  end subroutine read_cases

  !> Splits a line at blanks.
  subroutine split(line, tok, ntok)
    character(*), intent(in)  :: line
    character(*), intent(out) :: tok(:)
    integer, intent(out)      :: ntok
    integer :: i, j
    ntok = 0
    i = 1
    do while (i <= len_trim(line))
      if (line(i:i) == ' ') then; i = i + 1; cycle; endif
      j = i
      do while (j < len_trim(line) .and. line(j+1:j+1) /= ' ')
        j = j + 1
      enddo
      ntok = ntok + 1
      if (ntok <= size(tok)) tok(ntok) = line(i:j)
      i = j + 1
    enddo
  end subroutine split

end program test
