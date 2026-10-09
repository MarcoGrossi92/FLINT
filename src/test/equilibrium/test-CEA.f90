program test_eq
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  use FLINT_CEA_setup
  use FLINT_CEA_solver
  implicit none
  integer :: err
  real(kind=8), dimension(:), allocatable :: rhoi, y_eq, y0
  real(kind=8) :: T_, rho_, teq, blessed_Teq, blessed_y
  character(len=16) :: verdict
  integer :: i, N, nfail = 0
  real(kind=8) :: of, press

  T_ = 1000d0
  rho_ = 3.25d0

  N = 1000

  !-------------------------------------------------------------------------------------------------
  ! WD
  !-------------------------------------------------------------------------------------------------

  write(*,*)
  write(*,*) 'WD'

  call execute_command_line('mkdir -p WD/')
  err = read_idealgas_thermo('../../database/WD/')
  allocate(rhoi(ns))
  allocate(y_eq, mold=rhoi)

  ! Initialize the global variables for the CEA solver, which are needed to solve the equilibrium problem
  call CEA_initialize_global()

  write(*,*) 'Testing the equilibrium solver for different mixture ratios...'
  open(unit=313, file='WD/FLINT-CEA.txt', status='replace')
  do i = 1, N
    of = 0.01d0 * (100d0/0.01d0)**(real(i-1, kind=8)/real(N-1, kind=8))
    rhoi = 1d-20
    rhoi(2) = of/(of+1d0)   ! O2: database/WD is in the order of WD.f90 (CH4, O2, CO2, H2O, CO)
    rhoi(1) = 1d0/(of+1d0)
    rhoi = rhoi * rho_
    call CEA_solve(T_, rhoi, teq, y_eq)
    if (mod(i,100)==0) write(*,'(A,F10.3,A,F10.3,A)') ' Mixture ratio = ', of, ' -> equilibrium temperature = ', teq, ' K'
    write(313,*) of, teq
  end do
  close(313)
  call check_sweep('WD')

  write(*,*) 'Testing the equilibrium solver for a single point...'
  blessed_Teq = 4919.064001247616 
  blessed_y = 3.26511387e-01

  rhoi = 1d-20
  rhoi(2) = 0.8d0   ! O2
  rhoi(1) = 0.2d0

  call CEA_solve(T_, rhoi, teq, y_eq)

  write(*,*)'FLINT equilibrium temperature   = ', teq
  write(*,*)'Cantera equilibrium temperature = ', blessed_Teq
  write(*,*)'FLINT yCO   = ', y_eq(5)
  write(*,*)'Cantera yCO = ', blessed_y

  verdict = 'success'
  if (isnan(teq)) verdict = 'fail'
  if (abs(teq-blessed_Teq)/blessed_Teq*100d0>1d0) verdict = 'fail'
  if (abs(y_eq(5)-blessed_y)/blessed_y*100d0>1d0) verdict = 'fail'

  write(*,'(2A20)') 'Verdict -> ', verdict
  if (verdict == 'fail') nfail = nfail + 1

  deallocate(rhoi); deallocate(y_eq)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(species_names); deallocate(elements_names)
  deallocate(species_composition)

  !-------------------------------------------------------------------------------------------------
  ! ZK
  !-------------------------------------------------------------------------------------------------

  write(*,*)
  write(*,*) 'ZK'

  call execute_command_line('mkdir -p ZK/')
  err = read_idealgas_thermo('../../database/ZK/')
  allocate(rhoi(ns))
  allocate(y_eq, mold=rhoi)

  ! Initialize the global variables for the CEA solver, which are needed to solve the equilibrium problem
  call CEA_initialize_global()

  write(*,*) 'Testing the equilibrium solver for different mixture ratios...'
  open(unit=313, file='ZK/FLINT-CEA.txt', status='replace')
  do i = 1, N
    of = 0.01d0 * (100d0/0.01d0)**(real(i-1, kind=8)/real(N-1, kind=8))
    rhoi = 1d-20
    rhoi(6) = of/(of+1d0)
    rhoi(17) = 1d0/(of+1d0)
    rhoi = rhoi * rho_
    call CEA_solve(T_, rhoi, teq, y_eq)
    if (mod(i,100)==0) write(*,'(A,F10.3,A,F10.3,A)') ' Mixture ratio = ', of, ' -> equilibrium temperature = ', teq, ' K'
    write(313,*) of, teq
  end do
  close(313)
  call check_sweep('ZK')

  write(*,*) 'Testing the equilibrium solver for a single point...'
  blessed_Teq = 3611.151626300349
  blessed_y = 1.05464743e-01

  rhoi = 1d-20
  rhoi(6) = 0.8d0
  rhoi(17) = 0.2d0

  call CEA_solve(T_, rhoi, teq, y_eq)

  write(*,*)'FLINT equilibrium temperature   = ', teq
  write(*,*)'Cantera equilibrium temperature = ', blessed_Teq
  write(*,*)'FLINT yOH   = ', y_eq(7)
  write(*,*)'Cantera yOH = ', blessed_y

  verdict = 'success'
  if (isnan(teq)) verdict = 'fail'
  if (abs(teq-blessed_Teq)/blessed_Teq*100d0>1d0) verdict = 'fail'
  if (abs(y_eq(7)-blessed_y)/blessed_y*100d0>1d0) verdict = 'fail'

  write(*,'(2A20)') 'Verdict -> ', verdict
  if (verdict == 'fail') nfail = nfail + 1

  deallocate(rhoi); deallocate(y_eq)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(species_names); deallocate(elements_names)
  deallocate(species_composition)

  !-------------------------------------------------------------------------------------------------
  ! TSR-GP-24
  !-------------------------------------------------------------------------------------------------

  write(*,*)
  write(*,*) 'TSR-GP-24'

  call execute_command_line('mkdir -p TSR-GP-24/')
  err = read_idealgas_thermo('../../database/TSR-GP-24/')
  allocate(rhoi(ns))
  allocate(y_eq, mold=rhoi)

  ! Initialize the global variables for the CEA solver, which are needed to solve the equilibrium problem
  call CEA_initialize_global()

  write(*,*) 'Testing the equilibrium solver for different mixture ratios...'
  open(unit=313, file='TSR-GP-24/FLINT-CEA.txt', status='replace')
  do i = 1, N
    of = 0.01d0 * (100d0/0.01d0)**(real(i-1, kind=8)/real(N-1, kind=8))
    rhoi = 1d-20
    rhoi(19) = of/(of+1d0)
    rhoi(22) = 1d0/(of+1d0)
    rhoi = rhoi * rho_
    call CEA_solve(T_, rhoi, teq, y_eq)
    if (mod(i,100)==0) write(*,'(A,F10.3,A,F10.3,A)') ' Mixture ratio = ', of, ' -> equilibrium temperature = ', teq, ' K'
    write(313,*) of, teq
  end do
  close(313)
  call check_sweep('TSR-GP-24')

  write(*,*) 'Testing the equilibrium solver for a single point...'
  blessed_Teq = 3616.638618717671
  blessed_y = 1.00013716e-01

  rhoi = 1d-20
  rhoi(19) = 0.8d0
  rhoi(22) = 0.2d0

  call CEA_solve(T_, rhoi, teq, y_eq)

  write(*,*)'FLINT equilibrium temperature   = ', teq
  write(*,*)'Cantera equilibrium temperature = ', blessed_Teq
  write(*,*)'FLINT yOH   = ', y_eq(4)
  write(*,*)'Cantera yOH = ', blessed_y

  verdict = 'success'
  if (isnan(teq)) verdict = 'fail'
  if (abs(teq-blessed_Teq)/blessed_Teq*100d0>1d0) verdict = 'fail'
  if (abs(y_eq(4)-blessed_y)/blessed_y*100d0>1d0) verdict = 'fail'

  write(*,'(2A20)') 'Verdict -> ', verdict
  if (verdict == 'fail') nfail = nfail + 1

  deallocate(rhoi); deallocate(y_eq)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(species_names); deallocate(elements_names)
  deallocate(species_composition)

  !-------------------------------------------------------------------------------------------------
  ! Ecker
  !-------------------------------------------------------------------------------------------------

  write(*,*)
  write(*,*) 'ECKER'

  call execute_command_line('mkdir -p Ecker/')
  err = read_idealgas_thermo('../../database/Ecker/')
  allocate(rhoi(ns))
  allocate(y_eq, mold=rhoi)

  rhoi = 1d-20
  rhoi(7)  = 0.5d0
  rhoi(12) = 0.2d0
  rhoi(1)  = 0.2d0
  rhoi(2)  = 0.1d0
  T_ = 3000d0

  ! Initialize the global variables for the CEA solver, which are needed to solve the equilibrium problem
  call CEA_initialize_global()

  write(*,*) 'Testing the equilibrium solver for different mixture ratios...'
  open(unit=313, file='Ecker/FLINT-CEA.txt', status='replace')
  y0 = rhoi   ! mass fractions; the partial densities of each pressure are y0*rho_
  do i = 1, N
    press = 1d5 * 1d-5 * 10d0**(real(i-1, kind=8)*log10(100d0/1d-5)/real(N-1, kind=8))   ! 1e-5 .. 100 bar
    rho_ = press/(f_rtot(y0)*T_)
    call CEA_solve(T_, y0*rho_, teq, y_eq)
    if (mod(i,100)==0) write(*,'(A,F10.3,A,F10.3,A)') ' Pressure = ', press*1e-5, ' -> equilibrium temperature = ', teq, ' K'
    write(313,*) press, teq
  end do
  close(313)
  call check_sweep('Ecker')

  write(*,*) 'Testing the equilibrium solver for a single point...'
  blessed_Teq = 1571.8416518627969
  blessed_y = 1.93015081e-01

  T_ = 1000d0
  rhoi = 1d-20
  rhoi(12) = 0.18798856d0
  rhoi(1) = 0.00534534d0
  rhoi(14) = 0.8066661d0

  call CEA_solve(T_, rhoi, teq, y_eq)

  write(*,*)'FLINT equilibrium temperature   = ', teq
  write(*,*)'Cantera equilibrium temperature = ', blessed_Teq
  write(*,*)'FLINT yHCL   = ', y_eq(11)
  write(*,*)'Cantera yHCL = ', blessed_y

  verdict = 'success'
  if (isnan(teq)) verdict = 'fail'
  if (abs(teq-blessed_Teq)/blessed_Teq*100d0>1d0) verdict = 'fail'
  if (abs(y_eq(11)-blessed_y)/blessed_y*100d0>1d0) verdict = 'fail'

  write(*,'(2A20)') 'Verdict -> ', verdict
  if (verdict == 'fail') nfail = nfail + 1

  deallocate(rhoi); deallocate(y_eq)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(species_names); deallocate(elements_names)
  deallocate(species_composition)

  ! exit code 1 when a case fails (CTest)
  if (nfail > 0) stop 1

contains

  !> The sweep just written, <mech>/FLINT-CEA.txt, against the Cantera reference reference/<mech>.dat
  !> (test-equilCXX): the same swept values, equilibrium temperatures within sweep_tol (relative).
  subroutine check_sweep(mech)
    character(*), intent(in) :: mech
    real(kind=8), parameter :: sweep_tol = 2d-3   ! 8.4e-4 for TSR-GP-24 with the references of this commit
    real(kind=8) :: x, t, xr, tr, emax
    character(len=256) :: line
    integer :: uf, ur, ios, n
    logical :: same_x
    open(newunit=ur, file='reference/'//mech//'.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A)') ' [FAIL] no reference reference/'//mech//'.dat (test-equilCXX '//mech//')'
      nfail = nfail + 1
      return
    endif
    open(newunit=uf, file=mech//'/FLINT-CEA.txt', status='old', action='read')
    n = 0; emax = 0d0; same_x = .true.
    do
      read(ur,'(A)',iostat=ios) line
      if (ios /= 0) exit
      if (line(1:1) == '#') cycle
      read(line,*) xr, tr
      read(uf,*,iostat=ios) x, t
      if (ios /= 0) exit
      n = n + 1
      same_x = same_x .and. abs(x/xr - 1d0) < 1d-9
      emax = max(emax, abs(t/tr - 1d0))
      if (isnan(t)) emax = huge(1d0)
    enddo
    close(ur); close(uf)
    if (n == N .and. same_x .and. emax <= sweep_tol) then
      write(*,'(A,I0,A,ES9.2,A,ES8.1,A)') ' [ok]   '//mech//' sweep = Cantera at ', n, ' states (', emax, ' <= ', sweep_tol, ')'
    else
      write(*,'(A,I0,A,L1,A,ES9.2,A)') ' [FAIL] '//mech//' sweep: ', n, ' states, same swept values ', same_x, &
        ', max |T/T_ref - 1| = ', emax, ''
      nfail = nfail + 1
    endif
  end subroutine check_sweep

end program test_eq
