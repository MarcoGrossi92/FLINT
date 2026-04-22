program test_eq
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  use FLINT_CEA_setup
  use FLINT_CEA_solver
  implicit none
  integer :: err
  real(kind=8), dimension(:), allocatable :: rhoi, y_eq
  real(kind=8) :: T_, rho_, teq, blessed_Teq, blessed_y
  character(len=16) :: verdict
  integer :: i, N
  real(kind=8) :: of, press

  T_ = 1000d0
  rho_ = 3.25d0

  N = 1000

  !-------------------------------------------------------------------------------------------------
  ! WD
  !-------------------------------------------------------------------------------------------------

  write(*,*)
  write(*,*) 'WD'

  err = read_idealgas_thermo('WD/INPUT')
  allocate(rhoi(ns))
  allocate(y_eq, mold=rhoi)

  ! Initialize the global variables for the CEA solver, which are needed to solve the equilibrium problem
  call CEA_initialize_global()

  write(*,*) 'Testing the equilibrium solver for different mixture ratios...'
  open(unit=313, file='WD/OUTPUT/FLINT-CEA.txt', status='replace')
  do i = 1, N
    of = 0.01d0 * (100d0/0.01d0)**(real(i-1, kind=8)/real(N-1, kind=8))
    rhoi = 1d-20
    rhoi(5) = of/(of+1d0)
    rhoi(1) = 1d0/(of+1d0)
    rhoi = rhoi * rho_
    call CEA_solve(T_, rhoi, teq, y_eq)
    if (mod(i,100)==0) write(*,'(A,F10.3,A,F10.3,A)') ' Mixture ratio = ', of, ' -> equilibrium temperature = ', teq, ' K'
    write(313,*) of, teq
  end do
  close(313)

  write(*,*) 'Testing the equilibrium solver for a single point...'
  blessed_Teq = 4919.064001247616 
  blessed_y = 3.26511387e-01

  rhoi = 1d-20
  rhoi(5) = 0.8d0
  rhoi(1) = 0.2d0

  call CEA_solve(T_, rhoi, teq, y_eq)

  write(*,*)'FLINT equilibrium temperature   = ', teq
  write(*,*)'Cantera equilibrium temperature = ', blessed_Teq
  write(*,*)'FLINT yCO   = ', y_eq(2)
  write(*,*)'Cantera yCO = ', blessed_y

  verdict = 'success'
  if (isnan(teq)) verdict = 'fail'
  if (abs(teq-blessed_Teq)/blessed_Teq*100d0>1d0) verdict = 'fail'
  if (abs(y_eq(2)-blessed_y)/blessed_Teq*100d0>1d0) verdict = 'fail'

  write(*,'(2A20)') 'Verdict -> ', verdict

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

  err = read_idealgas_thermo('ZK/INPUT')
  allocate(rhoi(ns))
  allocate(y_eq, mold=rhoi)

  ! Initialize the global variables for the CEA solver, which are needed to solve the equilibrium problem
  call CEA_initialize_global()

  write(*,*) 'Testing the equilibrium solver for different mixture ratios...'
  open(unit=313, file='ZK/OUTPUT/FLINT-CEA.txt', status='replace')
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

  write(*,*) 'Testing the equilibrium solver for a single point...'
  blessed_Teq = 3611.151626300349
  blessed_y = 1.05464743e-01

  rhoi = 1d-20
  rhoi(6) = 0.8d0
  rhoi(17) = 0.2d0

  call CEA_solve(T_, rhoi, teq, y_eq)

  write(*,*)'FLINT equilibrium temperature   = ', teq
  write(*,*)'Cantera equilibrium temperature = ', blessed_Teq
  write(*,*)'FLINT yCO   = ', y_eq(7)
  write(*,*)'Cantera yCO = ', blessed_y

  verdict = 'success'
  if (isnan(teq)) verdict = 'fail'
  if (abs(teq-blessed_Teq)/blessed_Teq*100d0>1d0) verdict = 'fail'
  if (abs(y_eq(7)-blessed_y)/blessed_Teq*100d0>1d0) verdict = 'fail'

  write(*,'(2A20)') 'Verdict -> ', verdict

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

  err = read_idealgas_thermo('TSR-GP-24/INPUT')
  allocate(rhoi(ns))
  allocate(y_eq, mold=rhoi)

  ! Initialize the global variables for the CEA solver, which are needed to solve the equilibrium problem
  call CEA_initialize_global()

  ! write(*,*) 'Testing the equilibrium solver for different mixture ratios...'
  ! open(unit=313, file='TSR-GP-24/OUTPUT/FLINT-CEA.txt', status='replace')
  ! do i = 1, N
  !   of = 0.01d0 * (100d0/0.01d0)**(real(i-1, kind=8)/real(N-1, kind=8))
  !   rhoi = 1d-20
  !   rhoi(19) = of/(of+1d0)
  !   rhoi(22) = 1d0/(of+1d0)
  !   rhoi = rhoi * rho_
  !   call CEA_solve(T_, rhoi, teq, y_eq)
  !   if (mod(i,100)==0) write(*,'(A,F10.3,A,F10.3,A)') ' Mixture ratio = ', of, ' -> equilibrium temperature = ', teq, ' K'
  !   write(313,*) of, teq
  ! end do
  ! close(313)

  write(*,*) 'Testing the equilibrium solver for a single point...'
  blessed_Teq = 3616.638618717671
  blessed_y = 1.00013716e-01

  rhoi = 1d-20
  rhoi(19) = 0.8d0
  rhoi(22) = 0.2d0

  call CEA_solve(T_, rhoi, teq, y_eq)

  write(*,*)'FLINT equilibrium temperature   = ', teq
  write(*,*)'Cantera equilibrium temperature = ', blessed_Teq
  write(*,*)'FLINT yCO   = ', y_eq(4)
  write(*,*)'Cantera yCO = ', blessed_y

  verdict = 'success'
  if (isnan(teq)) verdict = 'fail'
  if (abs(teq-blessed_Teq)/blessed_Teq*100d0>1d0) verdict = 'fail'
  if (abs(y_eq(4)-blessed_y)/blessed_Teq*100d0>1d0) verdict = 'fail'

  write(*,'(2A20)') 'Verdict -> ', verdict

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

  err = read_idealgas_thermo('Ecker/INPUT')
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

  ! write(*,*) 'Testing the equilibrium solver for different mixture ratios...'
  ! open(unit=313, file='Ecker/OUTPUT/FLINT-CEA.txt', status='replace')
  ! do i = 1, N
  !   press = 1d5 * 10d0**(real(i-1, kind=8)*log10(100d0/1d-5)/real(N-1, kind=8))
  !   rho_ = press/(f_rtot(rhoi)*T_)
  !   rhoi = rhoi * rho_
  !   call CEA_solve(T_, rhoi, teq, y_eq)
  !   if (mod(i,100)==0) write(*,'(A,F10.3,A,F10.3,A)') ' Pressure = ', press*1e-5, ' -> equilibrium temperature = ', teq, ' K'
  !   write(313,*) press, teq
  ! end do
  ! close(313)

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
  write(*,*)'FLINT yOH   = ', y_eq(11)
  write(*,*)'Cantera yOH = ', blessed_y

  verdict = 'success'
  if (isnan(teq)) verdict = 'fail'
  if (abs(teq-blessed_Teq)/blessed_Teq*100d0>1d0) verdict = 'fail'
  if (abs(y_eq(11)-blessed_y)/blessed_Teq*100d0>1d0) verdict = 'fail'

  write(*,'(2A20)') 'Verdict -> ', verdict

  deallocate(rhoi); deallocate(y_eq)
  deallocate(wm_tab); deallocate(Ri_tab)
  deallocate(h_tab); deallocate(cp_tab); deallocate(dcpi_tab); deallocate(s_tab)
  deallocate(species_names); deallocate(elements_names)
  deallocate(species_composition)

end program test_eq
