! test-orders-warn: the WARNING of the general procedure for a chemistry-info.txt without the 'Reaction orders'
! block (a file written by an older table writer), on copies of the fixture test/orders/JLR-frassoldati made at
! run time (orders/noblock: no block; orders/zeroblock: a block with zero rows) and on a copy of database/WD
! without its block made at run time (orders/wd-noblock):
!  1. a block with zero rows: no WARNING, and the general procedure gives the same omegadot, bit for bit, as the
!     file without the block (the exponents of a block with no rows are the reactant coefficients);
!  2. no block, general procedure selected after the tables were loaded: WARNING (from Assign_Mechanism);
!  3. no block, general procedure selected before the tables were loaded: WARNING (from read_chemistry);
!  4. a hooked name on a folder without the block (orders/wd-noblock selected as WD): no WARNING;
!  5. in a child process (this program with the argument child-noblock): the WARNING is printed once on standard
!     output and once on the error unit for one load, also after several chemistry calls and a second
!     Assign_Mechanism.
! Needs no Cantera. Exit code 1 on failure.
program test
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  use FLINT_Load_chemistry
  use FLINT_Lib_Chemistry_data
  use FLINT_Lib_Chemistry_wdot
  implicit none
  character(32) :: mech_name
  character(len=512) :: self, arg
  integer :: err, nfail, k
  real(8), allocatable :: roi(:), w0(:), w1(:)
  real(8), parameter :: T = 1500.5d0
  real(8), parameter :: wm_ref(9) = [31.998d0, 16.043d0, 18.015d0, 28.010d0, 44.009d0, 2.016d0, 1.008d0, 15.999d0, 17.007d0]
  character(len=s_str_len), parameter :: names_ref(9) = [character(len=s_str_len) :: &
    'O2', 'CH4', 'H2O', 'CO', 'CO2', 'H2', 'H', 'O', 'OH']

  nfail = 0
  call get_command_argument(0, self)
  call get_command_argument(1, arg)
  ! species of phase.txt of the fixture (no thermo tables are needed: the general procedure is called directly)
  ns = 9
  allocate(wm_tab(ns), Ri_tab(ns), species_names(ns))
  wm_tab = wm_ref; Ri_tab = Runiv/wm_tab; species_names = names_ref
  allocate(roi(ns), w0(ns), w1(ns))
  roi = [0.20d0, 0.05d0, 0.02d0, 0.01d0, 0.03d0, 0.001d0, 1d-4, 2d-4, 3d-4]   ! kg/m3

  if (arg == 'child-noblock') then
    ! one load of a file without the block, general procedure selected first: one WARNING on each channel
    call Assign_Mechanism('JLR-Frassoldati')
    err = read_chemistry(folder='orders/noblock', mech_name=mech_name)
    if (err /= 0) then; write(*,'(A,I0)') '[FAIL] child-noblock: read_chemistry ios=', err; stop 2; endif
    do k = 1, 3
      call chemistry_source(roi, T, w0)
    enddo
    call Assign_Mechanism(mech_name)
    call chemistry_source(roi, T, w0)
    write(*,'(A)') 'child-noblock: done'
    stop
  endif

  call execute_command_line('mkdir -p orders/noblock orders/zeroblock && ' // &
    'cp orders/JLR-frassoldati/chemistry-Arrhenius.dat orders/noblock/ && ' // &
    'cp orders/JLR-frassoldati/chemistry-info-noblock.txt orders/noblock/chemistry-info.txt && ' // &
    'cp orders/JLR-frassoldati/chemistry-Arrhenius.dat orders/zeroblock/ && ' // &
    '{ cat orders/JLR-frassoldati/chemistry-info-noblock.txt; printf "\nReaction orders\n0\n"; } ' // &
    '> orders/zeroblock/chemistry-info.txt', exitstat=err)
  if (err /= 0) then; write(*,'(A)') '[FAIL] could not make orders/noblock and orders/zeroblock'; stop 1; endif

  ! 1. a block with zero rows, general procedure selected after the tables
  err = read_chemistry(folder='orders/zeroblock', mech_name=mech_name)
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_chemistry orders/zeroblock: ios=', err; stop 1; endif
  call Assign_Mechanism(mech_name)   ! JLR-Frassoldati is not hooked: general
  call verdict('1a block with zero rows: read as a block (have_orders)', have_orders)
  call verdict('1b block with zero rows: no WARNING', .not. orders_block_warned)
  call chemistry_source(roi, T, w1)
  call free_chemistry_data()

  ! 3. no block, general procedure already selected when the tables are loaded (the WARNING comes from read_chemistry)
  err = read_chemistry(folder='orders/noblock', mech_name=mech_name)
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_chemistry orders/noblock: ios=', err; stop 1; endif
  call verdict('3 no block, general selected before the tables: WARNING', orders_block_warned .and. .not. have_orders)
  call chemistry_source(roi, T, w0)
  call verdict('1c block with zero rows and no block: the same omegadot, bit for bit', &
    all(w0 == w1) .and. any(w0 /= 0d0))
  call free_chemistry_data()

  ! 2. no block, general procedure selected after the tables (the WARNING comes from Assign_Mechanism)
  general_selected = .false.         ! as at the start of a program: nothing selected yet
  err = read_chemistry(folder='orders/noblock', mech_name=mech_name)
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_chemistry orders/noblock: ios=', err; stop 1; endif
  call verdict('2a no block, nothing selected yet: no WARNING at the load', .not. orders_block_warned)
  call Assign_Mechanism(mech_name)
  call verdict('2b no block, general selected after the tables: WARNING', orders_block_warned)
  call free_chemistry_data()

  ! 4. a hooked name on a folder without the block: the compiled routine does not read the block, no WARNING
  !    (database/WD ends with the block, as a table writer writes it: the folder without it is a copy made here)
  call execute_command_line('mkdir -p orders/wd-noblock && ' // &
    'cp ../database/WD/*.txt ../database/WD/*.dat orders/wd-noblock/ && ' // &
    "sed '/^Reaction orders/,$d' ../database/WD/chemistry-info.txt > orders/wd-noblock/chemistry-info.txt", exitstat=err)
  if (err /= 0) then; write(*,'(A)') '[FAIL] could not make orders/wd-noblock'; stop 1; endif
  general_selected = .false.
  deallocate(species_names, wm_tab, Ri_tab)
  err = read_idealgas_thermo('orders/wd-noblock/')
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_idealgas_thermo orders/wd-noblock: ios=', err; stop 1; endif
  err = read_chemistry(folder='orders/wd-noblock/', mech_name=mech_name)
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_chemistry orders/wd-noblock: ios=', err; stop 1; endif
  call Assign_Mechanism(mech_name)
  call verdict('4 hooked name WD on a folder without the block: no WARNING', &
    trim(mech_name) == 'WD' .and. .not. have_orders .and. .not. orders_block_warned)
  call free_chemistry_data()

  ! 5. channels and count, in a child process
  call execute_command_line(trim(self)//' child-noblock > orders/child-noblock.out 2> orders/child-noblock.err', exitstat=err)
  call verdict('5a child with a file without the block: exit code 0', err == 0)
  call verdict('5b the WARNING once on standard output', count_lines('orders/child-noblock.out') == 1)
  call verdict('5c the WARNING once on the error unit', count_lines('orders/child-noblock.err') == 1)

  if (nfail > 0) then
    write(*,'(A,I0,A)') ' Verdict -> fail (', nfail, ' checks)'
    stop 1
  endif
  write(*,'(A)') ' Verdict -> pass'
contains
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
  !> number of lines of a file that carry the WARNING about the missing block
  integer function count_lines(path) result(n)
    character(*), intent(in) :: path
    character(len=1024) :: line
    integer :: u, ios
    n = 0
    open(newunit=u, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) return
    do
      read(u,'(A)',iostat=ios) line
      if (ios /= 0) exit
      if (index(line, "no 'Reaction orders' block") > 0) n = n + 1
    enddo
    close(u)
  end function count_lines
end program test
