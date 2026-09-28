! test-contract: the mechanism contract check of Assign_Mechanism.
! Loads database/WD (species in the routine order of WD.f90) and checks that
!  1. the loaded data pass the check and Assign_Mechanism hooks WD;
!  2. the alphabetical order of the old fixture (slots 2 and 5 swapped: CO where the routine reads O2) is refused;
!  3. a different number of Arrhenius reactions is refused;
!  4. species appended after the routine slots are accepted (prefix semantics);
!  5. a calibrated species under another name with the same composition is accepted by a hand-written routine;
!  6. without composition data the molecular weight of phase.txt stands in for the composition;
!  7. a name that is not hooked is not checked (general fallback);
!  8-9. the mechanism name is the whole first line of chemistry-info.txt;
!  11. a generated routine (Gerlinger9) also needs the slot NAME: H2OX in the H2O slot is refused;
!  10. in child processes (this program with the argument child-strict / child-refusal): the strict
!      fallback policy (FLINT_STRICT_MECHANISM=1 or Yes stops an unhooked name, =0 falls back to general,
!      an unrecognised value is reported on the error unit and ignored) and
!      the contract refusal stop the process, and their messages are on the error unit;
!  12. in child processes (child-early, child-early-refusal, child-early-unloaded): Assign_Mechanism
!      called before the tables are loaded prints a WARNING on both units and selects the routine; the
!      first chemistry call checks the contract (same omegadot as with the usual order), stops with
!      both lists on a mismatch, and stops with an [ERROR] line when the tables are still not loaded;
!  13. the committed src/lib/Lib_Chemistry_contract.f90 is the output of its generator
!      (python3 ../utils/mechanism_contract.py --check; skipped when python3 is not on the PATH).
! Needs no Cantera. Exit code 1 on failure.
program test
  use FLINT_Lib_Thermodynamic
  use FLINT_Load_ThermoTransport
  use FLINT_Load_chemistry
  use FLINT_Lib_Chemistry_data
  use FLINT_Lib_Chemistry_wdot
  use FLINT_Lib_Chemistry_contract
  implicit none
  character(32) :: mech_name
  character(len=s_str_len), allocatable :: names0(:), names1(:)
  real(8), allocatable :: comp0(:,:), comp1(:,:)
  integer :: err, ns0, nfail
  logical :: ok
  character(len=512) :: self, arg
  integer :: err11, i11
  real(8), allocatable :: r12(:), w12a(:), w12b(:)

  nfail = 0
  call get_command_argument(0, self)
  call get_command_argument(1, arg)
  if (arg == 'child-strict') then
    ! an unhooked name ('Nassini Original', written by case 9) through Assign_Mechanism
    err = read_idealgas_thermo('../database/WD/')
    err = read_chemistry(folder='tables/name-nassini', mech_name=mech_name)
    if (err /= 0) then; write(*,'(A,I0)') '[FAIL] child-strict: read_chemistry ios=', err; stop 2; endif
    call Assign_Mechanism(mech_name)
    write(*,'(A)') 'child-strict: fallback to the general procedure'
    stop
  else if (arg == 'child-refusal') then
    ! the old alphabetical order (slots 2 and 5 swapped) through Assign_Mechanism
    err = read_idealgas_thermo('../database/WD/')
    err = read_chemistry(folder='../database/WD/', mech_name=mech_name)
    if (err /= 0) then; write(*,'(A,I0)') '[FAIL] child-refusal: read_chemistry ios=', err; stop 2; endif
    names0 = species_names; comp0 = species_composition
    species_names(2) = names0(5); species_names(5) = names0(2)
    species_composition(:,2) = comp0(:,5); species_composition(:,5) = comp0(:,2)
    call Assign_Mechanism(mech_name)
    write(*,'(A)') 'child-refusal: the contract did not stop the process'
    stop
  else if (arg == 'child-early') then
    ! the mechanism is selected before the tables are loaded: WARNING, check at the first call
    call Assign_Mechanism('WD')
    err = read_idealgas_thermo('../database/WD/')
    err = read_chemistry(folder='../database/WD/', mech_name=mech_name)
    if (err /= 0) then; write(*,'(A,I0)') '[FAIL] child-early: read_chemistry ios=', err; stop 2; endif
    allocate(r12(ns), w12a(ns), w12b(ns))
    r12 = 0.05d0; w12a = 7.0d0
    call chemistry_source(r12, 1500.5d0, w12a)      ! first call: the contract is checked here
    call Assign_Mechanism('WD')                       ! the usual order: tables, then the mechanism
    r12 = 0.05d0; w12b = 7.0d0
    call chemistry_source(r12, 1500.5d0, w12b)
    if (any(w12a /= w12b) .or. all(w12a == 0d0)) then
      write(*,'(A)') 'child-early: the first call differs from the usual order'; stop 3
    endif
    write(*,'(A)') 'child-early: the first call checked the contract and gave the omegadot of the usual order'
    stop
  else if (arg == 'child-early-refusal') then
    call Assign_Mechanism('WD')
    err = read_idealgas_thermo('../database/WD/')
    err = read_chemistry(folder='../database/WD/', mech_name=mech_name)
    if (err /= 0) then; write(*,'(A,I0)') '[FAIL] child-early-refusal: read_chemistry ios=', err; stop 2; endif
    nrc_arrh = nrc_arrh + 1                           ! a reaction count the WD routine does not have
    allocate(r12(ns), w12a(ns)); r12 = 0.05d0
    call chemistry_source(r12, 1500.5d0, w12a)
    write(*,'(A)') 'child-early-refusal: the first call did not stop the process'
    stop
  else if (arg == 'child-early-unloaded') then
    call Assign_Mechanism('WD')
    allocate(r12(8), w12a(8)); r12 = 0.05d0
    call chemistry_source(r12, 1500.5d0, w12a)      ! still no tables at the first call
    write(*,'(A)') 'child-early-unloaded: the first call did not stop the process'
    stop
  endif
  err = read_idealgas_thermo('../database/WD/')
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_idealgas_thermo ../database/WD: ios=', err; stop 1; endif
  err = read_chemistry(folder='../database/WD/', mech_name=mech_name)
  if (err /= 0) then; write(*,'(A,I0)') '[FAIL] read_chemistry ../database/WD: ios=', err; stop 1; endif
  ns0 = ns; names0 = species_names; comp0 = species_composition

  ! 1. routine order: pass, and Assign_Mechanism hooks the routine
  call check_mechanism_contract(mech_name, ok)
  call verdict('1 database/WD in the WD.f90 order accepted', ok)
  call Assign_Mechanism(mech_name)
  call verdict('1b Assign_Mechanism(WD) hooked a routine', associated(chemistry_source))

  ! 2. old alphabetical order [CH4, CO, CO2, H2O, O2]: slot 2 is CO instead of O2 -> refused
  species_names(2) = names0(5); species_names(5) = names0(2)
  species_composition(:,2) = comp0(:,5); species_composition(:,5) = comp0(:,2)
  call check_mechanism_contract(mech_name, ok)
  call verdict('2 swapped slots 2/5 (old alphabetical fixture) refused', .not. ok)
  species_names = names0; species_composition = comp0

  ! 3. reaction count
  nrc_arrh = nrc_arrh + 1
  call check_mechanism_contract(mech_name, ok)
  call verdict('3 nrc_arrh + 1 refused', .not. ok)
  nrc_arrh = nrc_arrh - 1

  ! 4. appended species (prefix semantics): N2 after the 5 routine slots
  allocate(names1(ns0+1), comp1(size(comp0,1), ns0+1))
  names1(1:ns0) = names0; names1(ns0+1) = 'N2'
  comp1(:,1:ns0) = comp0; comp1(:,ns0+1) = 0.d0
  call move_alloc(names1, species_names); call move_alloc(comp1, species_composition); ns = ns0 + 1
  call check_mechanism_contract(mech_name, ok)
  call verdict('4 one inert species appended after the routine slots accepted', ok)
  species_names = names0; species_composition = comp0; ns = ns0

  ! 5. calibrated species: same composition, other name (hand-written routine)
  species_names(4) = 'H2OCalib'
  call check_mechanism_contract(mech_name, ok)
  call verdict('5 H2O renamed H2OCalib with the same composition accepted', ok)
  species_names = names0

  ! 5b. element symbols in another case ('CL', 'c') describe the same composition
  elements_names(1) = 'c'
  call check_mechanism_contract(mech_name, ok)
  call verdict('5b element symbol in lower case accepted', ok)
  elements_names(1) = 'C'

  ! 6. no composition data (INPUT folder written before the table writer wrote composition.txt): molecular weight
  deallocate(species_composition)
  species_names(4) = 'H2OCalib'
  call check_mechanism_contract(mech_name, ok)
  call verdict('6a without composition data, H2OCalib with the weight of H2O accepted', ok)
  species_names = names0
  wm_tab(4) = 34.014d0   ! H2O2 in the H2O slot
  call check_mechanism_contract(mech_name, ok)
  call verdict('6b without composition data, a species with another weight refused', .not. ok)
  wm_tab(4) = 18.015d0
  call check_mechanism_contract(mech_name, ok)
  call verdict('6c without composition data, exact names and weights accepted', ok)
  species_names(5) = 'N2'; wm_tab(5) = 28.014d0   ! N2 in the CO slot: weight within 0.05 kg/kmol
  call check_mechanism_contract(mech_name, ok)
  call verdict('6d without composition data, N2 in the CO slot (weight 28.014 vs 28.010) refused', .not. ok)
  species_names = names0; wm_tab(5) = 28.010d0
  species_names(2) = 'o2'
  call check_mechanism_contract(mech_name, ok)
  call verdict('6e without composition data, slot name in another case accepted', ok)
  species_names = names0
  species_composition = comp0

  ! 7. a name that is not hooked is not checked
  call check_mechanism_contract('nemo', ok)
  call verdict('7 unhooked name: nothing to check', ok)
  call free_chemistry_data()

  ! 8. the mechanism name is the whole first line of chemistry-info.txt ('WD 2026' is not 'WD');
  !    the rate table is the database/WD one (on the grid of the loaded thermo tables)
  call execute_command_line('mkdir -p tables/name-blank && cp ../database/WD/chemistry-Arrhenius.dat tables/name-blank/')
  open(newunit=err, file='tables/name-blank/chemistry-info.txt', status='replace', action='write')
  write(err,'(A)') '  WD 2026  '
  write(err,'(A)') 'N.ro species = 5'
  write(err,'(A)') 'N.ro reactions = 3'
  write(err,'(A)') ''
  write(err,'(A)') 'General loop info. They are not used if the mechanism is exlplicitly defined'
  write(err,'(A)') ''
  write(err,'(A)') 'Reaction type'
  write(err,'(A)') '1 Arrhenius'
  write(err,'(A)') '2 Arrhenius'
  write(err,'(A)') '3 Arrhenius'
  write(err,'(A)') ''
  write(err,'(A)') 'Reaction definition'
  close(err)
  err = read_chemistry(folder='tables/name-blank', mech_name=mech_name)
  call verdict('8a name with a blank: read_chemistry returns 0', err == 0)
  call verdict('8b name with a blank read whole and trimmed: '//trim(mech_name), mech_name == 'WD 2026')
  call check_mechanism_contract(mech_name, ok)
  call verdict('8c the name with a blank is not a hooked name (general fallback)', ok)
  call free_chemistry_data()
  ! 8d. TAB and CR are blanks, as they were for the list-directed read ('WD<TAB><CR>' is 'WD')
  open(newunit=err, file='tables/name-blank/chemistry-info.txt', status='replace', action='write')
  write(err,'(A)') achar(9)//'WD'//achar(9)//achar(13)
  write(err,'(A)') 'N.ro species = 5'
  write(err,'(A)') 'N.ro reactions = 3'
  write(err,'(A)') ''
  write(err,'(A)') 'General loop info. They are not used if the mechanism is exlplicitly defined'
  write(err,'(A)') ''
  write(err,'(A)') 'Reaction type'
  write(err,'(A)') '1 Arrhenius'
  write(err,'(A)') '2 Arrhenius'
  write(err,'(A)') '3 Arrhenius'
  write(err,'(A)') ''
  write(err,'(A)') 'Reaction definition'
  close(err)
  err = read_chemistry(folder='tables/name-blank', mech_name=mech_name)
  call verdict('8d TAB/CR around the name are blanks: ['//trim(mech_name)//']', err == 0 .and. mech_name == 'WD')
  call free_chemistry_data()

  ! 9. 'Nassini Original': a list-directed read gave 'Nassini' (hooked: Nassini_4, 3 slots O2/H2O/H2),
  !    the whole-line read gives an unhooked name; the WD data would not pass the Nassini_4 contract
  call execute_command_line('mkdir -p tables/name-nassini && cp ../database/WD/chemistry-Arrhenius.dat tables/name-nassini/')
  call execute_command_line("sed '1s/.*/Nassini Original/' tables/name-blank/chemistry-info.txt > tables/name-nassini/chemistry-info.txt")
  err = read_chemistry(folder='tables/name-nassini', mech_name=mech_name)
  call verdict('9a Nassini Original: read whole: '//trim(mech_name), err == 0 .and. mech_name == 'Nassini Original')
  call check_mechanism_contract(mech_name, ok)
  call verdict('9b Nassini Original is not a hooked name: nothing to check', ok)
  call check_mechanism_contract('Nassini', ok)
  call verdict('9c the truncated name Nassini is hooked and the WD data are refused by its contract', .not. ok)
  call free_chemistry_data()

  ! 11. a GENERATED routine needs the slot name too (database/Gerlinger: no composition.txt, so the slot test
  !     is weight + name prefix): 'H2OX' passes that test for the H2O slot but not the exact-name test
  call free_chemistry_data()
  if (allocated(species_names)) deallocate(species_names)
  if (allocated(wm_tab)) deallocate(wm_tab)
  if (allocated(Ri_tab)) deallocate(Ri_tab)
  if (allocated(h_tab)) deallocate(h_tab)
  if (allocated(s_tab)) deallocate(s_tab)
  if (allocated(cp_tab)) deallocate(cp_tab)
  if (allocated(dcpi_tab)) deallocate(dcpi_tab)
  err = read_idealgas_thermo('../database/Gerlinger/')
  err11 = read_chemistry(folder='../database/Gerlinger/', mech_name=mech_name)
  call check_mechanism_contract('Gerlinger-9', ok)
  ! (read_idealgas_thermo returns 5 when composition.txt is absent, as for this legacy folder)
  write(*,'(A,I0,A,I0)') '    database/Gerlinger: read_idealgas_thermo ios = ', err, ', read_chemistry ios = ', err11
  call verdict('11a database/Gerlinger accepted by the generated routine Gerlinger9', (err == 0 .or. err == 5) .and. err11 == 0 .and. ok)
  do i11 = 1, ns
    if (trim(species_names(i11)) == 'H2O') exit
  enddo
  if (i11 <= ns) species_names(i11) = 'H2OX'
  call check_mechanism_contract('Gerlinger-9', ok)
  call verdict('11b generated routine: a slot species under another name (H2OX for H2O) is refused', i11 <= ns .and. .not. ok)
  if (i11 <= ns) species_names(i11) = 'H2O'

  ! 10. child processes: strict fallback policy and the channels of the refusals
  call execute_command_line('FLINT_STRICT_MECHANISM=1 '//trim(self)// &
    ' child-strict > tables/child-strict.out 2> tables/child-strict.err', exitstat=err)
  call verdict('10a strict mode: an unhooked name stops the child process (exit code /= 0)', err /= 0)
  call execute_command_line("command grep -q 'strict mode is on' tables/child-strict.err", exitstat=err)
  call verdict('10b strict mode: the [ERROR] line is on the error unit', err == 0)
  call execute_command_line('FLINT_STRICT_MECHANISM=0 '//trim(self)// &
    ' child-strict > tables/child-nostrict.out 2> tables/child-nostrict.err', exitstat=err)
  call verdict('10c strict mode off: the child falls back to the general procedure (exit code 0)', err == 0)
  call execute_command_line("command grep -q 'defaulting to the general procedure' tables/child-nostrict.err", exitstat=err)
  call verdict('10d strict mode off: the fallback WARNING is on the error unit', err == 0)
  call execute_command_line('FLINT_STRICT_MECHANISM=Yes '//trim(self)// &
    ' child-strict > tables/child-strict-yes.out 2> tables/child-strict-yes.err', exitstat=err)
  call verdict('10h strict mode: the value is case-insensitive (Yes stops the child, exit code /= 0)', err /= 0)
  call execute_command_line('FLINT_STRICT_MECHANISM=maybe '//trim(self)// &
    ' child-strict > tables/child-strict-maybe.out 2> tables/child-strict-maybe.err', exitstat=err)
  call verdict('10i strict mode: an unrecognised value leaves the fallback (exit code 0)', err == 0)
  call execute_command_line("command grep -q 'FLINT_STRICT_MECHANISM=.maybe. is not one of' tables/child-strict-maybe.err", &
    exitstat=err)
  call verdict('10j strict mode: an unrecognised value is reported on the error unit', err == 0)
  call execute_command_line(trim(self)//' child-refusal > tables/child-refusal.out 2> tables/child-refusal.err', exitstat=err)
  call verdict('10e contract refusal: Assign_Mechanism stops the child process (exit code /= 0)', err /= 0)
  call execute_command_line("command grep -q 'expected: 5 species' tables/child-refusal.err", exitstat=err)
  call verdict('10f contract refusal: the two lists are on the error unit too', err == 0)
  call execute_command_line("command grep -q 'expected: 5 species' tables/child-refusal.out", exitstat=err)
  call verdict('10g contract refusal: the two lists are on standard output', err == 0)

  ! 12. Assign_Mechanism before the tables are loaded (child processes)
  call execute_command_line(trim(self)//' child-early > tables/child-early.out 2> tables/child-early.err', exitstat=err)
  call verdict('12a hook before the tables: the first call checks the contract and gives the omegadot of the usual order', &
    err == 0)
  call execute_command_line("command grep -q 'its contract is checked at the first chemistry call' tables/child-early.err", &
    exitstat=err)
  call verdict('12b hook before the tables: the WARNING is on the error unit', err == 0)
  call execute_command_line("command grep -q 'its contract is checked at the first chemistry call' tables/child-early.out", &
    exitstat=err)
  call verdict('12c hook before the tables: the WARNING is on standard output', err == 0)
  call execute_command_line(trim(self)//' child-early-refusal > tables/child-early-refusal.out 2> tables/child-early-refusal.err', &
    exitstat=err)
  call verdict('12d hook before the tables, data that do not match: the first call stops the process', err /= 0)
  call execute_command_line("command grep -q 'expected: 5 species' tables/child-early-refusal.err", exitstat=err)
  call verdict('12e hook before the tables, data that do not match: the two lists are on the error unit', err == 0)
  call execute_command_line(trim(self)//' child-early-unloaded > tables/child-early-unloaded.out 2> tables/child-early-unloaded.err', &
    exitstat=err)
  call verdict('12f hook, no tables at the first call: the process stops', err /= 0)
  call execute_command_line("command grep -q 'species and rate tables are not loaded' tables/child-early-unloaded.err", exitstat=err)
  call verdict('12g hook, no tables at the first call: the [ERROR] line is on the error unit', err == 0)

  ! 13. the committed contract module is the output of utils/mechanism_contract.py
  call execute_command_line('command -v python3 > /dev/null 2>&1', exitstat=err)
  if (err == 0) then
    call execute_command_line('python3 ../utils/mechanism_contract.py --check > tables/contract-generator.out 2>&1', &
      exitstat=err)
    call verdict('13 Lib_Chemistry_contract.f90 is the output of utils/mechanism_contract.py '// &
      '(diff in tables/contract-generator.out)', err == 0)
  else
    write(*,'(A)') ' [skip] 13 python3 not found: generator equality not checked'
  endif
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
end program test
