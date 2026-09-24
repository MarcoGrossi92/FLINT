! test-contract: the mechanism contract check of Assign_Mechanism.
! Loads database/WD (species in the routine order of WD.f90) and checks that
!  1. the loaded data pass the check and Assign_Mechanism hooks WD;
!  2. the alphabetical order of the old fixture (slots 2 and 5 swapped: CO where the routine reads O2) is refused;
!  3. a different number of Arrhenius reactions is refused;
!  4. species appended after the routine slots are accepted (prefix semantics);
!  5. a calibrated species under another name with the same composition is accepted by a hand-written routine;
!  6. without composition data the molecular weight of phase.txt stands in for the composition;
!  7. a name that is not hooked is not checked (general fallback).
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

  nfail = 0
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
