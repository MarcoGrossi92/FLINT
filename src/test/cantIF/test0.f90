program testN0
  use U_Lib_Thermodynamic
  use U_IO_Table
  use cantera
  use U_cantera_load
  implicit none
  type(phase_t) :: phase

  call Read_IdealGas_Properties()
  write(*,*) 'Native species number :',nsc

  call load_phase(phase, 'WD.yaml')
  write(*,*) 'Cantera species number:',nSpecies(phase)

end program testN0