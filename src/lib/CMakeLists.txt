find_library(ACCELERATE_FRAMEWORK Accelerate)

file(GLOB LIBSRC *.f90 Lib_ChemMech/*.f90)
add_library(${lib} ${LIBSRC})
target_link_libraries(${lib} PRIVATE ORION OSlo)
if (USE_TECIO)
  target_link_libraries(${lib} PRIVATE ${LINKLIBS})
endif()
if (USE_CANTERA)
  target_link_libraries(${lib} PRIVATE cantera cantera_fortran ${ACCELERATE_FRAMEWORK} c++)
endif()

target_include_directories(${lib}
  INTERFACE
      $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

IF(USE_OPENMP)
    SET_TARGET_PROPERTIES(${lib} PROPERTIES
                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
ENDIF(USE_OPENMP)