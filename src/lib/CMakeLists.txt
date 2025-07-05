find_library(ACCELERATE_FRAMEWORK Accelerate)

file(GLOB LIBSRC *.f90)
add_library(${lib} ${LIBSRC})
target_link_libraries(${lib} ORION OSlo)
if (USE_TECIO)
  target_link_libraries(${lib} ${LINKLIBS})
endif()
target_link_libraries(${lib} cantera cantera_fortran ${ACCELERATE_FRAMEWORK} c++)

target_include_directories(${lib}
  INTERFACE
      $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

IF(USE_OPENMP)
    SET_TARGET_PROPERTIES(${lib} PROPERTIES
                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
ENDIF(USE_OPENMP)