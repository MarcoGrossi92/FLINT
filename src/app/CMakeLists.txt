find_library(ACCELERATE_FRAMEWORK Accelerate)

add_executable(batch "batch.f90")
target_link_libraries(batch PRIVATE ${lib})

target_include_directories(${lib}
  INTERFACE
      $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

IF(USE_OPENMP)
    SET_TARGET_PROPERTIES(${lib} PROPERTIES
                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
ENDIF(USE_OPENMP)