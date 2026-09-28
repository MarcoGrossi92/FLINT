file(GLOB_RECURSE TESTS *.f90 *.F90)
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin/test")

foreach(TEST ${TESTS})
  get_filename_component(EXE_NAME ${TEST} NAME_WE)
  add_executable(${EXE_NAME} ${TEST} ${FLINT_SRC})
  # the drivers are Fortran programs: link them with the Fortran driver even when the TecIO
  # C++ library is in LINKLIBS (with ifx the CXX link dropped the Fortran main: undefined main)
  set_target_properties(${EXE_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
  target_include_directories(${EXE_NAME} INTERFACE $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

  # FLINT dependencies (see cmake/LinkDependencies.cmake): OSLO and Cantera (optional)
  link_FLINT_dependencies(${EXE_NAME})

  # ORION
  target_link_libraries(${EXE_NAME} PRIVATE ORION)
  if (USE_TECIO)
    target_link_libraries(${EXE_NAME} PRIVATE ${LINKLIBS})
  endif()

endforeach()
