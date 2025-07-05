
find_library(ACCELERATE_FRAMEWORK Accelerate)

file(GLOB_RECURSE TESTS *.f90 *.F90)
set(TESTS ${TESTS} PARENT_SCOPE)
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin/test")

foreach(TEST ${TESTS})
  get_filename_component(EXE_NAME ${TEST} NAME_WE)
  add_executable(${EXE_NAME} ${TEST})
  target_link_libraries(${EXE_NAME} cantera cantera_fortran ${ACCELERATE_FRAMEWORK} c++)
  set_target_properties(${EXE_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin/)
endforeach()
