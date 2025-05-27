set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin/test")
add_executable(demo demo.f90)
target_link_libraries(demo cantera_fortran cantera Threads::Threads stdc++)