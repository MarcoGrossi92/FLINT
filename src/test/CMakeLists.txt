# FLINT test suite: the drivers of src/test/<area>/ (built in bin/test/) and their CTest tests.
# A test <area>/<name> runs in test/<area>/, where its fixtures and outputs are, with two labels:
#   area  thermo, chemistry, equilibrium, batch               (ctest -L chemistry)
#   tier  quick: the subset run by the GitHub workflow       (ctest -L quick)
#         long:  the rest                                    (ctest runs both)
set(TEST_AREAS thermo chemistry equilibrium batch)
set(FLINT_TEST_DIR "${CMAKE_SOURCE_DIR}/test")
# another build (e.g. with Cantera) keeps its drivers apart with -DFLINT_TEST_BINDIR=<path>
set(FLINT_TEST_BINDIR "${CMAKE_SOURCE_DIR}/bin/test" CACHE PATH "Folder of the test drivers")
set(EXECUTABLE_OUTPUT_PATH "${FLINT_TEST_BINDIR}")

# -------- Drivers --------
# Every driver links the library FLINT, compiled once (a driver compiling the sources itself would write
# the same module files as the others: a parallel build collides on them)
set(DRIVER_LIBS FLINT)
if (USE_TECIO)
  list(APPEND DRIVER_LIBS ${LINKLIBS})
endif()

foreach(AREA ${TEST_AREAS})
  # the working directory of the tests of the area (git keeps it only while it holds fixtures)
  file(MAKE_DIRECTORY "${FLINT_TEST_DIR}/${AREA}")
  file(GLOB DRIVERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${AREA}/*.f90" "${CMAKE_CURRENT_SOURCE_DIR}/${AREA}/*.F90")
  foreach(DRIVER ${DRIVERS})
    get_filename_component(EXE_NAME ${DRIVER} NAME_WE)
    add_executable(${EXE_NAME} ${DRIVER})
    # the drivers are Fortran programs: link them with the Fortran driver even when the TecIO
    # C++ library is in LINKLIBS (with ifx the CXX link dropped the Fortran main: undefined main)
    set_target_properties(${EXE_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
    # FLINT and its dependencies (see cmake/LinkDependencies.cmake): OSLO, ORION, TecIO and Cantera (optional)
    target_link_libraries(${EXE_NAME} PRIVATE ${DRIVER_LIBS})
  endforeach()
endforeach()

# Cantera C++ batch reactor: the references of the batch tests
if (USE_CANTERA)
  add_executable(test-batchCXX batch/test-batchCXX.cpp)
  set_target_properties(test-batchCXX PROPERTIES CXX_STANDARD 17)
  # one SUNDIALS installation
  target_include_directories(test-batchCXX PRIVATE "${SUNDIALS_INSTALL_DIR}/include/"
                                                   "${canteraD}/include/")
  link_FLINT_dependencies(test-batchCXX)
  target_link_libraries(test-batchCXX PRIVATE ORION)

  # test/batch/reference/<case>.dat:  cmake --build . --target batch-reference
  # after a change to test/batch/cases.txt, to a mechanism of database/ or to Cantera; commit the files it writes.
  add_custom_target(batch-reference
    COMMAND test-batchCXX --reference
    WORKING_DIRECTORY "${FLINT_TEST_DIR}/batch"
    COMMENT "Writing the Cantera references of the batch tests in test/batch/reference")

  # Cantera C++ equilibrium: the references of the sweeps of test-CEA, test/equilibrium/reference/<mechanism>.dat
  #   cmake --build . --target equilibrium-reference
  add_executable(test-equilCXX equilibrium/test-equilCXX.cpp)
  set_target_properties(test-equilCXX PROPERTIES CXX_STANDARD 17)
  target_include_directories(test-equilCXX PRIVATE "${SUNDIALS_INSTALL_DIR}/include/"
                                                   "${canteraD}/include/")
  link_FLINT_dependencies(test-equilCXX)
  target_link_libraries(test-equilCXX PRIVATE ORION)
  add_custom_target(equilibrium-reference
    COMMAND test-equilCXX
    WORKING_DIRECTORY "${FLINT_TEST_DIR}/equilibrium"
    COMMENT "Writing the Cantera references of the equilibrium sweeps in test/equilibrium/reference")
endif()

# -------- CTest --------
# flint_add_test(<area> <name> <tier> <driver> [arguments...]): the test <area>/<name>
function(flint_add_test AREA NAME TIER)
  add_test(NAME ${AREA}/${NAME} COMMAND ${ARGN} WORKING_DIRECTORY "${FLINT_TEST_DIR}/${AREA}")
  set_tests_properties(${AREA}/${NAME} PROPERTIES LABELS "${AREA};${TIER}" TIMEOUT 600)
endfunction()

# Thermodynamics: the real-fluid (p, h) tables and the inverse (p, T) -> h
# (test-thermo-ideal: a Cantera benchmark, not a test)
flint_add_test(thermo      thermo-real  quick test-thermo-real)

# Chemistry: tables, mechanism contract, rate laws against Cantera references
# (test-wdot: Cantera at run time, not a test)
flint_add_test(chemistry   tables       quick test-tables)
flint_add_test(chemistry   ranges       quick test-ranges)
flint_add_test(chemistry   rhs-range    quick test-rhs-range)
flint_add_test(chemistry   contract     quick test-contract)
flint_add_test(chemistry   inert        quick test-inert)
flint_add_test(chemistry   falloff      quick test-falloff)
flint_add_test(chemistry   orders       quick test-orders)
flint_add_test(chemistry   orders-warn  quick test-orders-warn)
flint_add_test(chemistry   stoich       quick test-stoich)
flint_add_test(chemistry   andersen     quick test-andersen)
# both make orders/noblock
set_tests_properties(chemistry/orders chemistry/orders-warn PROPERTIES RESOURCE_LOCK orders-noblock)

# Equilibrium: CEA against Cantera values
flint_add_test(equilibrium CEA          quick test-CEA)

# Batch reactor: one test per case of test/batch/cases.txt, test-batchF check <case> against
# test/batch/reference/<case>.dat. Quick cases: a global mechanism with orders (WD), tables in
# TecIO's .szplt (Troyes, TSR-GP-24), Troe (TSR-GP-24), Troe and Lindemann (Pelucchi), H2/air (Gerlinger).
set(BATCH_QUICK WD Troyes TSR-GP-24 Pelucchi Gerlinger)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${FLINT_TEST_DIR}/batch/cases.txt")
file(STRINGS "${FLINT_TEST_DIR}/batch/cases.txt" BATCH_LINES REGEX "^[^# \t]")
foreach(LINE ${BATCH_LINES})
  string(REGEX MATCH "^[^ \t]+" CASE "${LINE}")
  if (CASE IN_LIST BATCH_QUICK)
    flint_add_test(batch ${CASE} quick test-batchF check ${CASE})
  else()
    flint_add_test(batch ${CASE} long  test-batchF check ${CASE})
  endif()
endforeach()

# Tables in TecIO's .szplt (Troyes, TSR-GP-24, TSR-Rich-31) are read with TecIO only
if (NOT USE_TECIO)
  file(GLOB SZPLT_FILES "${CMAKE_SOURCE_DIR}/database/*/*.szplt")
  foreach(FILE ${SZPLT_FILES})
    get_filename_component(DIR ${FILE} DIRECTORY)
    get_filename_component(CASE ${DIR} NAME)
    if (TEST batch/${CASE})
      set_tests_properties(batch/${CASE} PROPERTIES DISABLED TRUE)
    endif()
  endforeach()
  set_tests_properties(equilibrium/CEA PROPERTIES DISABLED TRUE)   # TSR-GP-24
endif()
