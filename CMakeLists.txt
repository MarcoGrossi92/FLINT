cmake_minimum_required(VERSION 3.23)

# -------- Say my name! --------
set(lib "FLINT")

project (FLINT LANGUAGES C CXX Fortran)

OPTION(USE_CANTERA "Use Cantera" OFF)
# Before the compile flags and the dependencies, which read them: OSLO's own options would turn OpenMP and MPI on in
# the whole build tree. A project that adds FLINT sets them before, and these do nothing
OPTION(USE_OPENMP "Use OpenMP for parallelization" OFF)
OPTION(USE_MPI "Use the MPI library for parallelization" OFF)

if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
  set(CMAKE_VERBOSE_MAKEFILE FALSE)

  # Search path for cmake modules
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

  # This INCLUDE statement executes code that sets the compile flags for DEBUG,
  # RELEASE, and TESTING.
  INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 

endif()

# -------- Dependencies linking --------
# By path: a module path of the project that adds FLINT comes first, and may have a file of the same name
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/LinkDependencies.cmake)

# -------- ORION and OSLO --------
# The newest releases that satisfy the requirements below. When FLINT is part of a larger project, the releases of its
# build tree, which must satisfy them too. See cmake/SemVerDeps.cmake: the releases stay those of the lock file
# deps.lock, -DSEMVER_UPDATE=ALL selects the newest ones again, -DFETCHCONTENT_SOURCE_DIR_<NAME>=<path> builds a local
# checkout instead. The sources are checked out in lib/OSLO and lib/ORION, out of the build folder, which keeps them
# from a build to the next
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SemVerDeps.cmake)

# Cantera is built on the SUNDIALS of OSLO
if (USE_CANTERA)
  set ( USE_SUNDIALS ON )
endif()

# Without the submodule lib/sundials, which OSLO clones only if it has to build SUNDIALS
semver_dependency(OSLO GIT_REPOSITORY https://github.com/MarcoGrossi92/OSLO.git VERSION 1.1 GIT_SUBMODULES ""
                  SOURCE_DIR lib/OSLO)
# 2.1: ORION_TECIO_DIR, the shared TecIO of install.sh and of the GitHub workflow
semver_dependency(ORION GIT_REPOSITORY https://github.com/MarcoGrossi92/ORION.git VERSION 2.1 SOURCE_DIR lib/ORION)
FetchContent_MakeAvailable(OSLO ORION)

# Cantera
if (USE_CANTERA)
  # The SUNDIALS installation of OSLO
  get_target_property(SUNDIALS_INSTALL_DIR OSLO OSLO_SUNDIALS_DIR)
  if (NOT SUNDIALS_INSTALL_DIR)
    message(FATAL_ERROR "Cantera needs SUNDIALS, but OSLO was added to the build tree with USE_SUNDIALS=OFF")
  endif()
  find_package(SUNDIALS REQUIRED PATHS ${SUNDIALS_INSTALL_DIR} NO_DEFAULT_PATH)

  set(canteraD "${CMAKE_CURRENT_SOURCE_DIR}/lib/cantera")
  # Check if cloned
  if(NOT EXISTS "${canteraD}/README.rst")
    message(STATUS "Cantera not found. Attempting to initialize...")
    execute_process(COMMAND git submodule update --init lib/cantera
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL 0)
      message(FATAL_ERROR "Failed to initialize git submodules.")
    endif()
  endif()
  # Check if already built
  if(NOT EXISTS "${canteraD}/build/lib/libcantera.a" OR NOT EXISTS "${canteraD}/build/lib/libcantera_fortran.a"
     OR NOT EXISTS "${canteraD}/build/src/fortran/cantera.mod")
    message(STATUS "Cantera not installed, building and installing it...")
    # Cantera needs the headers of Boost 1.70 or newer to build. If the system has no Boost, download the headers of 1.86.0, 
    # which is the latest release that Cantera 2.6.0 supports
    find_package(Boost 1.70)
    if (Boost_FOUND)
      set(boost_inc_dir ${Boost_INCLUDE_DIR})
    else()
      message(STATUS "Boost not found, downloading the headers of Boost 1.86.0...")
      # The release has no CMakeLists.txt: FetchContent_MakeAvailable only extracts it
      cmake_policy(PUSH)
      if (POLICY CMP0135)
        cmake_policy(SET CMP0135 NEW)
      endif()
      FetchContent_Declare(boost_headers
        URL https://github.com/boostorg/boost/releases/download/boost-1.86.0/boost-1.86.0-b2-nodocs.tar.xz
        URL_HASH SHA256=a4d99d032ab74c9c5e76eddcecc4489134282245fffa7e079c5804b92b45f51d)
      FetchContent_MakeAvailable(boost_headers)
      cmake_policy(POP)
      set(boost_inc_dir ${boost_headers_SOURCE_DIR})
    endif()
    # The library folder of the SUNDIALS installation: lib or lib64
    find_library(sundials_core_lib sundials_core PATHS ${SUNDIALS_INSTALL_DIR} PATH_SUFFIXES lib lib64
                 NO_DEFAULT_PATH NO_CACHE REQUIRED)
    get_filename_component(sundials_libdir ${sundials_core_lib} DIRECTORY)
    # Without CONDA_PREFIX: in an active conda environment Cantera would add its include and library paths, and run
    # its programs with the libstdc++ of the environment, which may be older than that of the compiler
    set(cantera_scons ${CMAKE_COMMAND} -E env --unset=CONDA_PREFIX scons -C ${canteraD})
    include(ExternalProject)
    ExternalProject_Add(
        cantera_world
        PREFIX cantera_build
        SOURCE_DIR ${canteraD}
        BINARY_DIR ${CMAKE_BINARY_DIR}/cantera_build
        CONFIGURE_COMMAND ""  # Cantera does not use configure step
        BUILD_COMMAND ${cantera_scons} build
            prefix=${CMAKE_BINARY_DIR}/cantera_install
            f90_interface=y
            FORTRAN=gfortran
            clib_legacy=y
            system_sundials=y
            sundials_include=${SUNDIALS_INSTALL_DIR}/include
            sundials_libdir=${sundials_libdir}
            boost_inc_dir=${boost_inc_dir}
            system_eigen=n
            system_fmt=n
            system_yamlcpp=n
        INSTALL_COMMAND ${cantera_scons} install
    )
  else()
    message(STATUS "Cantera already installed.")
  endif()
  link_directories("${canteraD}/build/lib")
  include_directories("${canteraD}/build/src/fortran/")
  add_compile_definitions(CANTERA)
  add_library(FLINT_Cantera INTERFACE)
  target_link_libraries(FLINT_Cantera INTERFACE cantera_fortran cantera
                        SUNDIALS::cvodes SUNDIALS::idas SUNDIALS::nvecserial SUNDIALS::core)
endif(USE_CANTERA)

# -------- Fortran directories --------
include(GNUInstallDirs)
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# -------- FLINT sources and library --------
# A consumer uses FLINT in one of two ways, the same one in the whole build tree (both ways in one tree compile the
# sources twice, into the same module files):
#   - the sources, FLINT_SRC: compiled in a target of the consumer, with link_FLINT_dependencies(<target>). FLINT_SRC
#     is set in the directory that adds FLINT only
#   - the static library FLINT (alias FLINT::FLINT): linked, from any directory of the tree. It is compiled once per
#     tree, however many targets link it, with the flags of the top-level project, which must set them before FLINT
#     is added: CMAKE_Fortran_FLAGS_<CONFIG>, CMAKE_EXE_LINKER_FLAGS_<CONFIG> and, for Intel IPO, the archiver
#     CMAKE_AR/CMAKE_RANLIB (llvm-ar/llvm-ranlib). Its objects are then those a consumer would compile from FLINT_SRC,
#     and with -ipo the link of each executable optimises them together with the consumer's own code
file(GLOB LIBSRC CONFIGURE_DEPENDS src/lib/*.f90 src/lib/Lib_ChemMech/*.f90)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
  set(FLINT_SRC ${LIBSRC})
else()
  set(FLINT_SRC ${LIBSRC} PARENT_SCOPE)
endif()

# Built only for the targets that link it (EXCLUDE_FROM_ALL): a consumer of FLINT_SRC does not build it as well
add_library(FLINT STATIC EXCLUDE_FROM_ALL ${LIBSRC})
add_library(FLINT::FLINT ALIAS FLINT)
target_include_directories(FLINT PUBLIC $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)
# PUBLIC: FLINT's modules use ORION's and OSLO's, and ifx reads those module files when it compiles a user of FLINT's
target_link_libraries(FLINT PUBLIC ORION OSLO)
link_FLINT_dependencies(FLINT)
# The OpenMP flag on FLINT's own objects: the threadprivate data of Lib_CEA_data and the critical section of
# Lib_Chemistry_wdot are ignored without it (a consumer compiling the sources gave them its own flag)
if (USE_OPENMP)
  if (NOT TARGET OpenMP::OpenMP_Fortran)
    find_package(OpenMP REQUIRED COMPONENTS Fortran)
  endif()
  target_link_libraries(FLINT PUBLIC OpenMP::OpenMP_Fortran)
endif()

# -------- Build FLINT test suite --------
# ctest in the build directory runs the tests registered in src/test/CMakeLists.txt
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
  enable_testing()
  add_subdirectory("src/test")
endif()