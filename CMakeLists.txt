cmake_minimum_required(VERSION 3.23)

# -------- Compiler overrides (for macOS) --------
# if (UNIX AND APPLE)
#     set(CMAKE_CXX_COMPILER "/opt/homebrew/bin/g++-15")
#     set(CMAKE_C_COMPILER   "/opt/homebrew/bin/gcc-15")
# endif ()

# -------- Say my name! --------
set(lib "FLINT")

project (FLINT LANGUAGES C CXX Fortran)

if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
  set(CMAKE_VERBOSE_MAKEFILE FALSE)

  # Search path for cmake modules
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

  # This INCLUDE statement executes code that sets the compile flags for DEBUG,
  # RELEASE, and TESTING.
  INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 

endif()

# -------- Dependencies linking --------
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(LinkDependencies)

# -------- Path fallbacks (standalone or when parent does not define them) --------
if (NOT DEFINED ORION_PATH)
  set(ORION_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/ORION/")
endif()
if (NOT DEFINED OSLO_PATH)
  set(OSLO_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/OSLO/")
endif()

# -------- Helper: initialize a git submodule from the correct git root --------
# Works whether the submodule belongs to this project or a parent project.
macro(ensure_submodule_initialized SUBMOD_PATH)
  if (NOT EXISTS "${SUBMOD_PATH}/CMakeLists.txt")
    message(STATUS "Submodule at ${SUBMOD_PATH} not initialized. Attempting git submodule update --init...")
    get_filename_component(_sm_parent "${SUBMOD_PATH}" DIRECTORY)
    execute_process(
      COMMAND git rev-parse --show-toplevel
      WORKING_DIRECTORY "${_sm_parent}"
      OUTPUT_VARIABLE _sm_git_root
      OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE _sm_git_result)
    if(_sm_git_result EQUAL 0)
      file(RELATIVE_PATH _sm_rel "${_sm_git_root}" "${SUBMOD_PATH}")
      execute_process(
        COMMAND git submodule update --init "${_sm_rel}"
        WORKING_DIRECTORY "${_sm_git_root}"
        RESULT_VARIABLE _sm_init_result)
      if(NOT _sm_init_result EQUAL 0)
        message(FATAL_ERROR "Failed to initialize submodule: ${SUBMOD_PATH}")
      endif()
    else()
      message(FATAL_ERROR "Could not find git root for submodule: ${SUBMOD_PATH}")
    endif()
  endif()
endmacro()

# -------- Git submodules --------

if (USE_CANTERA)
  set ( USE_SUNDIALS ON )
endif()

# OSLO (needed early for SUNDIALS before Cantera)
if (NOT TARGET OSLO)
  ensure_submodule_initialized("${OSLO_PATH}")
  add_subdirectory("${OSLO_PATH}" "${CMAKE_BINARY_DIR}/lib/OSLO")
endif()

# ORION
if (NOT TARGET ORION)
  ensure_submodule_initialized("${ORION_PATH}")
  add_subdirectory("${ORION_PATH}" "${CMAKE_BINARY_DIR}/lib/ORION")
endif()

# Cantera
# Cantera is optional and only used by FLINT
OPTION(USE_CANTERA "Use Cantera" ON)
if (USE_CANTERA)
  # Find SUNDIALS (created by OSLO subdirectory)
  set(SUNDIALS_INSTALL_DIR "${OSLO_PATH}/sundials-install")
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
  # Check if already installed
  set(CANTERA_CONF "${canteraD}/config.log")
  if(NOT EXISTS "${CANTERA_CONF}")
    message(STATUS "Cantera not installed, building and installing it...")
    find_package(Boost)
    if (NOT Boost_FOUND)
      message(FATAL_ERROR "Boost is required to build Cantera. Please install Boost and try again.")
    endif()
    if(APPLE)
      set (sundials_libdir ${OSLO_PATH}/sundials-install/lib)
    elseif(UNIX)
      set (sundials_libdir ${OSLO_PATH}/sundials-install/lib64)
    endif()
    include(ExternalProject)
    ExternalProject_Add(
        cantera_world
        PREFIX cantera_build
        SOURCE_DIR ${canteraD}
        BINARY_DIR ${CMAKE_BINARY_DIR}/cantera_build
        CONFIGURE_COMMAND ""  # Cantera does not use configure step
        BUILD_COMMAND scons -C ${canteraD} build
            prefix=${CMAKE_BINARY_DIR}/cantera_install
            f90_interface=y
            FORTRAN=gfortran
            system_sundials=y
            sundials_include=${OSLO_PATH}/sundials-install/include
            sundials_libdir=${sundials_libdir}
            boost_inc_dir=${Boost_INCLUDE_DIR}
            system_eigen=n
        INSTALL_COMMAND scons -C ${canteraD} install
    )
  else()
    message(STATUS "Cantera already installed.")
  endif()
  link_directories("${canteraD}/build/lib")
  include_directories("${canteraD}/build/src/fortran/")
  add_compile_definitions(CANTERA)
endif(USE_CANTERA)

# -------- Fortran directories --------
include(GNUInstallDirs)
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# -------- Define FLINT sources --------
file(GLOB LIBSRC src/lib/*.f90 src/lib/Lib_ChemMech/*.f90)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
  set(FLINT_SRC ${LIBSRC})
else()
  set(FLINT_SRC ${LIBSRC} PARENT_SCOPE)
endif()

# -------- Build FLINT test suite --------
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
  add_subdirectory("src/test/Fortran")
  if (USE_CANTERA)
    add_subdirectory("src/test/CXX")
  endif()
endif()