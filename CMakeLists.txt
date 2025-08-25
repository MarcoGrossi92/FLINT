cmake_minimum_required(VERSION 3.23)

# -------- Compiler overrides (for macOS) --------
# if (UNIX AND APPLE)
#     set(CMAKE_CXX_COMPILER "/opt/homebrew/bin/g++-14")
#     set(CMAKE_C_COMPILER   "/opt/homebrew/bin/gcc-14")
# endif ()

# -------- Say my name! --------
set(lib "cantIF")

project (cantIF LANGUAGES C CXX Fortran)

set(CMAKE_VERBOSE_MAKEFILE FALSE)

# -------- Git submodules --------
# ORION
if(MASTER STREQUAL None)
  set ( orionD ${CMAKE_CURRENT_SOURCE_DIR}/lib/ORION/ )
elseif(MASTER STREQUAL hydra)
  set ( orionD $ENV{HYDRADIR}/lib/ORION/ "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/ORION/")
endif()
add_subdirectory(${orionD})

# OSlo
if(MASTER STREQUAL None)
  set ( osloD ${CMAKE_CURRENT_SOURCE_DIR}/lib/OSlo/ )
elseif(MASTER STREQUAL hydra)
  set ( osloD $ENV{HYDRADIR}/lib/OSlo/ "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/Oslo/")
endif()
add_subdirectory(${osloD})

# Cantera
OPTION(USE_CANTERA "Use Canera" ON)
if (USE_CANTERA)
  set(canteraD "${CMAKE_CURRENT_SOURCE_DIR}/lib/cantera")
  # Check if cloned
  if(NOT EXISTS "${canteraD}/README.rst")
    message(STATUS "Cantera not found. Attempting to initialize...")
    execute_process(COMMAND git submodule update --init --recursive lib/cantera
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
            sundials_include="${osloD}"/sundials-install/include
            sundials_libdir="${osloD}"/sundials-install/lib
            boost_inc_dir=/opt/homebrew/opt/boost/include
        INSTALL_COMMAND scons -C ${canteraD} install
    )
  else()
    message(STATUS "Cantera already installed.")
  endif()
  link_directories("${canteraD}/build/lib")
  add_compile_definitions(CANTERA)
endif(USE_CANTERA)

# -------- Fortran directories --------
include(GNUInstallDirs)
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# -------- Compile  --------
add_subdirectory(src)