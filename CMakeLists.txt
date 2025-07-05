cmake_minimum_required(VERSION 3.10)

# -------- Compiler overrides (for macOS) --------
# if (UNIX AND APPLE)
#     set(CMAKE_CXX_COMPILER "/opt/homebrew/bin/g++-14")
#     set(CMAKE_C_COMPILER   "/opt/homebrew/bin/gcc-14")
# endif ()

# -------- Say my name! --------
set(lib "cantIF")

project (cantIF LANGUAGES CXX Fortran)

set(CMAKE_VERBOSE_MAKEFILE FALSE)

# -------- Git submodules --------
# ORION
add_subdirectory("lib/ORION")

# OSlo
add_subdirectory("lib/OSlo")

# Cantera
set(canteraDIR "~/Codici/cantera/")
link_directories("${canteraDIR}/build/lib")

# -------- Fortran directories --------
include(GNUInstallDirs)
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# -------- Compile  --------
add_subdirectory(src)