cmake_minimum_required(VERSION 3.10)
project (exe)

set(CMAKE_VERBOSE_MAKEFILE OFF)
enable_language(Fortran CXX)

find_package(Threads REQUIRED)

set(canteraDIR /home/marcogrossi/Codici/cantera)

include_directories(${canteraDIR}/build/src/fortran)
link_directories(${canteraDIR}/build/lib)

add_subdirectory(src/test)
