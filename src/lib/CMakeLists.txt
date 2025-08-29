# CMakeLists.txt for the FLINT library
file(GLOB LIBSRC *.f90 Lib_ChemMech/*.f90)
add_library(${lib} OBJECT ${LIBSRC})
target_include_directories(${lib} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
