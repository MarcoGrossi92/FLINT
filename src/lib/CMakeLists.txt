file(GLOB LIBSRC *.f90 Lib_ChemMech/*.f90)
add_library(${lib} ${LIBSRC})
if(APPLE)
    # On macOS use Accelerate framework
    find_library(ACCELERATE_FRAMEWORK Accelerate)
    if(ACCELERATE_FRAMEWORK)
        message(STATUS "Found Accelerate: ${ACCELERATE_FRAMEWORK}")
        target_link_libraries(${lib} PRIVATE ${ACCELERATE_FRAMEWORK})
    endif()
else()
    # On Linux, look for BLAS/LAPACK
    find_package(BLAS)
    find_package(LAPACK)

    if(BLAS_FOUND AND LAPACK_FOUND)
        target_link_libraries(${lib} PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    endif()
endif()
target_link_libraries(${lib} PRIVATE ORION OSlo)
if (USE_TECIO)
  target_link_libraries(${lib} PRIVATE ${LINKLIBS})
endif()
if (USE_CANTERA)
#   if(APPLE)
#     target_link_libraries(${lib} PRIVATE c++)
#   elseif(UNIX)
#     target_link_libraries(${lib} PRIVATE stdc++)
#   endif()
  target_link_libraries(${lib} PRIVATE cantera cantera_fortran)
endif()

target_include_directories(${lib}
  INTERFACE
      $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)
