# cmake/LinkDependencies.cmake

function(link_common_dependencies target_name)
  # Optimized math libraries
  if(APPLE)
    find_library(ACCELERATE_FRAMEWORK Accelerate)
    if(ACCELERATE_FRAMEWORK)
      message(STATUS "Found Accelerate: ${ACCELERATE_FRAMEWORK}")
      target_link_libraries(${target_name} PRIVATE ${ACCELERATE_FRAMEWORK})
    endif()
  else()
    find_package(BLAS)
    find_package(LAPACK)
    if(BLAS_FOUND AND LAPACK_FOUND)
      target_link_libraries(${target_name} PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    endif()
  endif()

  # Cantera
  if (USE_CANTERA)
    target_link_libraries(${target_name} PRIVATE cantera cantera_fortran)
  endif()

  # OSlo
  target_link_libraries(${target_name} PRIVATE OSlo)

endfunction()
