# cmake/LinkDependencies.cmake

function(link_FLINT_dependencies target_name)
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
    # cantera_fortran before cantera, which it calls: GNU ld searches a static library only where it is on the line
    target_link_libraries(${target_name} PRIVATE cantera_fortran cantera)
    # After the Cantera build (cantera_world, when this configure builds Cantera): the target uses its Fortran modules
    # and libraries, which a parallel build would otherwise look for too early
    if (TARGET cantera_world)
      add_dependencies(${target_name} cantera_world)
    endif()
    if (APPLE)
      target_link_libraries(${target_name} PRIVATE c++)
    elseif (UNIX)
      target_link_libraries(${target_name} PRIVATE stdc++)
    endif()
  endif()

  # OSLO
  target_link_libraries(${target_name} PRIVATE OSLO)

endfunction()
