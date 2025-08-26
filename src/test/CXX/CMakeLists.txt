
set(CMAKE_CXX_STANDARD 17)

set(exe "test-batchCXX")
add_executable(${exe} "test-batchCXX.cpp")

#Use one Sundials installation
#include_directories("${canteraDIR}/ext/sundials/include/")
include_directories("${CMAKE_SOURCE_DIR}/lib/OSlo/sundials-install/include/")
include_directories("${canteraD}/include/")

if(APPLE)
    # On macOS use Accelerate framework
    find_library(ACCELERATE_FRAMEWORK Accelerate)
    if(ACCELERATE_FRAMEWORK)
        message(STATUS "Found Accelerate: ${ACCELERATE_FRAMEWORK}")
        target_link_libraries(${exe} PRIVATE ${ACCELERATE_FRAMEWORK})
    endif()
else()
    # On Linux, look for BLAS/LAPACK
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)

    if(BLAS_FOUND AND LAPACK_FOUND)
        target_link_libraries(${exe} PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    else()
        message(FATAL_ERROR "BLAS/LAPACK not found")
    endif()
endif()

target_link_libraries(${exe} PRIVATE FLINT )#fmt::fmt Threads::Threads)

set_target_properties(${exe} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/test")