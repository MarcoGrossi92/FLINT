
set(CMAKE_CXX_STANDARD 17)

set(exe "test-batchCXX")

#find_package(Threads REQUIRED)
#find_package(fmt REQUIRED)
find_library(ACCELERATE_FRAMEWORK Accelerate)
include_directories("${canteraD}/include/")
#include_directories("${canteraDIR}/ext/sundials/include/")
include_directories("${CMAKE_SOURCE_DIR}/lib/OSlo/sundials-install/include/")

add_executable(${exe} "test-batchCXX.cpp")
target_link_libraries(${exe} PRIVATE cantIF ${ACCELERATE_FRAMEWORK} )#fmt::fmt Threads::Threads)

set_target_properties(${exe} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/test")