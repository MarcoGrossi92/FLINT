
set(CMAKE_CXX_STANDARD 17)

set(exe "test-batchCXX")

#find_package(Threads REQUIRED)
#find_package(fmt REQUIRED)
find_library(ACCELERATE_FRAMEWORK Accelerate)
include_directories("${canteraDIR}/include/")

add_executable(${exe} "test-batchCXX.cpp")
target_link_libraries(${exe} cantera ${ACCELERATE_FRAMEWORK} )#fmt::fmt Threads::Threads)

set_target_properties(${exe} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/test")