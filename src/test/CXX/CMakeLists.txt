
set(CMAKE_CXX_STANDARD 17)

set(exe "test-batchCXX")
add_executable(${exe} "test-batchCXX.cpp")

#Use one Sundials installation
include_directories("${CMAKE_SOURCE_DIR}/lib/OSlo/sundials-install/include/")
include_directories("${canteraD}/include/")

link_FLINT_dependencies(${exe})
target_link_libraries(${exe} PRIVATE ORION)

set_target_properties(${exe} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/test")