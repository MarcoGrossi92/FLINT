
#find_package(Threads REQUIRED)

include_directories("${canteraDIR}/build/src/fortran/")

add_subdirectory("lib")
add_subdirectory("test/cantIF")
