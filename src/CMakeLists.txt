
#find_package(Threads REQUIRED)

include_directories("${canteraDIR}/build/src/fortran/")

add_subdirectory("lib")
add_subdirectory("test/Fortran")
if (USE_CANTERA)
  add_subdirectory("test/CXX")
endif()
