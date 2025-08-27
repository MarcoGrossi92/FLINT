
#find_package(Threads REQUIRED)

if (USE_CANTERA)
  include_directories("${canteraD}/build/src/fortran/")
endif()

add_subdirectory("lib")
add_subdirectory("test/Fortran")
if (USE_CANTERA)
  add_subdirectory("test/CXX")
endif()
