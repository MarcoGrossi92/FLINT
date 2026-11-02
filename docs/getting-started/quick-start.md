# Quick Start

This page builds and runs a first program that uses FLINT. It loads a mechanism of FLINT's database, then computes
the properties of a gas mixture and its chemical source terms. It needs CMake 3.23 or newer, a Fortran compiler, Git
and the network (see [Prerequisites](installation.md#prerequisites)); FLINT does not have to be installed first.

## 1. The Project

The project is a folder with three files:

```text
quickstart/
├── CMakeLists.txt
├── cmake/SemVerDeps.cmake
└── main.f90
```

`cmake/SemVerDeps.cmake` is the CMake module with which a project requires FLINT by release. Copy it from FLINT's
repository:

```bash
mkdir -p quickstart/cmake && cd quickstart
curl -fL -o cmake/SemVerDeps.cmake https://raw.githubusercontent.com/MarcoGrossi92/FLINT/v2.0.0/cmake/SemVerDeps.cmake
```

`CMakeLists.txt` requires FLINT 2.0.0 or a later 2.x, builds the program, and copies a mechanism of FLINT's database,
Westbrook-Dryer (`WD`: CH4, O2, CO2, H2O, CO), to the folder `INPUT/` from which the program reads:

```cmake
cmake_minimum_required(VERSION 3.23)
project(quickstart LANGUAGES Fortran)

# FLINT 2.0.0 or a later 2.x, fetched at configure time with ORION and OSLO
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/SemVerDeps.cmake)
semver_dependency(FLINT GIT_REPOSITORY https://github.com/MarcoGrossi92/FLINT.git VERSION 2.0
                  GIT_SUBMODULES "")
FetchContent_MakeAvailable(FLINT)

add_executable(quickstart main.f90)
target_link_libraries(quickstart PRIVATE FLINT)

# The input of the program: the Westbrook-Dryer mechanism of FLINT's database, in INPUT/ of the build folder
file(COPY ${flint_SOURCE_DIR}/database/WD/ DESTINATION ${CMAKE_BINARY_DIR}/INPUT)
```

`main.f90`:

```fortran
program quickstart
  use FLINT_Load_ThermoTransport, only: read_idealgas_thermo
  use FLINT_Load_Chemistry,       only: read_chemistry
  use FLINT_Lib_Thermodynamic,    only: ns, species_names, Ri_tab, f_cp, f_gamma, f_ss
  use FLINT_Lib_Chemistry_wdot,   only: Assign_Mechanism, chemistry_source
  implicit none
  character(len=64) :: mech_name
  real(8), allocatable :: Y(:), rhoi(:), wdot(:)
  real(8) :: p, T, R, rho
  integer :: i

  ! Load the species, their thermodynamic tables and the reactions of the mechanism in INPUT/
  if (read_idealgas_thermo('INPUT') /= 0) error stop 'cannot load the thermodynamic data'
  if (read_chemistry('INPUT', mech_name) /= 0) error stop 'cannot load the chemistry'
  call Assign_Mechanism(mech_name)            ! the dedicated routine of the mechanism

  ! A CH4/O2 mixture at 1 atm and 1500 K, given to FLINT by its partial densities
  allocate(Y(ns), rhoi(ns), wdot(ns))
  Y = 0d0
  do i = 1, ns
    select case (trim(species_names(i)))
    case ('CH4'); Y(i) = 0.2d0
    case ('O2');  Y(i) = 0.8d0
    end select
  end do
  p = 101325d0
  T = 1500d0
  R = sum(Y*Ri_tab)                           ! gas constant of the mixture, J/(kg K)
  rho = p/(R*T)
  rhoi = rho*Y

  ! Properties of the mixture, and its chemical source terms
  print '(a,f10.2,a)', 'cp    = ', f_cp(rhoi, T, rho), ' J/(kg K)'
  print '(a,f10.4)',   'gamma = ', f_gamma(rhoi, p, rho, R)
  print '(a,f10.2,a)', 'c     = ', f_ss(rhoi, p, rho, R), ' m/s'
  call chemistry_source(rhoi, T, wdot)
  do i = 1, ns
    print '(a6,es12.4,a)', trim(species_names(i)), wdot(i), ' kg/(m3 s)'
  end do
end program quickstart
```

## 2. Build and Run

```bash
cmake -S . -B build -DUSE_TECIO=OFF -DUSE_SUNDIALS=OFF
cmake --build build
cd build && ./quickstart
```

The first configure fetches FLINT, ORION and OSLO, and writes the releases it selected in `deps.lock`. The options
leave out TecIO and SUNDIALS, which this program does not need. The program prints:

```text
cp    =    2046.64 J/(kg K)
gamma =     1.1795
c     =     742.42 m/s
   CH4 -2.1476E+02 kg/(m3 s)
    O2 -6.4252E+02 kg/(m3 s)
   CO2  0.0000E+00 kg/(m3 s)
   H2O  4.8232E+02 kg/(m3 s)
    CO  3.7496E+02 kg/(m3 s)
```

## 3. What the Program Does

- `read_idealgas_thermo` reads the species and their molecular weights (`phase.txt`) and their thermodynamic tables
  (`thermo.dat`). `read_chemistry` reads the reactions (`chemistry-info.txt`) and their rate tables, and returns
  the name of the mechanism, `WD`, with which `Assign_Mechanism` selects its dedicated routine.
- FLINT describes a mixture by its partial densities `rhoi` (kg/m³), here from the mass fractions `Y`, the pressure
  and the temperature.
- `f_cp`, `f_gamma` and `f_ss` return the specific heat, the ratio of specific heats and the speed of sound of the
  mixture. `chemistry_source` returns the net mass production rate of each species: CH4 and O2 are consumed, H2O and
  CO produced, and the rates add up to zero.

## Next Steps

- **[Using the Library](../user/using.md)**: the workflow of each feature: thermodynamics, chemistry, equilibrium
- **[Examples](../user/examples.md)**: a batch reactor and an equilibrium calculation
- **[API Reference](../user/api.md)**: every procedure, with its arguments and return codes
- **[Using FLINT in a CMake Project](../user/cmake.md)**: the static library or the sources, and the compile flags
- **[Available Mechanisms](../user/database.md)**: the other mechanisms of `database/`. To try one, replace `WD` in
  `CMakeLists.txt`, and the mixture in `main.f90` with species of that mechanism
