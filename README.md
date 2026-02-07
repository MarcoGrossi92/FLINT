# FLINT – Fortran Library for INtegrated Thermochemistry

A lightweight, modular Fortran library for accessing thermodynamic, transport, and chemical source term data with seamless integration into CFD solvers and reactive flow applications.

FLINT provides a clean interface layer that separates thermochemistry concerns from CFD complexity, enabling efficient and maintainable reactive flow simulations. Like the stone it's named after, FLINT is small but powerful — it provides the spark that ignites complex thermochemical simulations.

## Key Advantages

- **Decoupled Thermochemistry** – Focus on pure thermochemistry without CFD complexity
- **Modular Architecture** – Build and test chemical models independently
- **Optional Cantera Backend** – Leverage advanced kinetic models when needed
- **Production-Ready** – Built-in library of verified chemical models
- **Developer-Friendly** – Automated YAML-to-Fortran (YTF) source code generation

## Features

- Species thermodynamics (ideal gas models)
- Transport coefficients (viscosity, thermal conductivity, diffusion)
- Chemical source terms with finite-rate kinetic reactions
- Chemical source terms with equilibrium assumption

---

## Requirements

### System Requirements
- **Fortran Compiler:** Intel or GNU
- **C compiler**: For SUNDIALS support
- **CMake:** >= 3.23
- **Python:** 3.8+ (for utilities like YTF, optional)

### Optional Requirements
- **Cantera:** 2.5+ (for advanced kinetic models and validation)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## External Dependencies

This project is composed of multiple components maintained as separate
repositories and included here as git submodules:

- **ORION**  
  Input/output library for structured multi-block scientific data  
  https://github.com/MarcoGrossi92/ORION

- **OSlo**  
  Fortran-based ODE solver toolbox  
  https://github.com/MarcoGrossi92/OSlo

- **Cantera**  
  Chemical kinetics, thermodynamics, and transport tool suite  
  https://github.com/Cantera/cantera

## Acknowledgments

This project incorporates concepts and algorithms from:

- **NASA CEA (Chemical Equilibrium with Applications)**  
  Reference implementation for the chemical equilibrium solver  
  Apache License 2.0  
  https://github.com/nasa/cea

---

## Installation

The installation process provides:
- A Fortran static library linkable from external projects
- Test binaries

**WARNING!**  
Cantera installation and linking to the project may result a bit tricky. It is strongly recommended to try to install Cantera outside of this project and, if succesful, build FLINT with Cantera support.

### Quick Start

Clone the repository and run the installation script:

```bash
git clone https://github.com/MarcoGrossi92/FLINT.git
cd FLINT
./install.sh build --compiler=gnu
```

By default, Cantera support is not **enabled**. If requested, Cantera will be cloned, installed, and linked to the project (see warning above). 

### Build Options

The `install.sh` script provides several build options:

```bash
./install.sh [GLOBAL_OPTIONS] COMMAND [COMMAND_OPTIONS]
```

**Global Options:**
- `-v, --verbose`: Enable verbose output

**Commands:**

- **build**: Perform a full build
  - `--compilers=<name>`: Set compiler suite (intel, gnu)
  - `--include-orion=<path>`: Set external ORION path
  - `--include-oslo=<path>`: Set external OSlo path
  - `--use-cantera`: Use Cantera (Sundials required)
  - `--use-sundials`: Use Sundials (via OSlo)
  - `--use-tecio`: Use TecIO (via ORION)


### Manual CMake Build

If you prefer to use CMake directly:

```bash
# With Cantera and Sundials support
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DUSE_CANTERA=ON \
         -DUSE_SUNDIALS=ON

# Parallel build
make -j8
```

## Usage

### Basic Usage

```fortran
program simple_thermo
  use FLINT_Load_ThermoTransport
  use FLINT_Load_Chemistry
  implicit none
  integer :: err

  ! Load thermodynamic data
  err = read_idealgas_thermo(folder='WD/INPUT')
  if (err /= 0) print *, "Error loading thermo:", err

  ! Load transport properties
  err = read_idealgas_transport(folder='WD/INPUT')
  if (err /= 0) print *, "Error loading transport:", err

  ! Load chemical mechanism
  err = read_chemistry(folder='WD/INPUT')
  if (err /= 0) print *, "Error loading mechanism:", err

end program simple_thermo
```

### Example Programs

See the [src/test/](src/test/) directory for:
- Basic thermodynamic calculations
- Chemical source term computation
- 0D batch reactor
- Equilibrium chemistry

---

## Project Structure

```
FLINT/
├── src/              # Main source code
│   ├── lib/         # Library modules
│   └── test/        # Fortran tests
├── lib/             # External libraries (Cantera, OSlo, etc.)
├── utils/           # Utilities (YTF, etc.)
├── test/            # Integration and mechanism tests
├── cmake/           # CMake modules
├── docs/            # Documentation
└── build/           # Build directory (generated)
```

## Contributing

1. **Fork** the repository
2. **Create** a feature branch: `git checkout -b feature/my-feature`
3. **Commit** with clear messages
4. **Add and run tests** for new functionality
5. **Push** and create a pull request

## Quick Links

- [Documentation](docs/) – Full API documentation
- [Examples](src/test/) – Example programs and test cases
- [Issues](../../issues) – Report bugs or request features

## References & Resources

- **Cantera Documentation:** [cantera.org](https://cantera.org)

---
