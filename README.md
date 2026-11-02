# FLINT — Fortran Library for INtegrated Thermochemistry

FLINT is a lightweight, modular Fortran library for thermochemistry in reactive flow simulations.

It provides a clean interface to thermodynamic, transport, and chemical source term models, designed to integrate seamlessly into CFD solvers while keeping thermochemistry decoupled from CFD complexity.

Small but powerful — like the stone it’s named after — FLINT provides the spark for complex thermochemical simulations.

## Why FLINT?

- **Decoupled thermochemistry**  
  Develop and test chemical models independently of CFD solvers.

- **Modular & extensible**  
  Plug in different chemistry, transport, or equilibrium models with minimal coupling.

- **Optional Cantera backend**  
  Access advanced kinetic models and validation when needed.

- **Production-ready**  
  Ships with a library of verified chemical mechanisms.

- **Developer-friendly**  
  Automated YAML-to-Fortran (YTF) code generation for mechanisms and models.

## Capabilities

- Ideal-gas thermodynamics
- Real-fluid thermodynamics (pressure–enthalpy tables)
- Transport properties (viscosity, thermal conductivity, diffusion)
- Finite-rate chemical kinetics
- Equilibrium chemistry
- 0D reactors and source-term evaluation

## A Taste of FLINT

Thermochemistry setup in FLINT is explicit, modular, and solver-agnostic

```fortran
program taste_of_flint
  use FLINT_Load_ThermoTransport
  use FLINT_Load_Chemistry
  implicit none

  integer :: err

  ! Load ideal-gas thermodynamics
  err = read_idealgas_thermo(folder='WD/INPUT')
  if (err /= 0) stop "Error loading thermodynamics"

  ! Load chemical mechanism
  err = read_chemistry(folder='WD/INPUT')
  if (err /= 0) stop "Error loading chemistry"

end program taste_of_flint
```

## Quick Start

```bash
git clone https://github.com/MarcoGrossi92/FLINT.git
cd FLINT
./install.sh build --compilers=gnu
```

The first configure fetches ORION and OSLO from GitHub. Cantera support is **optional** and disabled by default:
`--use-cantera` builds Cantera 3.2.0 from the submodule `lib/cantera`.

## Requirements

* Fortran compiler (GNU or Intel)
* CMake ≥ 3.23
* Git, and network access at the first configure
* Python ≥ 3.8 (optional, for utilities)
* C compiler (for SUNDIALS support)

Optional:

* SCons, Boost and a C++ compiler, to build Cantera (`--use-cantera`)

## Documentation

📘 **Full documentation, theory, and examples:**
[https://marcogrossi92.github.io/FLINT/](https://marcogrossi92.github.io/FLINT/)

See the docs for:

* Installation details
* [API reference](https://marcogrossi92.github.io/FLINT/user/api/)
* Examples and test cases
* YAML-to-Fortran workflow

## Getting Help

Ask questions, report bugs and request features on the [issue tracker](https://github.com/MarcoGrossi92/FLINT/issues),
which has a template for each.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). The changes of each release are in [CHANGELOG.md](CHANGELOG.md); FLINT
follows [Semantic Versioning](https://marcogrossi92.github.io/FLINT/development/versioning/).

## License

GNU General Public License v3.0
See [LICENSE](LICENSE).

## Acknowledgments

FLINT incorporates concepts and algorithms from:

* **NASA CEA — Chemical Equilibrium with Applications**
  [https://github.com/nasa/cea](https://github.com/nasa/cea)
