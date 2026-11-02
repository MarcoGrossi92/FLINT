# Acknowledgements

FLINT is built upon several projects.

## ORION

**I/O Library for Fortran**

- **Repository:** [github.com/MarcoGrossi92/ORION](https://github.com/MarcoGrossi92/ORION)
- **License:** GPL v3.0

ORION provides built-in functions to read files in different formats. FLINT uses ORION to import thermodynamic, transport, and chemistry tabulated data.

## OSLO

**ODE solvers for Fortran**

- **Repository:** [github.com/MarcoGrossi92/OSLO](https://github.com/MarcoGrossi92/OSLO)
- **License:** GPL v3.0

OSLO provides a state-of-the-art Fortran framework to solve ODE systems. FLINT uses OSLO to build code drivers to test chemistry integration.

## CEA

**Chemical Equilibrium Solver**

- **Repository:** [github.com/nasa/cea](https://github.com/nasa/cea)
- **License:** Apache 2.0

Developed by NASA, CEA is the reference code for chemical equilibrium computations.

## Cantera

**Chemical kinetics, thermodynamics and transport**

- **Repository:** [github.com/Cantera/cantera](https://github.com/Cantera/cantera)
- **License:** BSD-3-Clause

Cantera is the reference against which FLINT's results are verified, through its Fortran interface (optional).

## SUNDIALS

**Suite of nonlinear and differential/algebraic equation solvers**

- **Repository:** [github.com/LLNL/sundials](https://github.com/LLNL/sundials)
- **License:** BSD-3-Clause

OSLO builds SUNDIALS and offers its integrators (optional); Cantera is built on the same SUNDIALS.

## Documentation Tools

### Zensical

**Static Site Generator**

- **Website:** [zensical.org](https://zensical.org)
- **License:** MIT
- **Contribution:** Documentation site generation and theme

Zensical, by the team of Material for MkDocs, turns FLINT's documentation into a searchable website. It reads the
configuration of MkDocs (`mkdocs.yml`), with which the documentation was first written.

## Institutional Support

While FLINT is an independent open-source project, it has benefited from:

- Academic research environments
- Access to various CFD datasets for validation

## License Compliance

FLINT respects all licenses of dependencies and foundations:

- **GPL v3.0** - FLINT, ORION and OSLO
- **Apache 2.0** - CEA
- **BSD-3-Clause** - Cantera, SUNDIALS
- **Boost Software License 1.0** - the Boost headers used to build Cantera
- **Tecplot License** - TecIO library

See the [License](license.md) page for FLINT's full license text.

---
