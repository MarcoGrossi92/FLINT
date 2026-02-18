# Acknowledgements

FLINT is built upon several projects.

## ORION

**I/O Library for Fortran**

- **Repository:** [github.com/MarcoGrossi92/ORION](https://github.com/MarcoGrossi92/ORION)
- **License:** GPL v3.0

ORION provides built-in functions to read files in different formats. FLINT uses ORION to import thermodynamic, transport, and chemistry tabulated data.

## OSlo

**ODE solvers for Fortran**

- **Repository:** [github.com/MarcoGrossi92/OSlo](https://github.com/MarcoGrossi92/OSlo)
- **License:** GPL v3.0

OSlo provides a state-of-the-art Fortran framework to solve ODE systems. FLINT uses OSlo to build code drivers to test chemistry integration.

## CEA

**Chemical Equilibrium Solver**

- **Repository:** [github.com/MarcoGrossi92/ORION](https://github.com/MarcoGrossi92/ORION)
- **License:** Apache 2.0

Developed by NASA, CEA is the reference code for chemical equilibrium computations.

## Documentation Tools

### MkDocs

**Static Site Generator**

- **Website:** [mkdocs.org](https://www.mkdocs.org)
- **License:** BSD-2-Clause
- **Contribution:** Documentation site generation

MkDocs transforms FLINT's documentation into a beautiful, searchable website.

### Material for MkDocs

**Modern Documentation Theme**

- **Website:** [squidfunk.github.io/mkdocs-material](https://squidfunk.github.io/mkdocs-material/)
- **License:** MIT
- **Contribution:** Professional documentation theme with advanced features

Material for MkDocs provides the sleek, modern interface for FLINT's documentation.

## Institutional Support

While FLINT is an independent open-source project, it has benefited from:

- Academic research environments
- Access to various CFD datasets for validation

## License Compliance

FLINT respects all licenses of dependencies and foundations:

- **GPL v3.0** - FLINT's license, compatible with OFF
- **BSD Licenses** - Lib_VTK_IO, NumPy, CMake
- **MIT License** - Material for MkDocs
- **Tecplot License** - TecIO library

See the [License](license.md) page for FLINT's full license text.

---
