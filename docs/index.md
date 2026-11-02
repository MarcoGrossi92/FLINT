<p align="center">
  <img src="flint-logo-concept-1.png" alt="FLINT logo" width="600">
</p>

# FLINT — Fortran Library for INtegrated Thermochemistry

FLINT is a lightweight, modular Fortran library for thermochemistry in reactive flow simulations. It gives a CFD
solver thermodynamic and transport properties, chemical source terms and chemical equilibrium, through an interface
that keeps the thermochemistry apart from the flow solver.

- Ideal-gas thermodynamics, and real-fluid thermodynamics on pressure–enthalpy tables
- Transport properties: viscosity, thermal conductivity, diffusion
- Finite-rate chemical kinetics, with dedicated routines for the mechanisms of its database and a general procedure
  for any other mechanism
- Chemical equilibrium, with the algorithm of NASA CEA
- 0-D constant-volume reactors
- An optional Cantera backend, against which FLINT is verified

## Where to Start

<div class="grid cards" markdown>

-   :material-rocket-launch:{ .lg .middle } __Getting Started__

    ---

    Install FLINT, then build and run a first program that uses it.

    [:octicons-arrow-right-24: Get started](getting-started/index.md)

-   :material-book-open-variant:{ .lg .middle } __User Guide__

    ---

    Embed FLINT in a solver: the workflow, the API, the input files and the mechanisms.

    [:octicons-arrow-right-24: Read the User Guide](user/index.md)

-   :material-check-decagram:{ .lg .middle } __Verification__

    ---

    FLINT's results and speed, compared with Cantera's.

    [:octicons-arrow-right-24: View the verification](verification/index.md)

-   :material-function-variant:{ .lg .middle } __Theoretical Guide__

    ---

    The models of FLINT and their equations.

    [:octicons-arrow-right-24: Read the theory](theory/index.md)

-   :material-file-tree:{ .lg .middle } __Developer Guide__

    ---

    The code, the tests, the release process, and how to contribute.

    [:octicons-arrow-right-24: Read the Developer Guide](development/index.md)

-   :material-information-outline:{ .lg .middle } __About__

    ---

    Acknowledgements and license.

    [:octicons-arrow-right-24: About FLINT](about/index.md)

</div>
