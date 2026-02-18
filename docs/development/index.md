# Developer Guide

This section provides resources for contributors and developers working on FLINT.

<div class="grid cards" markdown>

-   :material-file-tree:{ .lg .middle } __Code Structure__

    ---

    Overview of the repository layout, module organization,
    and build system.

    [:octicons-arrow-right-24: Explore structure](structure.md)

-   :material-atom:{ .lg .middle } __Generating New Mechanisms__

    ---

    How FLINT generates dedicated source files for new
    chemical mechanisms.

    [:octicons-arrow-right-24: Extend chemistry routines](chemistry_generation.md)

-   :material-test-tube:{ .lg .middle } __Testing Infrastructure__

    ---

    Regression tests, validation drivers, and numerical
    verification strategy.

    [:octicons-arrow-right-24: View testing framework](testing.md)

-   :material-source-branch:{ .lg .middle } __Contributing__

    ---

    Guidelines for pull requests, coding standards,
    and repository workflow.

    [:octicons-arrow-right-24: Contribution guidelines](contributing.md)

</div>

The Developer Guide is intended for contributors and advanced users who want to:

- Understand FLINT’s internal architecture
- Extend the present capabilities
- Add new mechanisms or dedicated routines
- Contribute improvements or fixes

FLINT is designed with modularity and reproducibility in mind.
Its structure separates:

- Source code library
- Testing and verification drivers

For basic usage and API access, see the [User Guide](../user/index.md).
