# User Guide

The User Guide describes how to configure and use FLINT in practical applications.


<div class="grid cards" markdown>

-   :material-cog-outline:{ .lg .middle } __Using the Library__

    ---

    How to use FLINT API.

    [:octicons-arrow-right-24: Embed FLINT in your code](using.md)


-   :material-file-document-outline:{ .lg .middle } __Input Format__

    ---

    How to define thermodynamic data, chemistry, and simulation parameters.
    Includes YAML and native formats.

    [:octicons-arrow-right-24: Configure input files](input/index.md)

-   :material-database:{ .lg .middle } __Chemistry Database__

    ---

    Available mechanisms and how to extend FLINT with new dedicated routines.

    [:octicons-arrow-right-24: Explore mechanisms](database.md)

</div>

This section is intended for:

* Developers embedding FLINT in a CFD or reacting-flow solver
* Users configuring thermodynamic and chemical mechanisms
* Engineers running simulations with dedicated or general chemistry routines

If you are new to FLINT, start from [Getting Started](../getting-started/index.md).

## Scope of This Section

The User Guide covers:

* API basics for all features
* Input formats (YAML and native)
* Available chemical mechanisms and how to generate optimized mechanism-specific routines

It does not describe:

* The mathematical background (see [Theoretical Guide](../theory/index.md))
* Internal architecture (see [Developer Guide](../development/index.md))
