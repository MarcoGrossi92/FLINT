# User Guide

The User Guide describes how to configure and use FLINT in practical applications.


<div class="grid cards" markdown>

-   :material-cog-outline:{ .lg .middle } __Using the Library__

    ---

    How to use FLINT API.

    [:octicons-arrow-right-24: Embed FLINT in your code](using.md)

-   :material-flask-outline:{ .lg .middle } __Examples__

    ---

    A batch reactor and an equilibrium calculation, in pseudo-code.

    [:octicons-arrow-right-24: Follow the examples](examples.md)

-   :material-api:{ .lg .middle } __API Reference__

    ---

    Every procedure and module variable, with its arguments,
    units and return codes.

    [:octicons-arrow-right-24: Look up the API](api.md)

-   :material-language-fortran:{ .lg .middle } __Using FLINT in a CMake Project__

    ---

    Fetch FLINT with CMake, and link the static library or compile the sources.

    [:octicons-arrow-right-24: Build FLINT in your project](cmake.md)

-   :material-file-document-outline:{ .lg .middle } __Input Format__

    ---

    How to define thermodynamic data, chemistry, and simulation parameters.
    Includes YAML and native formats.

    [:octicons-arrow-right-24: Configure input files](input/index.md)

-   :material-database:{ .lg .middle } __Chemistry Database__

    ---

    The mechanisms of the database, and how FLINT selects their dedicated routines or the general procedure.

    [:octicons-arrow-right-24: Explore mechanisms](database.md)

</div>

This section is intended for:

* Developers embedding FLINT in a CFD or reacting-flow solver
* Users configuring thermodynamic and chemical mechanisms
* Engineers running simulations with dedicated or general chemistry routines

If you are new to FLINT, start from [Getting Started](../getting-started/index.md).

The User Guide covers:

* API basics for all features, with examples
* Building FLINT in a CMake project
* Input formats (YAML and native)
* Available chemical mechanisms, and their dedicated routines

It does not describe:

* The mathematical background (see [Theoretical Guide](../theory/index.md))
* Internal architecture, and how to generate the routine of a new mechanism (see [Developer Guide](../development/index.md))
