# Verification

FLINT is verified against [Cantera](https://cantera.org), the reference: its results on the batch reactor cases and
the equilibrium sweeps of the test suite, and its speed on the batch reactor.

<div class="grid cards" markdown>

-   :material-check-decagram:{ .lg .middle } __Accuracy__

    ---

    Temperature histories of constant-volume batch reactors, and equilibrium temperatures, of FLINT and Cantera for
    the mechanisms of the database.

    [:octicons-arrow-right-24: View the comparisons](accuracy.md)

-   :material-speedometer:{ .lg .middle } __Performance__

    ---

    Execution times of the batch reactor with the dedicated routines, the general procedure and Cantera.

    [:octicons-arrow-right-24: View the timings](performance.md)

</div>

The test suite checks the same cases at every pull request, against Cantera references stored in the repository and
without Cantera: see [Testing](../development/testing.md).
