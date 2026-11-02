# Performance

The execution time of the constant-volume batch reactor was measured for the mechanisms of the database, with four
ways of computing the chemical source terms:

- Cantera (C++): Cantera's own reactor
- FLINT-Cantera: FLINT's integrator with the source terms of Cantera's Fortran interface
- FLINT-explicit: FLINT's integrator with the dedicated routine of the mechanism
- FLINT-general: FLINT's integrator with the general procedure

The dedicated (explicit) routines are the fastest for every mechanism tested. [Mechanism Routines](../user/chemistry_routines.md)
describes how FLINT selects a dedicated routine or the general procedure.

Benchmark environment:

- MacBook Air (Apple M1, 8 cores: 4P + 4E)
- 16 GB RAM
- macOS
- Homebrew GCC 15.2.0
- Apple clang 17.0.0 (clang-1700.6.3.2)
- Single-threaded execution

<figure>
  --8<-- "docs/verification/images/barplot.svg"
  <figcaption>Normalized execution times</figcaption>
</figure>

## Reproducing the timings

With a build with Cantera (see [Testing](../development/testing.md#building-with-cantera-and-regenerating-the-figures)),
from `test/batch`:

```bash
../../bin/test-cantera/test-batchCXX performance   # Cantera                          -> comp-batch-cantera.dat
../../bin/test-cantera/test-batchF performance     # FLINT-Cantera, -explicit, -general -> comp-batch-<backend>.dat
python3 batch-performance.py                       # docs/verification/images/barplot.svg
```

The chart leaves out the WD and Pelucchi cases.
