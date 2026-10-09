# Quick Start

This section verifies that FLINT is correctly installed by running one of the provided test binaries.

The test suite exercises the main FLINT capabilities and provides immediate feedback on numerical correctness and performance.

---

## Running the Tests

After a successful build, the test executables are in `bin/test` and registered with CTest. From the
build directory, run the quick tier of the suite (about 30 s):

```bash
cd build
ctest -L quick
```

or every test with `ctest -j8`. Each test prints its checks and a verdict, and exits with code 1 on failure; CTest summarizes them.

A single driver can also be run by hand from its area folder in `test/`, for example the universal gas constant
and the ideal-gas pressure against Cantera values:

```bash
cd test/thermo
../../bin/test/test-runiv
```

A successful run indicates that:

* FLINT is correctly linked
* Thermodynamic data are loaded properly
* The numerical kernels are functioning as expected

See [Testing](../development/testing.md) for the areas, the labels and the reference data.

---

## Notes on Cantera

If Cantera is available at build time, some tests will automatically compare
FLINT results against Cantera reference solutions.

If Cantera is not available:

* All FLINT native routines remain fully functional
* Test programs will run using FLINT-only paths
* This is the standard and recommended configuration for production use

---

## Next Steps

Once FLINT is built and verified, you may want to explore:

* **User Guide**  
  * Running simulations
  * Input formats
  * Chemistry databases

* **Examples**  
  * Verification and validation results

* **Developer Guide**  
  * Testing infrastructure
  * Generation of dedicated chemistry kernels
  * Extending FLINT

---
