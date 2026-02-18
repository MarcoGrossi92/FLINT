# Quick Start

This section verifies that FLINT is correctly installed by running one of the provided test binaries.

The test suite exercises the main FLINT capabilities and provides immediate feedback on numerical correctness and performance.

---

## Running a Test Binary

After a successful build, test executables are located in the `bin/test`
directory.

As a first check, run a thermodynamic verification test:

```bash
cd test
./../bin/test/test-thermo
```

This program compares FLINT native thermodynamic routines against reference
implementations (if enabled) and prints execution times and numerical errors.

A successful run indicates that:

* FLINT is correctly linked
* Thermodynamic data are loaded properly
* The numerical kernels are functioning as expected

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
