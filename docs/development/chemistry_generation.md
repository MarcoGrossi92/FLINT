# Chemistry Code Generation

FLINT includes a chemistry code generation tool that converts chemical mechanism files into optimized, standalone source code for reaction kinetics.

The generation process provides a Fortran chemistry source subroutine compatible with FLINT starting from a Cantera YAML mechanism. Given a mechanism name, it parses the corresponding YAML file using Cantera’s Python API and writes a Fortran module implementing the species source terms for that mechanism.

This procedure has the same objective of [KinetiX](https://github.com/bogdandanciu/KinetiX).

---

## Functionality

The script:

* Reads a Cantera YAML mechanism (`<mech>.yaml`)
* Extracts:

  * Species list
  * Reaction list
  * Stoichiometric coefficients
  * Third-body efficiencies
  * Falloff parameters (Lindemann and Troe)
* Generates a Fortran module:
  ```
  <mech>_mod
  ```
  containing a subroutine:
  ```
  subroutine <mech>(roi, temp, omegadot)
  ```

The resulting subroutine:

* Converts species densities to molar concentrations
* Evaluates forward and backward reaction rates
* Handles:

  * Elementary Arrhenius reactions
  * Three-body Arrhenius reactions
  * Lindemann falloff reactions
  * Troe falloff reactions
* Computes species mass production rates

---

## Supported Reaction Types

The following Cantera reaction types are supported:

* `Arrhenius`
* `three-body-Arrhenius`
* `falloff-Lindemann`
* `falloff-Troe`

* Reactions are treated as reversible by default
* To enforce irreversibility, the backward rate must be zero in the mechanism data
* Species efficiencies default to unity if not explicitly specified

---

## Input and Output

**Input**

* YAML mechanism file:
  ```
  <mech>.yaml
  ```
* Command-line argument:
  ```
  python YTF.py <mech>
  ```

**Output**

* Fortran source file:
  ```
  <mech>.f90
  ```

---

## Generated Fortran Interface

```fortran
subroutine <mech>(roi, temp, omegadot)
  real(8), intent(inout) :: roi(ns)      ! species densities [kg/m^3]
  real(8), intent(in)    :: temp         ! temperature [K]
  real(8), intent(out)   :: omegadot(ns) ! species source terms [kg/m^3/s]
end subroutine
```

The generated code depends on the following FLINT modules:

* `FLINT_Lib_Thermodynamic`
* `FLINT_Lib_Chemistry_data`
* `FLINT_Lib_Chemistry_falloff`

---

## Limitations and Caveats

* No custom reaction ordering is implemented
* The script assumes:

  * Consistent species ordering between FLINT and Cantera
  * Ideal-gas thermodynamics
* The script has been lightly tested only on mechanisms already included in FLINT
* Error handling is minimal
* The generated code has not been formally verified for all edge cases

⚠️ **Use at your own risk**: this tool is intended primarily for developer use and internal mechanism generation.

---

## Development History

* Initial prototype by [Marco Fabiani](https://github.com/marcofabiani)
* Refactored and adapted for the FLINT workspace by [Marco Grossi](https://github.com/MarcoGrossi92)

---
