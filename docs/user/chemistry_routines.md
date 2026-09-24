# Chemistry Code Generation

FLINT supports both general-purpose chemistry solver and mechanism-specific dedicated routines. The first one is a general subroutine that solves finite-rate chemistry using input data to define all necessary . On the other hand, dedicated routines are source files optimized for a specific chemical mechanism, hard-coding all required info directly in the source files.

The employment of mechanism-specific routines provides significant performance improvements and it is particularly important for production runs.

The generation of a new mechanism routine is always recommended. However, this operation requires the modification of the API, therefore it is suggested for developers or advanced users.

A python program is shipped with FLINT to automatically generate source files starting from yaml files. See the proper [Development](../development/chemistry_generation.md) section for more info.

!!! warning "API modifications required" 
    This procedure generates new FLINT source files and requires rebuilding the library. It is intended for advanced users and developers.

## Performance Comparison

A performance comparison among different methods of integration was performed by testing different chemical mechanisms on 0D batch reactor.

The results clearly show how explicit (hard-coded) routines outperform the other alternatives for each tested mechanism.

Benchmark environment:

- MacBook Air (Apple M1, 8 cores: 4P + 4E)
- 16 GB RAM
- macOS
- Homebrew GCC 15.2.0
- Apple clang 17.0.0 (clang-1700.6.3.2)
- Single-threaded execution

<figure>
  {% include "user/images/barplot.svg" %}
  <figcaption>Normalized execution times</figcaption>
</figure>


## Mechanism contract check

A mechanism name hooked in `Assign_Mechanism` selects a compiled routine whose species slots and
reaction tables are fixed in the source. `Assign_Mechanism`
compares the data loaded by `read_idealgas_thermo` / `read_chemistry` with the expectation record
of the routine (`src/lib/Lib_ChemMech/mechanism_contract.json`, turned into
`src/lib/Lib_Chemistry_contract.f90` by `utils/mechanism_contract.py`):

- the routine species must be the **first** `ns_r` loaded species, in the routine order
  (prefix semantics); species loaded after them are accepted and stay inert (the RHS and the
  Jacobian zero the source term before calling the routine);
- a slot matches by **elemental composition** (from `composition.txt`), and by **name** too for a
  generated routine or where two slots of the routine share a composition; a hand-written routine
  accepts a calibrated species under another name with the same composition (e.g. `H2ONassini`
  for `H2O`); without `composition.txt` (INPUT folders written before the table writer produced
  it) the molecular weight of `phase.txt` stands in for the composition (tolerance 0.05 kg/kmol)
  **and** the loaded name must begin with the slot name, case-insensitive (`H2ONassini` for `H2O`):
  the weight alone cannot tell CO from N2 or C2H4 (28.010, 28.014, 28.054 kg/kmol);
- the number of reactions per table type (Arrhenius, falloff-Troe, falloff-Lindemann) must match.

On a mismatch both lists are printed and the run stops (`error stop`): with the previous versions
the routine silently read the wrong species. The tables must therefore be loaded **before**
`Assign_Mechanism`. A name that is not hooked is not checked (it falls back to the general
procedure). `test-contract` exercises the rules; to add a hooked mechanism, add its record to the
JSON and re-run the generator.

## Fallback to the general procedure and strict mode

A mechanism name that is not hooked in `Assign_Mechanism` selects the data-driven `general`
procedure: FLINT prints `[WARNING] Explicit procedure for <name> not found, defaulting to the
general procedure` on standard output **and** on the error unit (standard error), so that a
solver log that captures only one of the two channels still records the fallback. Note that
`general` raises the concentrations to the integer stoichiometric coefficients (`nint`): a
mechanism with fractional reaction orders is a different model under `general`.

Strict mode turns the fallback into a refusal (`[ERROR] FLINT Assign_Mechanism: mechanism <name>
is not hooked and strict mode is on`, exit status 1): set the module flag
`FLINT_strict_mechanism = .true.` (module `FLINT_Lib_Chemistry_wdot`) before `Assign_Mechanism`,
or the environment variable `FLINT_STRICT_MECHANISM=1` (also `true`, `yes`, `on`). Default: off.
