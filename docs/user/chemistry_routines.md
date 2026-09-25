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
  Jacobian zero the source term before calling the routine); inert also as third bodies: a compiled
  routine sums its third-body concentration over its own slots only (e.g. `troyes.f90`, `ONERA-7.f90`),
  the general procedure over every loaded species with the efficiencies of `chemistry-info.txt`, so an
  appended collider with a non-zero efficiency counts in `general` only (an appended species with
  efficiency 0 in every reaction counts in neither, and both agree);
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
`general` applies the mass-action law with the real stoichiometric coefficients, as Cantera does
(`H2 + 0.5 O2 <=> H2O`: forward `[H2] [O2]^0.5`; fractional products enter the reverse rate the same
way); explicit yaml `orders:` reach it only through the optional `Reaction orders` block of the INPUT
folder (see *Native input*): without the block a mechanism with explicit orders is a different model
under `general`. Versions before `test-stoich` rounded the coefficients of the Arrhenius-type reactions
to the nearest integer (`[O2]^1` for `0.5 O2`): an INPUT folder with fractional coefficients and no
block gives different rates since then.

Strict mode turns the fallback into a refusal (`[ERROR] FLINT Assign_Mechanism: mechanism <name>
is not hooked and strict mode is on`, exit status 1): set the module flag
`FLINT_strict_mechanism = .true.` (module `FLINT_Lib_Chemistry_wdot`) before `Assign_Mechanism`,
or the environment variable `FLINT_STRICT_MECHANISM=1` (also `true`, `yes`, `on`, in any
case; `0`, `false`, `no`, `off` or an empty value leave it off; any other value is reported with a
`[WARNING]` on standard output and on the error unit and ignored). Default: off.

## Zero-concentration convention for negative reaction orders

One rule for every site of FLINT where a species enters a rate with a negative order: at zero
concentration the rate of that step is zero, as in Cantera 3.0.1 (its forward rate of progress is 0
when a species with a non-zero order has zero concentration, verified by execution), and the negative
power is never evaluated there (`0**(-0.25)` is +Infinity and traps under `-fpe0` / `-ffpe-trap=zero`).
The sites: `pow_order` of the general procedure (`Reaction orders` block), step 3 of `WD-Andersen`
([O2]^-0.25), the reverse term of H2 + 1/2 O2 <-> H2O of `Coronetti` ([H2]^-0.75, evaluated only above
the routine's 1e-10 kmol/m3 threshold, where it is used) and step 4 of the JLR family ([H2]^-0.75, the
forward-only branch below 1e-10 kmol/m3, unchanged). The thresholds of Coronetti and JLR are the
historical behaviour of those routines and are kept: for every state above them the results are
bit-identical to the previous version (`test-andersen` checks the convention, `test-orders` the helper).

## Temperature range of the tables

`rhs_native` and `jac_native` return the bail-out value (`F = -1`, zero Jacobian) when the
temperature is outside the thermo tables **or** outside the rate tables (`chemistry-*.dat`):
row T of every table is the value at T kelvin and a rate table must not be read below its first
row. The rate tables loaded by `read_chemistry` cover the thermo range (see the temperature grid of
the tables; a wider rate table is accepted), so the second guard is a defence for tables set by
another path. The tables of `database/WD` start at 1 K with a 1 K step; every table of
`database/TSR-Rich-31` (thermo and rates) starts at 500 K.

## Analytical Jacobian availability

`chemistry_jacobian` (module `FLINT_Lib_Chemistry_wdot`) is a null pointer by default: only the
`ONERA-7` (`ONERA_7_jac`) and `Frolov_nopressure` (`Frolov_nopressure_jac`) routines set it. For
every other mechanism `set_analytical_jacobian(.true.)` prints `[JAC] analytical Jacobian
requested but the active mechanism has none; falling back to finite differences` and the
integrator builds the Jacobian by finite differences.
