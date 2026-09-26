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
the routine silently read the wrong species. The check needs the loaded tables: the usual order is
`read_idealgas_thermo`, `read_chemistry`, then `Assign_Mechanism` (MOSE and Q2D). If
`Assign_Mechanism` is called with a hooked name **before** the tables are loaded, it selects the
routine as before, prints a `[WARNING]` (standard output and error unit) and the check runs at the
first call of `chemistry_source` or `chemistry_jacobian`, with the same refusal on a mismatch (and
an `[ERROR]` if the tables are still not loaded then). That first call may come from several threads
of a parallel region of the host program. With FLINT compiled with the OpenMP flags of the host, the
check runs in a critical section: one thread checks while the others wait, then all call the routine.
With FLINT compiled without OpenMP the critical section is only a comment: threads that make the first
call together may each run the check (a refusal is then printed once per thread before the run stops),
and every thread still ends on the selected routine, because the procedure pointers are assigned only
from local copies of it. (The CMake option `USE_OPENMP` of FLINT adds no OpenMP compile flag.) A name
that is not hooked is not checked (it falls back to the general procedure). `test-contract` exercises the rules; to add a hooked mechanism, add its record to the
JSON and re-run the generator (`python3 utils/mechanism_contract.py`). The Fortran text of the
checker is kept in the generator too: `python3 utils/mechanism_contract.py --check` (run by
`test-contract`) exits with 1 and prints the difference when the committed module is not the
output of the generator.

### The contract file

`src/lib/Lib_ChemMech/mechanism_contract.json` is FLINT's published description of its hooked routines:
a table writer can compare the species and reactions of the folder it writes for a hooked name with what
the routine expects, without reading FLINT's sources. It is maintained with FLINT:
`python3 utils/mechanism_contract.py --fingerprints` writes `n_reactions`, `fingerprint` and
`zeroes_omegadot` from the routine sources; `python3 utils/mechanism_contract.py` writes the checker
module from the records; `--check` (run by `test-contract`) fails when a fingerprint or
`zeroes_omegadot` differs from the sources, when the table counts of a generated routine differ from
its `nrc`, or when the module is not the output of the records.

Top level: `schema` (one-line summary of the format), `origin` (how the records are made), `cases` (one
record per hooked name; the key is the mechanism name of line 1 of `chemistry-info.txt`). A record:

| field | content |
|-------|---------|
| `routine`, `file` | the subroutine selected by `Assign_Mechanism` and its source file in `src/lib/Lib_ChemMech` |
| `kind` | `generated` (written by `utils/YTF.py`, one `! reac n. <n>: <equation>` comment per reaction) or `hand-written` |
| `ns` | number of species slots |
| `nrc` | reactions per table type: `arrhenius` (every type whose name contains neither Troe nor Lindemann), `troe`, `lindemann` |
| `species` | the slots in routine order: `slot` (from 1), `name`, `composition` (element: count) |
| `tables_used` | the routine reads the rate tables |
| `zeroes_omegadot` | the routine sets the whole `omegadot` to zero first, so species after its slots get a zero source |
| `jacobian` | the analytical Jacobian routine, or `null` |
| `source` | where the record comes from |
| `n_reactions` | generated routines only: the number of reactions |
| `fingerprint` | generated routines only: the structure of every reaction (below) |

The fingerprint is `{"M": [...], "reactions": {...}, "omegadot": {...}}`:

- `M`: the distinct third-body efficiency sets of the routine, in order of first use; a set is a list of
  `[efficiency, [slots]]` pairs; a slot with efficiency 0 is not listed (`M=sum(coi(1:12))` gives
  `[[1.0, [1, 2, ..., 12]]]`).
- `reactions`: the key is the reaction number n of the `! reac n.` comment (a string); the value is
  `{"eq": <equation of the comment>, "f": <forward rate>, "b": <reverse rate>}`, and a rate is
  `{"tab": [<type>, <index>], "fac": {"<slot>": <count>}, "usesM": <bool>, "Mi": <index>}`: `tab` is the
  table the rate reads (`["arrhenius", i]`: column i of the Arrhenius-type tables through `f_kf`/`f_kb`;
  `["troe", j]` or `["lindemann", j]`: falloff table j through `k(1)`/`k(2)`); `fac` counts the factors
  `coi(slot)` of the rate (the generated routines write integer powers as repeated factors); `usesM` says
  whether the rate is multiplied by the third-body concentration `M`; `Mi` (from 0) points to the set of
  `M` given by the `M=` line of the reaction (the third body of `usesM`, or the bath gas of a falloff rate),
  and is absent when the reaction has no `M=` line.
- `omegadot`: the key is the slot s (a string); the value `{"<n>": <nu'' - nu'>}` gives the net
  stoichiometric coefficient of slot s in reaction n as written in `omegadot(s)`; reactions with a zero
  net coefficient are not listed.

## Fallback to the general procedure and strict mode

A mechanism name that is not hooked in `Assign_Mechanism` selects the data-driven `general`
procedure: FLINT prints `[WARNING] Explicit procedure for <name> not found, defaulting to the
general procedure` on standard output **and** on the error unit (standard error), so that a
solver log that captures only one of the two channels still records the fallback. Note that
`general` applies the mass-action law with the real stoichiometric coefficients, as Cantera does
(`H2 + 0.5 O2 <=> H2O`: forward `[H2] [O2]^0.5`; fractional products enter the reverse rate the same
way); explicit yaml `orders:` reach it only through the `Reaction orders` block that ends
`chemistry-info.txt` (see *Native input*). A file without the block comes from an older table writer:
`general` then takes the reactant coefficients as orders (the exponents of a block with no rows) and,
once per load, prints on standard output and on the error unit `[WARNING] FLINT: chemistry-info.txt has
no 'Reaction orders' block (the file comes from an older table writer): ... regenerate the chemistry
tables ...`; for a mechanism with explicit orders it is then a different model. A hooked name does not
warn: the compiled routines do not read the block. Versions before `test-stoich` rounded the coefficients of the Arrhenius-type reactions
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

## Rate laws coded in the routines

A compiled routine fixes its species slots, its stoichiometry and the exponents of the concentrations in
its source; the INPUT folder supplies the species data and the rate constants (column i of the rate
tables is reaction i of the routine). A mechanism file with the same reactions but other orders, or
without explicit orders, is a different model: under the general procedure (or Cantera) its rates follow
the mass-action law of the file, not the exponents below. The generated routines raise the concentrations
to the stoichiometric coefficients of the file they were generated from (`utils/YTF.py` implements no
custom orders); their reaction structure is published in the contract file (see *The contract file*
below). Concentrations in kmol/m3, `kf`/`kb` = the forward/backward columns of the table of the reaction.

- **`CoronettiC4H6` selects `Coronetti`** (`coronetti.f90`); slots O2, C4H6, H2O, CO, CO2, H2, O, H, OH;
  negative partial densities of the argument `roi` are set to zero:
    1. C4H6 + 2 O2 → 4 CO + 3 H2: kf1 [C4H6]^0.5 [O2]^1.25, zero when [C4H6] or [O2] < 1e-10;
    2. C4H6 + 4 H2O → 4 CO + 7 H2: kf2 [C4H6] [H2O] (the mass-action law of this equation gives
       [C4H6] [H2O]^4);
    3. CO + H2O ⇌ CO2 + H2: kf3 [CO] [H2O] − kb3 [CO2] [H2];
    4. H2 + ½ O2 ⇌ H2O: kf4 [H2]^0.25 [O2]^1.5 − kb4 [H2O] [O2] [H2]^-0.75 (the mass-action law gives
       [H2] [O2]^0.5 and [H2O]; forward minus reverse orders equal the stoichiometric coefficients, so
       kf4/kb4 is the equilibrium constant of the step); the forward term is zero when [H2] or [O2] <
       1e-10, the reverse term when [H2O], [O2] or [H2] < 1e-10;
    5. O2 ⇌ 2 O: kf5 [O2] − kb5 [O]^2;
    6. H2O ⇌ OH + H: kf6 [H2O] − kb6 [OH] [H].
- **`JLR-Nasuti` selects `JLR`** (`JLR.f90`); slots O2, CH4, H2O, CO, CO2, H2, H, O, OH; concentrations
  below 1e-10 are taken as zero:
    1. CH4 + ½ O2 → CO + 2 H2: kf1 [CH4]^0.5 [O2]^1.25;
    2. CH4 + H2O → CO + 3 H2: kf2 [CH4] [H2O];
    3. CO + H2O ⇌ CO2 + H2, 5. O2 ⇌ 2 O, 6. H2O ⇌ H + OH, 7. OH + H2 ⇌ H + H2O: mass-action law;
    4. H2 + ½ O2 ⇌ H2O: kf4 [H2]^0.25 [O2]^1.5 − kb4 [H2O] [O2] [H2]^-0.75, forward term only when
       [H2] < 1e-10.
- **`Frassoldati` selects `Frassoldati`** (`JLR.f90`); the slots and steps 1-6 of `JLR`, with kf1
  [CH4]^0.5 [O2]^1.3 and kf4 [H2]^0.3 [O2]^1.55 − kb4 [H2O] [O2] [H2]^-0.75: the forward minus reverse
  orders of step 4 are 1.05 (H2) and 0.55 (O2) instead of the stoichiometric 1 and 0.5, so kf4/kb4 is not
  the equilibrium constant of the step. Only the name `Frassoldati` selects this routine: a mechanism named
  otherwise (e.g. `JLR-Frassoldati` of `test-orders`) runs on the general procedure, with the orders of its
  `Reaction orders` block.
- **`Frolov` selects `Frolov`** (`global-H2.f90`); slots O2, H2O, H2 (+ inert species): the routine is the
  model. The progress rate of 2 H2 + O2 → 2 H2O is coded as q = 0.5 × 8e11 × (p/101325)^-1.15 × [H2]^2
  [O2] × exp(−10000/T), with p = Σ roi Ri T from the state and concentrations below 1e-12 taken as zero;
  the rate tables are not read (the INPUT folder supplies the species, their molecular weights and
  thermodynamics), so a pressure dependence written in the mechanism file (e.g. as a PLOG rate) plays no
  part. The H2O slot accepts a calibrated water species under another name with the composition of H2O.
- **`Frolov_nopressure` selects `Frolov_nopressure`** (`global-H2.f90`); slots O2, H2O, H2, N2 (the N2
  row is zero): 2 H2 + O2 ⇌ 2 H2O with kf [H2]^2 [O2] − kb [H2O]^2 from reaction 1 of the tables, no
  pressure factor; analytical Jacobian `Frolov_nopressure_jac`.
- **`Nassini` selects `Nassini_4`** (`global-H2.f90`); slots O2, H2O, H2 (+ inert species),
  concentrations below 1e-12 taken as zero; two irreversible reactions: 1. H2 + ½ O2 → H2O: kf1 [H2] [O2];
  2. H2O → H2 + ½ O2: kf2 [H2O], the backward step (the kb columns are not read). These exponents hold
  whatever orders the mechanism file gives.
- **`ONERA-7` selects `ONERA_7`** (`ONERA-7.f90`, generated); slots O2, H2O, H2, H, O, OH, N2: 14
  Arrhenius-type reactions, each of the 7 steps written as two irreversible reactions (1. H2 + O2 ⇒ 2 OH,
  2. 2 OH ⇒ H2 + O2, ...; third body in 11-14); analytical Jacobian `ONERA_7_jac`. An INPUT folder with the
  7 steps as reversible reactions has 7 Arrhenius-type reactions: the contract check refuses it under the
  name `ONERA-7` (reaction count), and it runs on the general procedure under a name that is not hooked.
- **`FFCMy-12` selects `FFCMy_12`** and **`SanDiego` selects `sandiego20161214`** (generated): 13 slots
  (12 reacting species and N2, which enters as a third body only) and 38 reactions (34 Arrhenius-type including three-body ones, 3
  falloff-Troe, 1 falloff-Lindemann); 57 slots and 268 reactions (245 Arrhenius-type, 23 falloff-Troe).

The rate tables carry Arrhenius-type, falloff-Troe and falloff-Lindemann reactions only: a reaction type
without tables of its own (e.g. falloff-SRI) is counted with the Arrhenius-type ones and `read_chemistry`
refuses the folder (`ios = 4`: fewer zones in `chemistry-Arrhenius.dat` than Arrhenius-type reactions).
