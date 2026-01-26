# Native File Format Specification

Complete reference for FLINT native file formats (non-Cantera builds).

## Table of Contents

1. [Overview](#overview)
2. [File Structure](#native-file-structure)
3. [Detailed Format Specifications](#detailed-format-specifications)
   - [phase.txt](#1-phasetxt--species-definitions)
   - [input.ini](#2-inputini--configuration)
   - [chemistry-info.txt](#3-chemistry-infotxt--mechanism-metadata)
   - [chemistry-Arrhenius.dat](#4-chemistry-arrheniusdat--reaction-rate-tables)
   - [chemistry-Troe.dat](#5-chemistry-troedat--fall-off-correction-tables)
   - [thermo.dat](#6-thermodat--thermodynamic-property-tables)
4. [Complete Example](#complete-example-coria-mechanism)
5. [Usage Workflow](#native-input-path-usage)

---

## Overview

When FLINT is compiled **without Cantera support**, input data is loaded from a set of pre-generated native files. This path is designed for **high-performance production use** where pre-computed tables eliminate the need for runtime polynomial evaluation or mechanism parsing.

---

## Native File Structure

The complete input consists of up to **seven files** in the same directory:

```
INPUT/
  ├── phase.txt               # Species list
  ├── chemistry-info.txt      # Chemistry - Metadata & stoichiometry
  ├── chemistry-Arrhenius.dat # Chemistry - Reaction rates
  ├── chemistry-Troe.dat      # Chemistry - Fall-off parameters (optional)
  ├── chemistry-Lindemann.dat # Chemistry - Fall-off parameters (optional)
  ├── transport.dat           # Transport properties (optional)
  └── thermo.dat              # Thermodynamic properties
```

**File Sizes (CORIA 17-species / 44-reaction example):**
- phase.txt: ~250 bytes (trivial)
- input.ini: ~40 bytes (trivial)
- chemistry-info.txt: ~15 KB (metadata)
- chemistry-Arrhenius.dat: ~43 MB (15000 temperature points × 44 reactions)
- chemistry-Troe.dat: ~8 MB (fall-off reactions only)
- thermo.dat: ~15 MB (15000 temperature points × 17 species × 5 properties)

**Total: ~81 MB** for a medium-sized mechanism. Pre-tabulation trades memory for speed.

---

## Detailed Format Specifications

### 1. `phase.txt` – Species Definitions

**Purpose:** Define species names and molecular weights.

**Format:** ASCII text, one species per line.

```
ideal-gas phase
<species_name> <molecular_weight_g_mol>
<species_name> <molecular_weight_g_mol>
...
```

**Example (CORIA mechanism):**

```
ideal-gas phase
H2      2.016000
H       1.008000
O      15.999000
O2     31.998000
OH     17.007000
H2O    18.015000
HO2    33.006000
CH3    15.035000
CH4    16.043000
CO     28.010000
CO2    44.009000
H2O2   34.014000
HCO    29.018000
CH2O   30.026000
CH3O   31.034000
C2H6   30.070000
CH3O2  47.033000
```

**Parsing Notes:**
- Line 1 is always `ideal-gas phase` (header)
- Subsequent lines: space-separated `<name>` and `<MW>`
- MW values are floats in g/mol
- Species **order must match** `chemistry-info.txt` and `.dat` files
- Typical precision: 6 decimal places

---

### 2. `chemistry-info.txt` – Mechanism Metadata

**Purpose:** Provide resume of the chemical mechanism: species/reaction counts, reaction types, and stoichiometry.

**Format:** ASCII text with three sections.

```
<mechanism_name>
N.ro species = <count>
N.ro reactions = <count>

General loop info. They are not used if the mechanism is exlplicitly defined

Reaction type
<index> <type>
<index> <type>
...

Reaction definition
<rxn_id> <species_name> <coeff_reactant> <coeff_product> <coeff_thirdbody>
<rxn_id> <species_name> <coeff_reactant> <coeff_product> <coeff_thirdbody>
...
```

**Example (CORIA, header section):**

```
CORIA-CNRS
N.ro species = 17
N.ro reactions = 44

General loop info. They are not used if the mechanism is exlplicitly defined

Reaction type
1 Arrhenius
2 Arrhenius
3 Arrhenius
4 Arrhenius
5 three-body-Arrhenius
6 Arrhenius
7 falloff-Troe
...
44 Arrhenius

Reaction definition
1 H2 1.0 0.0 0.0
1 H 0.0 1.0 0.0
1 O 1.0 0.0 0.0
1 O2 0.0 0.0 0.0
...
```

**Section 1: Mechanism Header**

```
<mechanism_name>
N.ro species = <n_species>
N.ro reactions = <n_reactions>
```

- `<mechanism_name>`: Arbitrary string (e.g., "CORIA-CNRS", "GRI-MECH3")
- `n_species`: Total number of distinct species
- `n_reactions`: Total number of chemical reactions

**Section 2: Reaction Type List**

```
Reaction type
<reaction_index> <reaction_type>
...
```

Reaction types:
- `Arrhenius` – Standard elementary reaction (2-body, 3-body forms handled by coefficients)
- `three-body-Arrhenius` – Three-body reaction (includes third-body collision partner)
- `falloff-Troe` – Pressure-dependent reaction with Troe fall-off correction

**Section 3: Reaction Definition (Stoichiometry Matrix)**

```
<reaction_id> <species_name> <n_reactant> <n_product> <n_thirdbody>
```

For each reaction `<reaction_id>`, there are `n_species` lines (one per species):

| Column | Type | Description |
|--------|------|-------------|
| `<reaction_id>` | int | Reaction index (1–44 in CORIA) |
| `<species_name>` | string | Species name (must match `phase.txt`) |
| `<n_reactant>` | float | Stoichiometric coefficient as reactant (0 if product/inert) |
| `<n_product>` | float | Stoichiometric coefficient as product (0 if reactant/inert) |
| `<n_thirdbody>` | float | Third-body efficiency (0 if not involved, 1.0 if collision partner) |

**Example:** Reaction 5 (three-body-Arrhenius):
```
5 H2 0.0 0.0 0.0    # H2 is neither reactant nor third body
5 H 0.0 1.0 0.0     # H is product
5 O 1.0 0.0 0.0     # O is reactant
5 O2 0.0 0.0 1.0    # O2 is third body (efficiency 1.0)
5 OH 0.0 1.0 0.0    # OH is product
5 H2O 0.0 0.0 1.0   # H2O is third body
...
5 M 0.0 0.0 1.0     # Generic third body (efficiency 1.0)
```

**Parsing Notes:**
- 846 lines total in CORIA (lines 1–4 header + 44 reactions × 18 species ≈ 844 lines)
- Species appear in **same order** as in `phase.txt`
- Third-body coefficients: `0.0` (not involved), `1.0` (default), other floats (enhanced efficiency)

---

### 3. `chemistry-Arrhenius.dat` – Reaction Rate Tables

**Purpose:** Pre-computed forward and backward Arrhenius rates for all reactions across temperature range.

**Format:** Tecplot ASCII (binary-like structure, not human-readable).

```
TITLE = "Arrhenius data"
VARIABLES = "Temperature", "kf", "kb"
ZONE T=Reaction1
I=15000, F=POINT
<T_1>  <kf_1>  <kb_1>
<T_2>  <kf_2>  <kb_2>
...
<T_15000>  <kf_15000>  <kb_15000>

ZONE T=Reaction2
I=15000, F=POINT
<T_1>  <kf_1>  <kb_1>
...
```

**Structure:**

| Component | Description |
|-----------|-------------|
| TITLE | File description (fixed: "Arrhenius data") |
| VARIABLES | Column headers (fixed: "Temperature", "kf", "kb") |
| ZONE | Zone header (one per reaction; `T=Reaction<N>`) |
| I=N | Number of temperature points (e.g., 15000) |
| F=POINT | Data format (always POINT = single values per line) |
| Data rows | `Temperature  kf  kb` (space-separated floats) |

**Example (first reaction, first 6 points):**

```
ZONE T=Reaction1
I=15000, F=POINT
1.0           0.00000000000000000000E+00    0.00000000000000000000E+00
2.0           0.00000000000000000000E+00    0.00000000000000000000E+00
3.0           0.00000000000000000000E+00    0.00000000000000000000E+00
4.0           0.00000000000000000000E+00    0.00000000000000000000E+00
5.0           4.31541883433268038949E-272    4.26390547350824394621E-194
6.0           4.65732315911476184290E-226    6.07657342283880940567E-161
...
```

**Data Format:**

- **Temperature**: 1.0 to 15000.0 K (in 1 K increments)
- **kf**: Forward rate constant (scientific notation, units depend on reaction order)
  - For 2-body Arrhenius: cm³/mol/s
  - For 3-body: cm⁶/mol²/s
- **kb**: Backward rate constant (same units as kf, from reverse reaction)
- Precision: 20 significant decimal places (double precision)

**Temperature Indexing:**

The file contains **one zone per reaction**, with all reactions using the **same temperature grid**:
- Zone 1 (Reaction 1): T = 1 to 15000 K
- Zone 2 (Reaction 2): T = 1 to 15000 K
- ...
- Zone 44 (Reaction 44): T = 1 to 15000 K

To access rate for reaction N at temperature T(K):
```
line_index = (T - 1)  # Since T starts at 1.0
value = ZONE[N].data[line_index]
```

**File Size:** ~43 MB (CORIA: 44 reactions × 15000 points × ~65 bytes/line)

---

### 4. `chemistry-Troe.dat` – Fall-off Correction Tables

**Purpose:** Pre-computed Troe fall-off correction factors for pressure-dependent reactions.

**Format:** Tecplot ASCII (same as `chemistry-Arrhenius.dat`).

```
TITLE = "Fall-off data"
VARIABLES = "Temperature", "k_inf", "k_0", "k_c", "F_cent"
ZONE T=Reaction7
I=15000, F=POINT
<T_1>  <k_inf_1>  <k_0_1>  <k_c_1>  <F_cent_1>
...
```

**Structure:**

| Component | Description |
|-----------|-------------|
| TITLE | File description (fixed: "Fall-off data") |
| VARIABLES | Column headers (5 values: T, k_inf, k_0, k_c, F_cent) |
| ZONE | One zone per **fall-off reaction** (not all reactions) |
| I=N | Number of temperature points (e.g., 15000) |
| Data rows | Pre-computed Troe fall-off parameters |

**Troe Parameters:**

| Column | Description | Formula |
|--------|-------------|---------|
| `Temperature` | Temperature (K) | 1.0 to 15000.0 |
| `k_inf` | High-pressure limit rate | Computed from Arrhenius A, b, Ea |
| `k_0` | Low-pressure limit rate | Third-order rate coefficient |
| `k_c` | Blending coefficient | Intermediate value for fall-off |
| `F_cent` | Centering factor | Troe fall-off correction factor |

**Troe Fall-off Theory:**

The Troe correction extends the 2-body rate to handle pressure effects:

$$k(P) = k_{\infty} \left( \frac{k_0 [M]}{1 + k_0 [M]} \right) F_c^x$$

where:
- $k_{\infty}$ = high-pressure limit (provided in `chemistry-Arrhenius.dat`)
- $k_0$ = low-pressure limit (from Troe parameters)
- $F_c$ = centering factor = $(1-a) \exp(-T/T^*) + a \exp(-T/T^{**}) + \exp(-T^{***}/T)$
- $[M]$ = third-body concentration
- $x$ = fall-off function exponent

This file provides the **pre-computed** Troe parameters to avoid runtime calculation.

**Parsing Notes:**
- Only **4 zones** in CORIA (reactions 7, 16, 19, 37 are fall-off types)
- Zones are labeled with **original reaction index** (e.g., `T=Reaction7`, `T=Reaction16`)
- Zone order in file **may not match** reaction order (must parse zone header)
- `F_cent` values ≈ 1.0 at low T, < 1.0 at moderate T (indicates correction strength)

**File Size:** ~8 MB (CORIA: 4 fall-off reactions × 15000 points × ~130 bytes/line)

---

### 5. `thermo.dat` – Thermodynamic Property Tables

**Purpose:** Pre-computed thermodynamic properties (Cp, H, S, dCp) for all species across temperature range.

**Format:** Tecplot ASCII (human-readable header, dense numeric data).

```
TITLE = "Mass Thermodynamic Properties"
VARIABLES = "Temperature", "Cp", "Enthalpy", "Entropy", "dCp"
ZONE T=H2
I=15000, F=POINT
<T_1>  <Cp_1>  <H_1>  <S_1>  <dCp_1>
<T_2>  <Cp_2>  <H_2>  <S_2>  <dCp_2>
...

ZONE T=H
I=15000, F=POINT
...

ZONE T=O
I=15000, F=POINT
...
```

**Structure:**

| Component | Description |
|-----------|-------------|
| TITLE | File description (fixed: "Mass Thermodynamic Properties") |
| VARIABLES | Column headers (5 values: T, Cp, Enthalpy, Entropy, dCp) |
| ZONE | One zone per species; `T=<species_name>` |
| I=N | Number of temperature points (e.g., 15000) |
| F=POINT | Data format (always POINT) |
| Data rows | Temperature and properties (space-separated floats) |

**Properties (per unit mass, not molar):**

| Column | Unit (SI) | Description |
|--------|-----------|-------------|
| `Temperature` | K | 1.0 to 15000.0 K (1 K increments) |
| `Cp` | J/(kg·K) | Heat capacity at constant pressure |
| `Enthalpy` | J/kg | Specific enthalpy (H - H_ref) |
| `Entropy` | J/(kg·K) | Specific entropy |
| `dCp` | J/(kg·K) | Temperature derivative of Cp (usually ~0 for ideal gas) |

**Example (H2 species, first 10 points):**

```
ZONE T=H2
I=15000, F=POINT
1.0 13654.381351 -4093931.342254 -13115.572143 0.000000
2.0 13654.381351 -4080276.960903 -3651.076207 0.000000
3.0 13654.381351 -4066622.579552 1885.299003 0.000000
4.0 13654.381351 -4052968.198202 5813.419728 0.000000
5.0 13654.381351 -4039313.816851 8860.306874 0.000000
6.0 13654.381351 -4025659.435500 11349.794939 0.000000
7.0 13654.381351 -4012005.054150 13454.627107 0.000000
8.0 13654.381351 -3998350.672799 15277.915664 0.000000
9.0 13654.381351 -3984696.291448 16886.170149 0.000000
10.0 13654.381351 -3971041.910098 18324.802809 0.000000
```

**Observations:**
- **Cp is constant** for each species (ideal gas approximation)
- **H varies linearly** with T (for ideal gas: $H = H_{ref} + C_p (T - T_{ref})$)
- **S varies logarithmically** with T (for ideal gas: $S = S_{ref} + C_p \ln(T/T_{ref})$)
- **dCp ≈ 0** (indicates ideal gas model, no temperature dependence of Cp)

**Zone Ordering:**

Species zones appear in **same order** as in `phase.txt`:

1. H2
2. H
3. O
4. O2
5. OH
6. ...
7. CH3O2 (last)

**File Size:** ~15 MB (CORIA: 17 species × 15000 points × ~60 bytes/line)

---

## Complete Example: CORIA Mechanism

### Files

**1. phase.txt:**

```
ideal-gas phase
H2      2.016000
H       1.008000
O      15.999000
O2     31.998000
OH     17.007000
H2O    18.015000
HO2    33.006000
CH3    15.035000
CH4    16.043000
CO     28.010000
CO2    44.009000
H2O2   34.014000
HCO    29.018000
CH2O   30.026000
CH3O   31.034000
C2H6   30.070000
CH3O2  47.033000
```

**2. chemistry-info.txt (excerpt):**

```
CORIA-CNRS
N.ro species = 17
N.ro reactions = 44

Reaction type
1 Arrhenius
2 Arrhenius
...
7 falloff-Troe
...
44 Arrhenius

Reaction definition
1 H2 1.0 0.0 0.0
1 H 0.0 1.0 0.0
1 O 1.0 0.0 0.0
1 O2 0.0 0.0 0.0
... (repeated for each species)
```

**3. chemistry-Arrhenius.dat:**

44 Tecplot zones with 15000 temperature points each.

**4. chemistry-Troe.dat:**

4 Tecplot zones for fall-off reactions (Reactions 7, 16, 19, 37).

**5. thermo.dat:**

17 Tecplot zones for species with properties at 15000 temperature points each.

---

## Native Input Path Usage

To use native files in FLINT:

1. **Ensure Cantera is NOT compiled in** (configure build without Cantera)
2. **Place all files in the same directory** (e.g., `INPUT/`)
3. **Configure FLINT** to point to this directory
4. **Load mechanism** at runtime:

```fortran
! Pseudo-code
call load_phase_native("INPUT/")          ! Reads phase.txt
call load_mechanism_native("INPUT/")      ! Reads chemistry-info.txt, .dat files
call load_thermo_native("INPUT/")         ! Reads thermo.dat
```

---