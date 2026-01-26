# YAML Input Format Specification

Complete reference for FLINT YAML format (Cantera-enabled builds).

## Table of Contents

1. [Overview](#overview)
2. [Top-Level Sections](#top-level-sections)
   - [Units](#units)
   - [Phases](#phases)
   - [Elements](#elements)
   - [Species](#species)
   - [Reactions](#reactions)
3. [Complete Example](#complete-example)
4. [Tips & Best Practices](#tips--best-practices)
5. [Troubleshooting](#troubleshooting)

---

## Overview

A FLINT input file (`.yaml`) defines a complete chemical system:
- Chemical elements present
- Species composition and properties
- Thermodynamic data (enthalpy, entropy)
- Transport coefficients (viscosity, diffusion)
- Reaction mechanisms (if applicable)
- Operating conditions (temperature, pressure, phase state)

### File Structure

```yaml
units:         # Unit system definitions
  ...
phases:        # Gas/liquid/solid phase definitions
  ...
elements:      # Chemical elements (H, C, O, etc.)
  - H
  - C
  - O
species:       # Species definitions with thermo and transport
  - name: ...
    thermo: ...
    transport: ...
reactions:     # Chemical reactions (optional)
  - reactants: ...
    products: ...
```

---

## Top-Level Sections

### 1. `units` (Required)

Defines the unit system used throughout the file.

**Parameters:**

| Key | Values | Description | Example |
|-----|--------|-------------|---------|
| `length` | `m`, `cm`, `mm`, `in`, `ft` | Length/diameter units | `cm` |
| `quantity` | `mol`, `kmol` | Molar quantity units | `mol` |
| `energy` | `J`, `cal`, `eV` | Energy units for thermodynamics | `cal` |
| `activation-energy` | `J/mol`, `cal/mol`, `K` | Activation energy units | `cal/mol` |
| `temperature` | `K` | Temperature (always Kelvin) | `K` |

**Example:**

```yaml
units:
  length: cm
  quantity: mol
  activation-energy: cal/mol
  energy: cal
  temperature: K
```

**Important Notes:**
- All energies must be in the same unit system
- Temperature is always in Kelvin (cannot be changed)
- Diameter and well-depth (transport) use the `length` and `energy` units specified

---

### 2. `phases` (Required)

Defines the thermodynamic phases present in the system.

**Phase Definition:**

```yaml
phases:
  - name: <phase_name>
    thermo: <model>
    species: <selection>
    reactions: <selection>
    kinetics: <solver>
    transport: <model>
    state:
      T: <temperature_K>
      P: <pressure_Pa>
      X: {<species>: <mole_fraction>, ...}  # Optional initial composition
```

**Parameters:**

| Key | Allowed Values | Description |
|-----|---|---|
| `name` | String | Unique phase identifier |
| `thermo` | `ideal-gas`, `stoichiometric-solid`, `ideal-condensed` | Thermodynamic model |
| `species` | `all` or list | Species included in phase (usually `all`) |
| `reactions` | `all` or list | Reactions occurring in phase |
| `kinetics` | `gas-phase`, `surface`, `edge` | Kinetics type |
| `transport` | `mixture-averaged`, `multi-component` | Transport model |

**Example:**

```yaml
phases:
  - name: gas
    thermo: ideal-gas
    species: all
    reactions: all
    kinetics: gas-phase
    transport: mixture-averaged
    state:
      T: 300.0
      P: 101325.0  # 1 atm in Pa
```

---

### 3. `elements` (Required)

Declares the chemical elements used in the mechanism.

**Format:**

```yaml
elements:
  - H
  - C
  - O
  - N
  - Ar
```

**Supported Elements:**

All periodic table elements are supported. Common ones:
- **Common:** H, C, N, O, S, Cl, Ar, He
- **Metals:** Fe, Cu, Ni, Al, Na
- **Halogens:** F, Cl, Br, I

**Note:** Every atom used in any species must be declared here.

---

### 4. `species` (Required)

Defines individual chemical species with composition, thermodynamic properties, and transport data.

**Species Definition:**

```yaml
species:
  - name: <species_name>
    composition: {<element>: <count>, ...}
    thermo:
      model: <thermo_model>
      temperature-ranges: [<T1>, <T2>, <T3>, ...]
      data: [...]  # Model-specific data
    transport:
      model: <transport_model>
      <transport_parameters>
```

#### Composition

Specifies the atomic makeup of the species.

```yaml
composition: {H: 2, O: 1}  # H₂O
composition: {C: 1, O: 2}   # CO₂
composition: {H: 1}         # H (atom)
```

#### Thermodynamic Models

##### **NASA7** (Most Common)

7-coefficient polynomial model for enthalpy and entropy.

```yaml
thermo:
  model: NASA7
  temperature-ranges: [200.0, 1000.0, 3500.0]  # Low, mid, high
  data:
    # Low temperature (200–1000 K)
    - [a0, a1, a2, a3, a4, a5, a6]
    # High temperature (1000–3500 K)
    - [a0, a1, a2, a3, a4, a5, a6]
```

**Data Format (NASA7 coefficients):**
- Coefficients: `[a₀, a₁, a₂, a₃, a₄, a₅, a₆]`
- Used in polynomials:
  - $C_p^° = R(a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4)$
  - $H^° = RT(a_0 + \frac{a_1}{2} T + \frac{a_2}{3} T^2 + \frac{a_3}{4} T^3 + \frac{a_4}{5} T^4 + \frac{a_5}{T})$
  - $S^° = R(a_0 \ln T + a_1 T + \frac{a_2}{2} T^2 + \frac{a_3}{3} T^3 + \frac{a_4}{4} T^4 + a_6)$

**Example:**

```yaml
- name: H2O
  composition: {H: 2, O: 1}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 5000.0]
    data:
      - [4.19864056E+00, -2.03643537E-03, 6.52040211E-06, -5.48797062E-09, 1.77197367E-12, -3.02937267E+04, -8.49032283E-01]
      - [3.03399649E+00, 2.17691804E-03, -1.64072518E-07, -9.70419870E-11, 1.68200992E-14, -3.00042971E+04, 4.96677010E+00]
```

##### **NASA9** (Alternative)

9-coefficient polynomial (higher accuracy over wider range).

```yaml
thermo:
  model: NASA9
  temperature-ranges: [200.0, 1000.0, 6000.0]
  data:
    - [a0, a1, ..., a8]  # Low T
    - [a0, a1, ..., a8]  # High T
```

#### Transport Models

##### **Gas Phase Transport**

Chapman-Enskog kinetic theory parameters.

```yaml
transport:
  model: gas
  geometry: <atom|linear|nonlinear>
  diameter: <sigma_Angstrom>           # Lennard-Jones diameter (Å)
  well-depth: <epsilon_Kelvin>         # Well-depth (K)
  polarizability: <alpha_Angstrom3>    # Polarizability (Ų, optional)
  rotational-relaxation: <Zrot>        # Rotational relaxation number (optional)
  dipole: <mu_Debye>                   # Dipole moment (Debye, optional)
```

**Parameters:**

| Parameter | Units | Description | Typical Range |
|-----------|-------|-------------|---|
| `geometry` | — | Molecular shape: `atom`, `linear`, `nonlinear` | — |
| `diameter` | Å | Lennard-Jones collision diameter | 2.0–4.0 |
| `well-depth` | K | Lennard-Jones well-depth | 50–500 |
| `polarizability` | Ų | Molecular polarizability | 0.5–2.0 |
| `rotational-relaxation` | — | Zrot (dimensionless) | 0.5–2.0 |
| `dipole` | Debye | Dipole moment for polar molecules | 0–3 |

**Example:**

```yaml
- name: O2
  composition: {O: 2}
  transport:
    model: gas
    geometry: linear
    diameter: 3.46
    well-depth: 107.4
    polarizability: 1.60
    rotational-relaxation: 3.8
```

---

### 5. `reactions` (Optional)

Defines chemical reactions and rate parameters.

**Reaction Definition:**

```yaml
reactions:
  - reactants: <reactant_string>
    products: <product_string>
    type: <elementary|three-body|falloff|...>
    rate-coefficients:
      A: <Arrhenius_prefactor>
      b: <temperature_exponent>
      Ea: <activation_energy>
    duplicate: <true|false>        # For duplicate reactions
    negative-A: <true|false>       # For reverse reaction (advanced)
```

#### Reaction String Format

Reactions use chemical notation with optional stoichiometric coefficients:

```
"2 H2 + O2 => 2 H2O"
"H + O2 => OH + O"
"M + H2 => H + H + M"  # Third-body (M = any species)
```

**Rules:**
- Use `=>` for forward reaction (not `→` or `<=>`)
- Stoichiometric coefficients default to 1 if omitted
- Third-body reactions use `M` as placeholder

#### Arrhenius Rate Expression

Forward reaction rate coefficient:

$$k_f = A T^b \exp\left(-\frac{E_a}{RT}\right)$$

**Parameters:**

| Key | Units (depends on `units.activation-energy`) | Description |
|-----|---|---|
| `A` | cm³/(mol·s) or equivalent | Pre-exponential factor |
| `b` | — | Temperature exponent |
| `Ea` | As specified in `units` | Activation energy |

**Example:**

```yaml
reactions:
  - reactants: "H2 + O2"
    products: "2 OH"
    rate-coefficients:
      A: 1.7e13
      b: 0.0
      Ea: 48000    # in cal/mol (per units.activation-energy)

  - reactants: "H + O2 + M"
    products: "HO2 + M"
    rate-coefficients:
      A: 2.1e15
      b: -0.6
      Ea: 0.0      # No activation energy barrier
    third-body-efficiencies:
      H2: 2.0
      H2O: 6.0
      AR: 0.7
```

#### Reaction Types

##### **Elementary Reaction**

Single-step reaction with Arrhenius parameters.

```yaml
- reactants: "H + O2"
  products: "OH + O"
  rate-coefficients:
    A: 5.06e04
    b: 2.67
    Ea: 6292.0
```

##### **Three-Body Reaction**

Requires a third-body collision. Use `M` in reactants/products.

```yaml
- reactants: "H + H2 + M"
  products: "H2 + H + M"
  rate-coefficients:
    A: 1.0e18
    b: -1.0
    Ea: 0.0
  third-body-efficiencies:
    H2: 2.5
    H2O: 12.0
    AR: 0.0
```

##### **Falloff Reaction** (Troe)

Pressure-dependent rate with low- and high-pressure limits.

```yaml
- reactants: "H + O2 + M"
  products: "HO2 + M"
  type: falloff
  low-pressure-coefficients:
    A: 6.366e20
    b: -1.72
    Ea: 524.8
  high-pressure-coefficients:
    A: 1.475e12
    b: 0.6
    Ea: 0.0
  Troe:
    a: 0.5
    T3: 1e-30
    T1: 1.0e+100
```

---

## Complete Example

### CORIA Mechanism (Simplified)

```yaml
units:
  length: cm
  quantity: mol
  activation-energy: cal/mol
  energy: cal
  temperature: K

phases:
  - name: gas
    thermo: ideal-gas
    species: all
    reactions: all
    kinetics: gas-phase
    transport: mixture-averaged
    state:
      T: 300.0
      P: 101325.0

elements:
  - H
  - C
  - O

species:
  - name: H2
    composition: {H: 2}
    thermo:
      model: NASA7
      temperature-ranges: [200.0, 1000.0, 3500.0]
      data:
        - [2.34433112E+00, 7.98052075E-03, -1.94781510E-05, 2.01572094E-08, -7.37611761E-12, -9.17935173E+02, 6.83010238E-01]
        - [3.33727920E+00, -4.94024731E-05, 4.99456778E-07, -1.79566394E-10, 2.00255376E-14, -9.50158922E+02, -3.20502331E+00]
    transport:
      model: gas
      geometry: linear
      diameter: 2.920
      well-depth: 38.0
      polarizability: 0.79
      rotational-relaxation: 280.0

  - name: O2
    composition: {O: 2}
    thermo:
      model: NASA7
      temperature-ranges: [200.0, 1000.0, 3500.0]
      data:
        - [3.78245636E+00, -2.99673416E-03, 9.84730201E-06, -9.68129509E-09, 3.24372837E-12, -1.06394356E+03, 3.65767573E+00]
        - [3.28253784E+00, 1.48308754E-03, -7.57966669E-07, 2.09470555E-10, -2.16717794E-14, -1.08845772E+03, 5.45323129E+00]
    transport:
      model: gas
      geometry: linear
      diameter: 3.458
      well-depth: 107.4
      polarizability: 1.60
      rotational-relaxation: 3.8

  - name: OH
    composition: {O: 1, H: 1}
    thermo:
      model: NASA7
      temperature-ranges: [200.0, 1000.0, 3500.0]
      data:
        - [3.99201543E+00, -2.40131752E-03, 4.61793841E-06, -3.88113333E-09, 1.26418070E-12, 3.61508056E+03, -1.03925458E+00]
        - [3.09288767E+00, 5.48429716E-04, 1.26505228E-07, -8.79461556E-11, 2.38336365E-15, 3.65767573E+03, 4.47669610E+00]
    transport:
      model: gas
      geometry: linear
      diameter: 2.750
      well-depth: 80.0

reactions:
  - reactants: "H2 + O2"
    products: "2 OH"
    rate-coefficients:
      A: 1.7e13
      b: 0.0
      Ea: 48000.0

  - reactants: "H + O2"
    products: "OH + O"
    rate-coefficients:
      A: 5.06e04
      b: 2.67
      Ea: 6292.0

  - reactants: "H + H2O"
    products: "H2 + OH"
    rate-coefficients:
      A: 1.00e08
      b: 1.60
      Ea: 13630.0
```

---

## Tips & Best Practices

### 1. Unit Consistency
- Define all energy units upfront in the `units` section
- Don't mix cal/mol and J/mol
- Use consistent diameter units (Ångstroms is standard)

### 2. Thermodynamic Data Source
- NASA polynomials: Use data from NASA CEA code or Cantera database
- Verify temperature ranges don't exceed 200–5000 K for NASA7
- Reference papers for custom data

### 3. Transport Data Source
- Lennard-Jones parameters: NIST, Cantera, or literature
- Check if polarizability is needed (usually for polar molecules like H₂O)
- Typical dipole moments: 0 (nonpolar) to 3 D (polar)

### 4. Reaction Mechanisms
- Elementary steps only (no multi-step approximations)
- Use reliable mechanisms (e.g., GRI-Mech 3.0, KAUST Surrogate, Pelucchi)
- Document reaction source in comments

### 5. Validation
- Test at known equilibrium points (CEA module)
- Compare computed properties (Cp, μ, λ) to literature or Cantera
- Verify flame speeds match experimental data

---

## Troubleshooting

### YAML Format Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `Species not found` | Species used in reaction not defined | Add species to `species` section |
| `Element mismatch` | Composition uses undefined elements | Add element to `elements` section |
| `Temperature out of range` | Evaluation outside NASA7 range | Extend `temperature-ranges` or provide additional data |
| `Negative concentration` | Numerical instability | Reduce time step or improve initial guess |
| `Invalid YAML syntax` | Malformed YAML structure | Check indentation, quotes, and key names |

---

## See Also

- [Cantera Documentation](https://cantera.org) – Official Cantera reference
