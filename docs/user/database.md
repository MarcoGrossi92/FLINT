# Chemical Mechanism Database

FLINT includes a curated collection of chemical reaction mechanisms for combustion and high-temperature chemistry. These mechanisms are implemented as dedicated Fortran routines and range from simple global models (for fast simulations) to reduced/skeletal mechanisms (for moderate accuracy). The database covers a variety of fuels and applications, including hydrogen, methane, hydrocarbons, hybrid rocket fuels, and solid rocket motor exhaust chemistry.

---

## Quick Reference

| Mechanism | Type | Species | Reactions | Primary Application |
|-----------|------|---------|-----------|---------------------|
| **Frolov** | Global | 3 | 1 | Hydrogen combustion (ultra-fast) |
| **Nassini** | Global | 3 | 1 | Hydrogen combustion (ultra-fast) |
| **WD** | Global | 5 | 3 | CH₄ global reaction, CFD |
| **JLR** | Global | 9 | 7 | CH₄ rocket engines |
| **ONERA-7** | Reduced | 7 | 14 | H₂/air scramjet |
| **Smooke** | Reduced | 16 | 35 | CH₄ premixed flames |
| **CORIA-CNRS** | Reduced | 18 | 44 | CH₄ high-pressure combustion |
| **ZK** | Skeletal | 25 | 51 | CH₄ rocket engines (high-pressure) |
| **TSR-CDF-13** | Skeletal | 13 | 46 | CH₄ diffusion flames |
| **TSR-PSR-11** | Skeletal | 11 | 22 | CH₄ well-stirred reactor |
| **TSR-GP-24** | Skeletal | 24 | 110 | CH₄ general purpose |
| **TSR-Rich-31** | Skeletal | 31 | 197 | CH₄ rich combustion |
| **FFCMy_12** | Detailed | 38 | 291 | CH₄ detailed chemistry |
| **San Diego** | Detailed | – | – | Hydrocarbon detailed mechanism |
| **Coronetti** | Global | 9 | 6 | C₄H₆ HTPB hybrid rockets |
| **Singh** | Global | 10 | 11 | C₃₂H₆₆ paraffin wax hybrid rockets |
| **Cross** | Reduced | 20 | 33 | SRM plume (HCl/HCN) |
| **Ecker** | Reduced | 14 | 28 | SRM plume (HCl) |
| **Troyes** | Reduced | 12 | 17 | SRM plume (HCl/HCN) |
| **Pelucchi** | Detailed | 25 | 103 | Chlorine combustion (HCl/Cl₂) |

---

## Hydrogen Mechanisms

### ONERA-7

A reduced H₂/O₂ mechanism with 7 species developed for supersonic combustion applications.

**Characteristics:**  
- **Species / Reactions**: 7 / 14  
- **Temperature**: High-temperature hydrogen combustion   
- **Application**: Scramjet, supersonic combustor, hypersonic flow  
- **Accuracy**: Reduced mechanism, suitable for hypersonic simulations  
- **File**: `ONERA-7.f90`

---

## Methane Mechanisms

### Global Mechanisms

#### Westbrook-Dryer

A simplified global reaction model for methane combustion with minimal species.

**Characteristics:**
- **Species / Reactions**: 5 / 3
- **Temperature range**: 1–30 atm
- **Application**: Large-scale CFD, RANS/LES turbulent combustion modeling
- **Accuracy**: Global reaction, ultra-fast computation
- **File**: `WD.f90`

**Reference:**  
Westbrook, C.K., and Dryer, F.L. "Chemical Kinetic Modeling of Hydrocarbon Combustion." *Progress in Energy and Combustion Science*, 10(1), 1–57, 1984.

---

#### JLR (Rocket Engine Global)

A global mechanism for methane-oxygen combustion in rocket engines.

**Characteristics:**
- **Species / Reactions**: 9 / 7
- **Application**: Rocket engines, liquid propellant combustion, rapid estimation
- **Accuracy**: Global mechanism, fast computation
- **File**: `JLR.f90`

---

### Reduced Mechanisms (Moderate Complexity)

#### Smooke

A reduced kinetic mechanism widely used for premixed flame structure studies.

**Characteristics:**
- **Species / Reactions**: 16 / 35
- **Temperature range**: 300–2500 K
- **Pressure**: Atmospheric conditions
- **Application**: Laminar premixed flames, flame structure, pollutant precursors
- **Accuracy**: Reduced, good balance between speed and detail
- **File**: `smooke.f90`

**Reference:**  
Smooke, M.D. (ed.) *Reduced Kinetic Mechanisms and Asymptotic Approximations for Methane-Air Flames: A Topical Volume*. Springer-Verlag, 1991.

---

#### CORIA-CNRS

A RAMEC-based reduced mechanism for high-pressure methane combustion.

**Characteristics:**
- **Species / Reactions**: 18 / 44
- **Equivalence ratio**: 0.2–14 (ultra-lean to ultra-rich)
- **Pressure range**: 1–100 bar
- **Application**: High-pressure rocket engines, gas turbines, supercritical combustion
- **Accuracy**: Reduced, validated for high pressure
- **File**: `coria.f90`

**Reference:**  
Monnier, F., and Ribert, G. "Simulation of High-Pressure Methane-Oxygen Combustion with a New Reduced Chemical Mechanism." *Combustion and Flame*, 235, 111735, 2022.

---

### Skeletal Mechanisms (Temperature-Sensitive-Reduction)

#### ZK (Zhukov-Kong)

A skeletal mechanism for high-pressure methane-oxygen combustion in rocket engines.

**Characteristics:**
- **Species / Reactions**: 25 / 51
- **Pressure range**: 10–300 bar
- **Oxidizer**: Pure O₂ or oxygen-enriched mixtures
- **Application**: Liquid rocket engines, high-pressure combustors
- **Accuracy**: Skeletal, validated for rocket conditions
- **File**: `ZK.f90`

**Reference:**  
Zhukov, V.P., and Kong, A.F. "A Compact Reaction Mechanism of Methane Oxidation at High Pressures." *Progress in Reaction Kinetics and Mechanism*, 43(1), 62–78, 2018.

---

#### TSR-CDF-13 (Diffusion Flame)

A skeletal mechanism developed via Temperature-Sensitive-Reduction for counterflow diffusion flame simulations.

**Characteristics:**
- **Species / Reactions**: 13 / 46
- **Application**: Counterflow diffusion flames, flame sheet models
- **Accuracy**: Skeletal, optimize for diffusion flame chemistry
- **File**: `TSR-CDF-13.f90`

---

#### TSR-GP-24 (General Purpose)

A general-purpose skeletal mechanism for methane combustion across multiple flame types.

**Characteristics:**
- **Species / Reactions**: 24 / 110
- **Pressure range**: 1–40 bar
- **Application**: General combustion simulation, good balance of accuracy and speed
- **Accuracy**: Skeletal, versatile
- **File**: `TSR-GP-24.f90`

---

#### TSR-Rich-31 (Rich Combustion)

A skeletal mechanism optimized for rich methane combustion regimes.

**Characteristics:**
- **Species / Reactions**: 31 / 197
- **Application**: Rich mixtures, fuel-lean combustion, NO formation
- **Accuracy**: Skeletal, detailed intermediate chemistry
- **File**: `TSR-Rich-31.f90`

---

### Detailed Mechanisms

#### FFCMy_12 (Flamelet-Generated Manifold)

A detailed methane mechanism for comprehensive combustion modeling.

**Characteristics:**
- **Species / Reactions**: 38 / 291
- **Fuel**: Methane (CH₄)
- **Application**: Detailed prediction, engineering simulations with moderate complexity
- **Accuracy**: Detailed, computationally more intensive
- **File**: `FFCMy_12.f90`

---

## Hydrocarbon Mechanisms

### sandiego20161214

A comprehensive hydrocarbon mechanism for multi-fuel combustion.

**Characteristics:**
- **Fuels**: Hydrocarbons (C₁–C₄ and beyond)
- **Application**: Detailed hydrocarbon chemistry, multiple fuel types
- **Accuracy**: Detailed mechanism
- **File**: `sandiego20161214.f90`

---

## Hybrid Rocket Fuel Mechanisms

### Coronetti (HTPB - Global)

A global mechanism for HTPB (hydroxyl-terminated polybutadiene) combustion.

**Characteristics:**
- **Species / Reactions**: 9 / 6
- **Fuel**: C₄H₆ (HTPB energy release surrogate)
- **Application**: Quick HTPB simulations, hybrid rocket motors, preliminary design
- **Accuracy**: Global mechanism, ultra-fast
- **File**: `coronetti.f90`

---

### Singh (Paraffin Wax - Global)

A global mechanism for paraffin wax combustion in hybrid rocket motors.

**Characteristics:**
- **Species / Reactions**: 10 / 11
- **Fuel**: C₃₂H₆₆ (paraffin wax surrogate)
- **Application**: Paraffin-based hybrid rockets, rapid analysis
- **Accuracy**: Global mechanism
- **File**: `singh.f90`

---

## Solid Rocket Motor (SRM) Plume Mechanisms

Mechanisms for HCl and halocarbon chemistry in aluminized solid rocket motor exhaust.

### Troyes

A mechanism for SRM plume chemistry with HCl and HCN formation.

**Characteristics:**
- **Species / Reactions**: 12 / 17
- **Chemical focus**: HCl, HCN formation pathways
- **Application**: SRM plume simulation, nozzle exhaust, plume chemistry
- **File**: `troyes.f90`

---

### Ecker

A simplified mechanism for HCl chemistry in SRM plumes.

**Characteristics:**
- **Species / Reactions**: 14 / 28
- **Chemical focus**: HCl formation, simplified sub-mechanism
- **Application**: Fast SRM plume estimates, HCl formation analysis
- **File**: `ecker.f90`

---

### Cross

A mechanism combining HCl and HCN formation in SRM exhaust.

**Characteristics:**
- **Species / Reactions**: 20 / 33
- **Chemical focus**: HCl, HCN formation, combined pathways
- **Application**: Detailed SRM plume chemistry, exhaust analysis
- **File**: `cross.f90`

### Pelucchi (Chlorine / HCl Oxidation)

A detailed mechanism for HCl and Cl₂ chemistry at high temperatures.

**Characteristics:**
- **Species / Reactions**: 25 / 103
- **Chemical focus**: Chlorine combustion, HCl oxidation, detailed sub-mechanisms
- **Temperature range**: High-temperature conditions
- **Application**: Chlorine chemistry studies, safety analysis, specialized combustion
- **Accuracy**: Detailed mechanism
- **File**: `pelucchi.f90`

---

## Selection Guide

**Choose based on your application:**

| Need | Recommended | Why |
|------|-------------|-----|
| **Ultra-fast 3D CFD** | WD, JLR, global-H2 | Minimal species (1–9) |
| **Premixed flame structure** | Smooke, CORIA-CNRS | Proven for laminar flames, reduced detail |
| **Rocket engine (CH₄)** | ZK, CORIA-CNRS, JLR | Validated at high pressure |
| **General combustion** | TSR-GP-24, Smooke | Balanced accuracy and speed |
| **Diffusion flames** | TSR-CDF-13 | Optimized for diffusion-dominated regimes |
| **Rich combustion** | TSR-Rich-31 | Better intermediates for fuel-rich conditions |
| **HTPB hybrid rockets** | Coronetti, Singh | Fast hybrid rocket simulation |
| **SRM exhaust** | Troyes, Ecker, Cross | HCl/HCN chemistry tailored to SRM conditions |
| **Chlorine chemistry** | Pelucchi | Specialized high-temperature halogen chemistry |
| **Detailed hydrocarbons** | FFCMy_12, sandiego | More comprehensive species set |

---