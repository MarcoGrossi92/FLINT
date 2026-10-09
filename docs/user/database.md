# Chemical Mechanism Database

FLINT includes a curated collection of chemical reaction mechanisms for combustion and high-temperature chemistry. These mechanisms are implemented as dedicated Fortran routines and range from simple global models (for fast simulations) to reduced/skeletal mechanisms (for moderate accuracy). The database covers a variety of fuels and applications, including hydrogen, methane, hydrocarbons, hybrid rocket fuels, and solid rocket motor exhaust chemistry.

---

## Quick Reference

The **mechanism name** is the name that `Assign_Mechanism` hooks: line 1 of `chemistry-info.txt` must
hold it exactly for the compiled routine to be used (any other name runs the general procedure). The
species and reaction counts are those of the routine, as recorded in
`src/lib/Lib_ChemMech/mechanism_contract.json`; the reactions are counted per table type
(Arrhenius-type, including three-body, / falloff-Troe / falloff-Lindemann).

| Mechanism name | Routine (file) | Type | Species | Reactions (Arrh. / Troe / Lind.) | Primary application | Reference |
|----------------|----------------|------|---------|----------------------------------|---------------------|-----------|
| `Frolov` | `Frolov` (`global-H2.f90`) | Global | 3 (+ inert) | 1 (rate hard-coded, tables ignored) | Hydrogen combustion (ultra-fast) | Frolov, Dubrovskii, Ivanov, *Progress in Propulsion Physics* 4 (2013) 467-488, doi:10.1051/eucass/201304467, Eq. (11) |
| `Frolov_nopressure` | `Frolov_nopressure` (`global-H2.f90`) | Global | 4 | 1 (reversible, tables; CFD++ variant) | Hydrogen combustion (ultra-fast), analytical Jacobian | variant of Frolov et al. 2013 (no pressure factor, A doubled) |
| `Nassini` | `Nassini_4` (`global-H2.f90`) | Global | 3 (+ inert) | 2 irreversible (tables at 1 atm) | Hydrogen combustion (ultra-fast) | Nassini, "High-fidelity Numerical Investigations of a Hydrogen Rotating Detonation Combustor", PhD thesis, University of Florence (XXXIV cycle, 2018-2021); Nassini, Andreini, Bohon, *Combust. Flame* 258 (2023) 113050, doi:10.1016/j.combustflame.2023.113050 |
| `ONERA-7` | `ONERA_7` (`ONERA-7.f90`) | Reduced | 7 | 14 / 0 / 0 | H₂/air scramjet, analytical Jacobian | Davidenko, Gökalp, Dufour, Magre, AIAA 2006-7913, doi:10.2514/6.2006-7913, Table A.1 |
| `Gerlinger-9` | `Gerlinger9` (`Gerlinger-9.f90`) | Reduced | 9 | 19 / 0 / 0 | H₂/air supersonic combustion | Gerlinger, Möbus, Brüggemann, *J. Comput. Phys.* 167 (2001) 247-276, doi:10.1006/jcph.2000.6671 (modified Jachimowski 1988) |
| `WD` | `WD` (`WD.f90`) | Global | 5 | 3 / 0 / 0 | CH₄ global reaction, CFD | Westbrook, Dryer, *Prog. Energy Combust. Sci.* 10 (1984) 1-57 |
| `WD-Andersen` | `Andersen` (`WD.f90`) | Global | 5 | 3 / 0 / 0 | CH₄ global reaction with the Andersen CO/CO₂ steps | Andersen, Rasmussen, Giselsson, Glarborg, *Energy Fuels* 23(3) (2009) 1379-1389, doi:10.1021/ef8003619 |
| `OSK` | `OSK` (`WD.f90`) | Global, one step | 4 | 1 / 0 / 0 | CH₄ one-step ([CH₄]^0.7 [O₂]^0.8) | to be confirmed |
| `JLR-Nasuti` | `JLR` (`JLR.f90`) | Global | 9 | 7 / 0 / 0 | CH₄ rocket engines | Jones, Lindstedt, *Combust. Flame* 73 (1988) 233-249; Betti et al., *AIAA J.* 54(5) (2016) 1693-1703 |
| `Frassoldati` | `Frassoldati` (`JLR.f90`) | Global | 9 | 6 / 0 / 0 | CH₄ rocket engines | Jones-Lindstedt scheme with other rate parameters; to be confirmed |
| `CKJLR-10sp` | `CKJLR10sp` (`JLR.f90`) | Global | 10 | 8 / 0 / 0 | C₁₂H₂₄ fuel | to be confirmed |
| `Smooke` | `smooke` (`smooke.f90`) | Reduced | 16 | 35 / 0 / 0 | CH₄ premixed flames | Smooke (ed.), Springer, 1991 |
| `CORIA-CNRS` | `coria` (`coria.f90`) | Reduced | 17 | 40 / 4 / 0 | CH₄ high-pressure combustion | Monnier, Ribert, *Combust. Flame* 235 (2022) 111735 |
| `ZK` | `ZK` (`ZK.f90`) | Skeletal | 25 | 45 / 6 / 0 | CH₄ rocket engines (high pressure) | Zhukov, Kong, *Prog. React. Kinet. Mech.* 43(1) (2018) 62-78, doi:10.3184/146867818X15066862094914 |
| `TSR-CDF-13` | `TSRCDF13` (`TSR-CDF-13.f90`) | Skeletal | 13 | 43 / 3 / 0 | CH₄ diffusion flames | Liberatori et al., *J. Propul. Power* 40(2) (2024) 303-319, doi:10.2514/1.B39283 |
| `TSR-GP-24` | `TSRGP24` (`TSR-GP-24.f90`) | Skeletal | 24 | 102 / 8 / 0 | CH₄ general purpose | Liberatori et al. 2024 (as above) |
| `TSR-Rich-31` | `TSRRich31` (`TSR-Rich-31.f90`) | Skeletal | 31 | 185 / 12 / 0 | CH₄ rich combustion | Liberatori et al. 2024 (as above) |
| `FFCMy-12` | `FFCMy_12` (`FFCMy_12.f90`) | Reduced | 13 (12 + N₂) | 34 / 3 / 1 | CH₄ combustion | Xu et al., *Combust. Flame* 263 (2024) 113380, doi:10.1016/j.combustflame.2024.113380 |
| `SanDiego` | `sandiego20161214` (`sandiego20161214.f90`) | Detailed | 57 | 245 / 23 / 0 | Hydrocarbon detailed mechanism | San Diego Mechanism, University of California at San Diego, version 2016-12-14 |
| `CoronettiC4H6` | `Coronetti` (`coronetti.f90`) | Global | 9 | 6 / 0 / 0 | C₄H₆ HTPB hybrid rockets | Coronetti, Sirignano, *J. Propul. Power* 29(2) (2013) 371-384, doi:10.2514/1.B34760 |
| `Singh` | `singh` (`singh.f90`) | Quasi-global | 9 | 10 / 0 / 0 | C₂H₄ combustion | Singh, Jachimowski, *AIAA J.* 32(1) (1994) 213-216, doi:10.2514/3.11972 |
| `Singh-WC32` | `Singh_WC32` (`singh.f90`) | Quasi-global | 10 | 11 / 0 / 0 | C₃₂H₆₆ paraffin wax hybrid rockets | Migliorino, Bianchi, Nasuti, *J. Propul. Power* 36(6) (2020) 806-819, doi:10.2514/1.B37914; Singh, Jachimowski 1994 |
| `Cross` | `cross` (`cross.f90`) | Reduced | 19 | 33 / 0 / 0 | SRM plume (HCl/HCN) | to be confirmed |
| `Ecker` | `ecker` (`ecker.f90`) | Reduced | 14 | 28 / 0 / 0 | SRM plume (HCl) | Ecker, Karl, Hannemann, EUCASS 2019 |
| `Troyes` | `troyes` (`troyes.f90`) | Reduced | 12 | 17 / 0 / 0 | SRM plume (HCl/HCN) | Troyes et al., AIAA 2006-4414, doi:10.2514/6.2006-4414 |
| `Pelucchi` | `pelucchi` (`pelucchi.f90`) | Detailed | 25 | 97 / 2 / 4 | Chlorine combustion (HCl/Cl₂) | Pelucchi et al., *Combust. Flame* 162(6) (2015) 2693-2704, doi:10.1016/j.combustflame.2015.04.002 |

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

### Global hydrogen schemes (`global-H2.f90`)

- **Frolov** (`Frolov`): one global step 2 H₂ + O₂ → 2 H₂O with the rate law **hard-coded** in the
  routine (A = 8e11 on the progress rate, (p/101325)^-1.15, Ea/R = 10000 K); the kf/kb tables of
  `chemistry-Arrhenius.dat` are **ignored**, the INPUT folder only supplies the species, their
  molecular weights and thermodynamics. Species slots `O2, H2O, H2` (+ inert species after them);
  this layout exists since commit 6e12db2 (the older routine addressed slots 2/3/5).
- **Frolov_nopressure** (`Frolov_nopressure`): the **CFD++ variant** of the global step (the reaction
  panel/file used in CFD++), not the paper formula: 2 H₂ + O₂ ⇌ 2 H₂O reversible, kf/kb read from the
  tables (forward A = 8e11 on the progress rate, i.e. twice the published Frolov rate at 1 atm, no
  pressure dependence), slots `O2, H2O, H2, N2`; carries an analytical Jacobian.
- **Nassini** (`Nassini_4`): **two irreversible reactions** read from the Arrhenius tables at 1 atm
  (temperature-only tables): reaction 1 H₂ + ½ O₂ → H₂O with rate kf₁[H₂][O₂], reaction 2
  H₂O → H₂ + ½ O₂ with rate kf₂[H₂O]. Since commit dc1caa3 (merged in 6e12db2) reaction 2 is the
  backward step (kb is not used); the older routine had a single reaction with its reverse from
  kb, which the tables of an irreversible reaction leave at zero. Slots `O2, H2O, H2` (+ inert).

The species slots and the reaction counts of every hooked routine are checked at
`Assign_Mechanism` (see *Mechanism contract check* in the chemistry routines page).

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

#### Westbrook-Dryer with the Andersen closure (`WD-Andersen`)

The Westbrook-Dryer steps 1-2 with the CO2 dissociation step written as the explicit inverse of the
CO oxidation step: rate = k3(T) [CO2] [H2O]^0.5 [O2]^-0.25 (the former FLINT law was [CO2]^1.25), so
that the pair 2/3 reaches the equilibrium of CO + 0.5 O2 <-> CO2 (k2/k3 is within 0.5 % of Kc between
1100 and 2000 K and 2 % at 3000 K on tables written from the mechanism's yaml). A species with a
negative order at zero
concentration gives a zero rate of that step: the convention of Cantera, one rule for every FLINT site
with a negative order (see the chemistry routines page). In Cantera's yaml format the step carries
`orders: {CO2: 1.0, H2O: 0.5, O2: -0.25}` with `negative-orders: true` and `nonreactant-orders: true`.

**Characteristics:**
- **Species / Reactions**: 5 / 3 (slots CH4, O2, CO2, H2O, CO)
- **File**: `WD.f90` (routine `Andersen`, mechanism name `WD-Andersen`)
- **Test**: `chemistry/andersen` (`test-andersen`: Cantera references on the tables of `test/chemistry/andersen/WD-Andersen`)

**Reference:**  
Andersen, J., Rasmussen, C.L., Giselsson, T., Glarborg, P. *Energy & Fuels*, 23(3), 1379–1389, 2009, DOI 10.1021/ef8003619.

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
- **Species / Reactions**: 17 / 44
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

### Reduced Mechanism from FFCM-2

#### FFCMy-12 (Foundational Fuel Chemistry Model, reduced)

A reduced methane mechanism derived from an early version (FFCMy) of FFCM-2, the Foundational Fuel
Chemistry Model.

**Characteristics:**
- **Species / Reactions**: 13 (12 + N₂) / 38 (34 Arrhenius-type including three-body, 3 falloff-Troe, 1 falloff-Lindemann)
- **Fuel**: Methane (CH₄)
- **Application**: Engineering simulations that need more than a global scheme at moderate cost
- **Accuracy**: Reduced
- **File**: `FFCMy_12.f90` (mechanism name `FFCMy-12`)

**Reference:**  
Xu, R., et al. *Combustion and Flame*, 263, 113380, 2024, DOI 10.1016/j.combustflame.2024.113380.

---

## Hydrocarbon Mechanisms

### sandiego20161214

A comprehensive hydrocarbon mechanism for multi-fuel combustion.

**Characteristics:**
- **Species / Reactions**: 57 / 268 (245 Arrhenius-type, 23 falloff-Troe), version of 2016-12-14
- **Fuels**: Hydrocarbons (C₁–C₄ and beyond)
- **Application**: Detailed hydrocarbon chemistry, multiple fuel types
- **Accuracy**: Detailed mechanism
- **File**: `sandiego20161214.f90` (mechanism name `SanDiego`)

**Reference:**  
San Diego Mechanism, University of California at San Diego, version 2016-12-14.

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

### Singh and Singh-WC32 (Ethylene and Paraffin Wax - Quasi-global)

Two routines of `singh.f90`:

- **Singh** (routine `singh`): the quasi-global ethylene scheme, 9 species / 10 reactions.
- **Singh-WC32** (routine `Singh_WC32`): the same scheme plus the cracking of the paraffin wax surrogate
  C₃₂H₆₆, 10 species / 11 reactions.

**Characteristics (Singh-WC32):**
- **Species / Reactions**: 10 / 11
- **Fuel**: C₃₂H₆₆ (paraffin wax surrogate)
- **Application**: Paraffin-based hybrid rockets, rapid analysis
- **Accuracy**: Quasi-global mechanism
- **File**: `singh.f90`

**References:**  
Singh, D.J., and Jachimowski, C.J. *AIAA Journal*, 32(1), 213–216, 1994, DOI 10.2514/3.11972.  
Migliorino, M.T., Bianchi, D., and Nasuti, F. *Journal of Propulsion and Power*, 36(6), 806–819, 2020, DOI 10.2514/1.B37914.

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
- **Species / Reactions**: 19 / 33
- **Chemical focus**: HCl, HCN formation, combined pathways
- **Application**: Detailed SRM plume chemistry, exhaust analysis
- **File**: `cross.f90`
- **Reference**: to be confirmed

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
| **HTPB / paraffin hybrid rockets** | CoronettiC4H6, Singh-WC32 | Fast hybrid rocket simulation |
| **SRM exhaust** | Troyes, Ecker, Cross | HCl/HCN chemistry tailored to SRM conditions |
| **Chlorine chemistry** | Pelucchi | Specialized high-temperature halogen chemistry |
| **Detailed hydrocarbons** | SanDiego (FFCMy-12 as a reduced CH₄ alternative) | More comprehensive species set |

---