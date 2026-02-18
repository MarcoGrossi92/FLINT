# Chemical Mechanism Database

FLINT provides a curated collection of chemical reaction mechanisms for combustion and high-temperature chemistry applications. These mechanisms range from simple global models (suitable for fast qualitative simulations) to detailed mechanisms (required for accurate prediction of minor species and pollutants).

Note that these mechanisms already have dedicated, explicit Fortran routines implemented in FLINT. Nevertheless, FLINT also provides a general subroutine that allows any mechanism to be executed, provided that FLINT uses the Cantera interface or that the appropriate input is supplied.

---

## Methane (CH₄) Mechanisms

Methane is the primary component of natural gas and a key fuel for power generation, rocket propulsion, and industrial processes.

---

### Smooke

A reduced kinetic mechanism for methane-air flames, widely used for premixed combustion modeling.

**Key Features:**

- **Developer**: Yale University (M. D. Smooke)
- **Number of species**: 16
- **Number of reactions**: 35
- **Primary fuel**: Methane (CH₄)
- **Applications**: Laminar premixed flames, flame structure studies
- **Temperature range**: 300–2500 K
- **Pressure range**: Atmospheric conditions

**Reference:**  
Smooke, M. D., ed. *Reduced Kinetic Mechanisms and Asymptotic Approximations for Methane-Air Flames: A Topical Volume*. Springer-Verlag, 1991.

---

### CORIA-CNRS

A RAMEC-based reduced mechanism designed for high-pressure methane combustion over a wide range of equivalence ratios.

**Key Features:**

- **Developer**: CORIA - CNRS, Normandie Université (France)
- **Number of species**: 17
- **Number of reactions**: 44
- **Primary fuel**: Methane (CH₄)
- **Equivalence ratio**: 0.2–14 (ultra-lean to ultra-rich)
- **Pressure range**: 1–100 bar (high-pressure validated)
- **Applications**: Rocket engines, gas turbines, supercritical combustion

**Reference:**  
Monnier, F., and Ribert, G. "Simulation of High-Pressure Methane-Oxygen Combustion with a New Reduced Chemical Mechanism." *Combustion and Flame*, vol. 235, 2022, 111735.

---

### Zhukov-Kong

A skeletal mechanism specifically developed for high-pressure methane-oxygen combustion in rocket engines.

**Key Features:**

- **Developer**: Institute of Space Propulsion, German Aerospace Centre (DLR)
- **Number of species**: 23
- **Number of reactions**: 51
- **Primary fuel**: Methane (CH₄)
- **Oxidizer**: Pure oxygen (O₂) or oxygen-enriched mixtures
- **Pressure range**: 10–300 bar (rocket engine conditions)
- **Applications**: Liquid rocket engines, high-pressure combustors

**Reference:**  
Zhukov, V. P., and Kong, A. F. "A Compact Reaction Mechanism of Methane Oxidation at High Pressures." *Progress in Reaction Kinetics and Mechanism*, vol. 43, no. 1, 2018, pp. 62-78.

---

### Westbrook-Dryer

A simplified kinetic model providing a computationally efficient representation of hydrocarbon combustion with minimal species.

**Key Features:**

- **Developer**: Lawrence Livermore National Laboratory, Princeton University
- **Number of species**: 5 (Fuel, O₂, CO₂, H₂O, N₂)
- **Number of reactions**: 2–4 (depending on variant)
- **Primary fuels**: Methane (CH₄), Ethane (C₂H₆), Propane (C₃H₈)
- **Pressure range**: 1–30 atm
- **Applications**: Large-scale CFD, turbulent combustion modeling (RANS, LES)

**Reference:**  
Westbrook, C. K., and Dryer, F. L. "Chemical Kinetic Modeling of Hydrocarbon Combustion." *Progress in Energy and Combustion Science*, vol. 10, no. 1, 1984, pp. 1-57.

---