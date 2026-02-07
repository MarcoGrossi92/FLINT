# FLINT Theoretical Guide

FLINT is a comprehensive framework for modeling chemically reacting flows with accurate thermodynamic, transport, and kinetic properties. This guide provides the theoretical foundation for FLINT's capabilities in multi-species gas dynamics and chemical kinetics.

---

## Document Structure

This theoretical guide is organized into the following sections:

1. **[Thermodynamic and Transport Properties](theory_thermo.md)**: Mixture rules, equation of state, property evaluation
2. **[Finite-Rate Kinetics](theory_kinetics.md)**: Arrhenius reactions, Lindemann falloff, Troe formulation
3. **[Chemical Equilibrium](theory_equilibrium.md)**: NASA CEA algorithm for UV problems

For implementation details and API documentation, please refer to the FLINT code.

---

## Thermodynamic and Transport Properties

FLINT computes mixture properties from individual species data using established mixing rules. Key capabilities include:

**Thermodynamic Properties**

- Mixture density, gas constant, heat capacities ($c_p$, $c_v$, $\gamma$)
- Enthalpy, internal energy, entropy (absolute and sensible)
- Speed of sound
- Temperature from energy (Newton-Raphson inversion)

**Transport Properties**

- Dynamic viscosity via Wilke's mixing rule
- Thermal conductivity via Wilke's mixing rule
- Species diffusion coefficients (optional)

**For detailed formulations, see:** [Thermodynamic and Transport Properties](theory_thermo.md)

---

## Chemical Kinetics

FLINT provides two approaches to modeling chemical composition:

**Finite-Rate Kinetics**

Compute mass source terms for each species, accounting for:

- Elementary reactions: Arrhenius kinetics with modified temperature dependence
- Three-body reactions: Collision partners with species-specific efficiencies
- Pressure-dependent reactions: Lindemann and Troe falloff

**For detailed formulations, see:** [Finite-Rate Kinetics](theory_kinetics.md)

**Chemical Equilibrium**

Computes equilibrium compositions by thermodynamic optimization:

- UV problems, constant internal energy and volume (or density) 
- NASA CEA methodology

**For detailed formulations, see:** [Chemical Equilibrium](theory_equilibrium.md)

---
