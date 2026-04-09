# Theoretical Guide

FLINT is a comprehensive framework for modeling chemically reacting flows with accurate thermodynamic, transport, and kinetic properties. This guide provides the theoretical foundation for FLINT's capabilities in multi-species gas dynamics and chemical kinetics.

<div class="grid cards" markdown>

-   :material-thermometer:{ .lg .middle } __Thermodynamic and Transport Properties__

    ---

    Mixture rules, equation of state, property evaluation

    [:octicons-arrow-right-24: Explore thermodynamic models](thermo.md)

-   :material-chart-bell-curve:{ .lg .middle } __Finite-Rate Kinetics__

    ---

    Arrhenius reactions, Lindemann falloff, Troe formulation

    [:octicons-arrow-right-24: Explore kinetics models](kinetics.md)

-   :material-scale-balance:{ .lg .middle } __Chemical Equilibrium__

    ---

    NASA CEA algorithm for UV problems

    [:octicons-arrow-right-24: Explore equilibrium solver](equilibrium.md)

</div>


## Thermodynamic and Transport Properties

FLINT computes mixture properties from individual species data using established mixing rules. Key capabilities include:

**Thermodynamic Properties**

- Mixture density, gas constant, heat capacities ($c_p$, $c_v$, $\gamma$)
- Enthalpy, internal energy, entropy (absolute and sensible)
- Speed of sound
- Partial derivatives with respect to pressure and enthalpy for real fluids

**Transport Properties**

- Dynamic viscosity via Wilke's mixing rule
- Thermal conductivity via Wilke's mixing rule
- Species diffusion coefficients (optional)

**For detailed formulations, see:** [Thermodynamic and Transport Properties](thermo.md)

---

## Chemical Kinetics

FLINT provides two approaches to modeling chemical composition:

**Finite-Rate Kinetics**

Compute mass source terms for each species, accounting for:

- Elementary reactions: Arrhenius kinetics with modified temperature dependence
- Three-body reactions: Collision partners with species-specific efficiencies
- Pressure-dependent reactions: Lindemann and Troe falloff

**For detailed formulations, see:** [Finite-Rate Kinetics](kinetics.md)

**Chemical Equilibrium**

Computes equilibrium compositions by thermodynamic optimization:

- UV problems, constant internal energy and volume (or density) 
- NASA CEA methodology

**For detailed formulations, see:** [Chemical Equilibrium](equilibrium.md)

---
