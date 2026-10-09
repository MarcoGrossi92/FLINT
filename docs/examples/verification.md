# Verification Outcomes

This section presents verification results for FLINT obtained by comparing native implementations with reference solutions computed using Cantera.

The objective is to assess numerical consistency across thermodynamics, chemical kinetics, reactor integration, and equilibrium calculations for a wide range of chemical mechanisms.

## Batch Reactor Integration

Constant-volume batch reactor simulations were performed for multiple chemical mechanisms. For each case, temperature evolution obtained with FLINT was compared against Cantera reference solutions.

Four solutions are compared:

- Cantera: Cantera's own reactor (C++ `IdealGasReactor`, rtol 1e-10, atol 1e-15), the reference
- FLINT Cantera: FLINT's integrator with the source terms of the Cantera interface
- FLINT Explicit: FLINT's integrator with the dedicated explicit chemistry routine
- FLINT General: FLINT's integrator with the general chemistry subroutine

**Results**

Verification results for the different chemical mechanisms are presented in the figures below. As shown, the numerical results exhibit perfect agreement across all cases.

<div class="grid">

<figure>
  {% include "examples/images/WD.svg" %}
  <figcaption>Westbrook-Dryer (no general procedure: its orders are not in chemistry-info.txt)</figcaption>
</figure>

<figure>
  {% include "examples/images/Troyes.svg" %}
  <figcaption>Troyes</figcaption>
</figure>

<figure>
  {% include "examples/images/Ecker.svg" %}
  <figcaption>Ecker</figcaption>
</figure>

<figure>
  {% include "examples/images/Cross.svg" %}
  <figcaption>Cross</figcaption>
</figure>

<figure>
  {% include "examples/images/Pelucchi.svg" %}
  <figcaption>Pelucchi</figcaption>
</figure>

<figure>
  {% include "examples/images/TSR-CDF-13.svg" %}
  <figcaption>TSR-CDF-13</figcaption>
</figure>

<figure>
  {% include "examples/images/TSR-GP-24.svg" %}
  <figcaption>TSR-GP-24</figcaption>
</figure>

<figure>
  {% include "examples/images/TSR-Rich-31.svg" %}
  <figcaption>TSR-Rich-31</figcaption>
</figure>

<figure>
  {% include "examples/images/Smooke.svg" %}
  <figcaption>Smooke</figcaption>
</figure>

<figure>
  {% include "examples/images/CORIA.svg" %}
  <figcaption>CORIA-CNRS</figcaption>
</figure>

<figure>
  {% include "examples/images/ZK.svg" %}
  <figcaption>Zhukov-Kong</figcaption>
</figure>

<figure>
  {% include "examples/images/Gerlinger.svg" %}
  <figcaption>Gerlinger</figcaption>
</figure>

</div>

## Chemical Equilibrium

Constant-volume equilibrium (constant internal energy and volume) was computed for the species of several mechanisms: over a sweep of the O2/CH4 mixture ratio at 1000 K and 3.25 kg/m³ (Westbrook-Dryer, Zhukov-Kong, TSR-GP-24), and over a sweep of pressure from 1e-5 to 100 bar at 3000 K (Ecker). The equilibrium temperature of FLINT is compared with Cantera's (`equilibrate("UV")`).

**Results**

Verification results are presented in the figures below. As shown, the numerical outcomes exhibit perfect agreement across all cases.

<div class="grid">

<figure>
  {% include "examples/images/WD-eq.svg" %}
  <figcaption>Westbrook-Dryer</figcaption>
</figure>

<figure>
  {% include "examples/images/ZK-eq.svg" %}
  <figcaption>Zhukov-Kong</figcaption>
</figure>

<figure>
  {% include "examples/images/TSR-GP-24-eq.svg" %}
  <figcaption>TSR-GP-24</figcaption>
</figure>

<figure>
  {% include "examples/images/Ecker-eq.svg" %}
  <figcaption>Ecker (pressure sweep, 3000 K)</figcaption>
</figure>

</div>

---

## Reproducibility

All cases presented in this section are based on configurations included in the FLINT repository. Test programs and input data can be executed directly after installation.

For details on the testing infrastructure and framework, see the [Developer Guide](../development/testing.md) section.