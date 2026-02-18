# Verification Outcomes

This section presents verification results for FLINT obtained by comparing native implementations with reference solutions computed using Cantera.

The objective is to assess numerical consistency across thermodynamics, chemical kinetics, reactor integration, and equilibrium calculations for a wide range of chemical mechanisms.

## Batch Reactor Integration

Constant-volume batch reactor simulations were performed for multiple chemical mechanisms. For each case, temperature evolution obtained with FLINT was compared against Cantera reference solutions.

Three modeling approaches were evaluated:

- Cantera interface
- General chemistry subroutine
- Dedicated explicit chemistry routines

**Results**

Verification results for the different chemical mechanisms are presented in the figures below. As shown, the numerical results exhibit perfect agreement across all cases.

<div class="grid">

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

</div>

## Chemical Equilibrium

Constant-volume equilibrium simulations were performed for multiple set of species defined by the chemical mechanisms. For each case, the resulting temperature obtained with FLINT was compared against Cantera reference solutions.

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

</div>

---

## Reproducibility

All cases presented in this section are based on configurations included in the FLINT repository. Test programs and input data can be executed directly after installation.

For details on the testing infrastructure and framework, see the [Developer Guide](../development/testing.md) section.