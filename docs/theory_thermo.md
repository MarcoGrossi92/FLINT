# Thermodynamic and Transport Properties

## Overview

FLINT supports the definition of thermodynamic and transport properties for a mixture of $N_s$ thermally perfect gases following the thermal equation of state:

$$
p = \rho R_\text{mix} T
$$

where $p$ is pressure, $\rho$ is density, $R_\text{mix}$ is the mixture gas constant, and $T$ is temperature.

**Property evaluation methods:**

- **Cantera integration**: Theoretical background for properties computed via Cantera can be found in the official Cantera documentation [1].

- **Native tabulated data**: The native structure uses temperature-varying tabulated properties for each species: specific heat capacity at constant pressure $c_{p,\text{tab}}$, specific enthalpy $h_{\text{tab}}$, specific entropy $s_{\text{tab}}$, dynamic viscosity $\mu_{\text{tab}}$, and thermal conductivity $k_{\text{tab}}$.

In addition to thermodynamic and transport properties, FLINT provides conversions between primitive and conservative variables for computational fluid dynamics (CFD) applications.

---

## Table Interpolation

For tabulated properties, linear interpolation is applied between integer temperatures $T_i$ and $T_{i+1}$ surrounding the requested temperature $T$ to retrieve any property $f$ for each species:

$$ 
f(T) = f(T_i) + \frac{f(T_{i+1}) - f(T_i)}{T_{i+1} - T_i} \cdot (T - T_i) 
$$

For integer-spaced tables where $T_{i+1} - T_i = 1$, this simplifies to:

$$ 
f(T) = f(T_i) + (f(T_{i+1}) - f(T_i)) \cdot (T - T_i) 
$$

---

## Thermodynamic Properties

Once individual species properties are known, mixture thermodynamic quantities are computed using mass-weighted averaging [2,3].

**Mixture density:**
$$ 
\rho = \sum_{s=1}^{N_s}\rho_s 
$$

**Mixture gas constant:**
$$
R_\text{mix} = \sum_{s=1}^{N_s} Y_s R_s = \sum_{s=1}^{N_s} \frac{\rho_s}{\rho} R_s
$$

where $Y_s = \rho_s/\rho$ is the mass fraction of species $s$, and $R_s = R_u/M_s$ is the specific gas constant ($R_u = 8314.46$ J/(kmol·K) is the universal gas constant and $M_s$ is the molecular weight).

**Mixture specific heat capacity at constant pressure:**
$$
c_{p,\text{mix}} = \sum_{s=1}^{N_s} Y_s c_{p,s}
$$

**Specific heat capacity at constant volume:**
$$
c_{v,\text{mix}} = c_{p,\text{mix}} - R_\text{mix}
$$

**Heat capacity ratio (specific heat ratio):**
$$
\gamma = \frac{c_{p,\text{mix}}}{c_{v,\text{mix}}} = \frac{c_{p,\text{mix}}}{c_{p,\text{mix}} - R_\text{mix}}
$$

**Speed of sound:**
$$ 
a = \sqrt{\gamma R_\text{mix} T} 
$$

**Mixture specific enthalpy:**
$$ 
h(\rho_s, T) = \sum_{s=1}^{N_s} Y_s h_s(T)
$$

**Total specific enthalpy (including kinetic energy):**
$$
h_0(\rho_s, T, \mathbf{u}) = h(\rho_s, T) + \frac{1}{2} |\mathbf{u}|^2
$$

where $\mathbf{u}$ is the velocity vector.

**Mixture specific internal energy:**
$$ 
e(\rho_s, T) = h(\rho_s, T) - \frac{p}{\rho} = \sum_{s=1}^{N_s} Y_s h_s(T) - R_\text{mix} T
$$

**Total specific internal energy:**
$$
e_0(\rho_s, T, \mathbf{u}) = e(\rho_s, T) + \frac{1}{2} |\mathbf{u}|^2
$$

---

## Transport Properties

Once individual species transport properties are known, mixture transport properties are computed using **Wilke's mixing rule** [4,5], which accounts for molecular interactions between different species.

**Wilke's interaction parameter:**
$$ 
\phi_{ij} = \frac{\left[1 + \left(\mu_i/\mu_j\right)^{1/2} \left(M_j/M_i\right)^{1/4}\right]^2}{\sqrt{8 \left(1 + M_i/M_j\right)}} 
$$

where $\mu_i$ is the dynamic viscosity of species $i$, $M_i$ is the molecular weight of species $i$, and $X_i$ is the mole fraction of species $i$.

**Mixture dynamic viscosity:**
$$ 
\mu_\text{mix} = \sum_{i=1}^{N_s} \frac{X_i \mu_i}{\sum_{j=1}^{N_s} X_j \phi_{ij}} 
$$

**Mixture thermal conductivity:**
$$ 
k_\text{mix} = \sum_{i=1}^{N_s} \frac{X_i k_i}{\sum_{j=1}^{N_s} X_j \phi_{ij}} 
$$

where $k_i$ is the thermal conductivity of species $i$.

**Note on mole fractions:** The mole fraction $X_s$ is related to mass fraction $Y_s$ by:
$$
X_s = \frac{Y_s / M_s}{\sum_{j=1}^{N_s} Y_j / M_j}
$$

---

## Primitive ↔ Conservative Variable Conversion

For CFD applications, conversions between primitive variables (density, velocity, pressure) and conservative variables (conserved mass, momentum, energy) are essential [6].

### Primitive → Conservative

$$
\begin{aligned}
(\rho_s)^\text{cons} &= (\rho_s)^\text{prim} \\
(\rho \mathbf{u})^\text{cons} &= \rho^\text{prim} \mathbf{u}^\text{prim} \\
(\rho e_0)^\text{cons} &= \rho^\text{prim} \left[ e(\rho_s^\text{prim}, T^\text{prim}) + \frac{1}{2} |\mathbf{u}^\text{prim}|^2 \right]
\end{aligned}
$$

### Conservative → Primitive

The conversion from conservative to primitive variables requires solving for temperature $T$ iteratively:

1. Extract species densities: $\rho_s = (\rho_s)^\text{cons}$
2. Compute total density: $\rho = \sum_{s=1}^{N_s} \rho_s$
3. Compute velocity: $\mathbf{u} = (\rho \mathbf{u})^\text{cons} / \rho$
4. Compute specific internal energy: $e = (\rho e_0)^\text{cons}/\rho - \frac{1}{2}|\mathbf{u}|^2$
5. **Solve for temperature** using Newton–Raphson iteration:
   $$
   e(\rho_s, T) = \sum_{s=1}^{N_s} Y_s h_s(T) - R_\text{mix}(Y_s) T
   $$
   Find $T$ such that $e(\rho_s, T) = e_\text{known}$
6. Compute pressure: $p = \rho R_\text{mix} T$

---

## References

[1] Cantera Development Team. *Cantera: An Object-oriented Software Toolkit for Chemical Kinetics, Thermodynamics, and Transport Processes*. https://cantera.org, Version 3.0+, 2024.

[2] Poinsot, T., and Veynante, D. *Theoretical and Numerical Combustion*, 3rd edition. Published by the authors, 2012.

[3] Kee, R. J., Coltrin, M. E., and Glarborg, P. *Chemically Reacting Flow: Theory and Practice*, 2nd edition. John Wiley & Sons, 2003.

[4] Wilke, C. R. "A Viscosity Equation for Gas Mixtures." *The Journal of Chemical Physics*, vol. 18, no. 4, 1950, pp. 517-519. https://doi.org/10.1063/1.1747673

[5] Bird, R. B., Stewart, W. E., and Lightfoot, E. N. *Transport Phenomena*, 2nd edition. John Wiley & Sons, 2002.

[6] Blazek, J. *Computational Fluid Dynamics: Principles and Applications*, 3rd edition. Butterworth-Heinemann, 2015.