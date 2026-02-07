# Chemistry - Equilibrium

## Overview

FLINT calculates chemical equilibrium compositions for multi-component gas mixtures using entropy maximization. The equilibrium solver supports the **UV problem**: constant internal energy and volume (or equivalently, constant specific internal energy and density).

This approach is based on the **NASA CEA (Chemical Equilibrium with Applications)** methodology developed by Gordon and McBride [1,2], which is the standard for aerospace and combustion equilibrium calculations.

---

## Mathematical Formulation

For a closed system at constant internal energy $U$ and volume $V$, chemical equilibrium corresponds to the state of **maximum entropy** $S$ [3]. The optimization problem is:

$$
\max_{\mathbf{n}, T} \, S(\mathbf{n}, T) = -\frac{1}{T} \sum_{i=1}^{N_s} n_i \left[\bar{h}_i^{\circ}(T) - T\bar{s}_i^{\circ}(T) - R_u T \ln\left(\frac{x_i P}{P^{\circ}}\right)\right]
$$

subject to:

**Element conservation constraints:**
$$
\sum_{i=1}^{N_s} a_{ij} \, n_i = b_j \quad \text{for } j = 1, 2, \ldots, N_e
$$

**Internal energy constraint:**
$$
\sum_{i=1}^{N_s} n_i \bar{u}_i(T) = U
$$

**Equation of state (for pressure):**
$$
P V = \left(\sum_{i=1}^{N_s} n_i\right) R_u T
$$

**Non-negativity constraints:**
$$
n_i \geq 0 \quad \text{for all } i = 1, 2, \ldots, N_s
$$

where:

- $\mathbf{n} = [n_1, n_2, \ldots, n_{N_s}]^T$ = vector of species mole numbers
- $T$ = temperature (unknown, to be determined)
- $\bar{h}_i^{\circ}(T)$ = standard-state molar enthalpy of species $i$ (J/mol)
- $\bar{s}_i^{\circ}(T)$ = standard-state molar entropy of species $i$ (J/(mol·K))
- $\bar{u}_i(T) = \bar{h}_i^{\circ}(T) - R_u T$ = molar internal energy of species $i$
- $x_i = n_i / \sum_k n_k$ = mole fraction of species $i$
- $P^{\circ}$ = standard pressure (typically 1 bar or 1 atm)
- $a_{ij}$ = number of atoms of element $j$ in one molecule of species $i$
- $b_j$ = total number of moles of element $j$ (conserved quantity)
- $U$ = total internal energy (specified)
- $V$ = total volume (specified)
- $N_s$ = number of species, $N_e$ = number of elements
- $R_u$ = universal gas constant (8.314 J/(mol·K))

Equivalently, this can be formulated as a **Gibbs free energy minimization** at the unknown equilibrium temperature and pressure [1].

---

## Solution Method

The constrained optimization problem is solved using the **method of Lagrange multipliers** combined with **Newton-Raphson iteration** on an extended system that includes temperature as an unknown [1,4].

### Lagrangian Formulation

Introduce Lagrange multipliers $\boldsymbol{\lambda} = [\lambda_1, \lambda_2, \ldots, \lambda_{N_e}]^T$ for element conservation and $\lambda_U$ for the energy constraint:

$$
\mathcal{L}(\mathbf{n}, T, \boldsymbol{\lambda}, \lambda_U) = G(\mathbf{n}, T, P) - \sum_{j=1}^{N_e} \lambda_j \left(\sum_{i=1}^{N_s} a_{ij} n_i - b_j\right) - \lambda_U \left(\sum_{i=1}^{N_s} n_i \bar{u}_i(T) - U\right)
$$

where $G$ is the Gibbs free energy and $P$ is determined from the equation of state.

At equilibrium, the stationarity conditions are:

$$
\frac{\partial \mathcal{L}}{\partial n_i} = \mu_i(T, P, \mathbf{x}) - \sum_{j=1}^{N_e} \lambda_j a_{ij} - \lambda_U \bar{u}_i(T) = 0 \quad \text{for all } i
$$

$$
\frac{\partial \mathcal{L}}{\partial T} = \sum_{i=1}^{N_s} n_i \frac{\partial \mu_i}{\partial T} - \lambda_U \sum_{i=1}^{N_s} n_i \bar{c}_{v,i}(T) = 0
$$

$$
\sum_{i=1}^{N_s} a_{ij} n_i - b_j = 0 \quad \text{for all } j
$$

$$
\sum_{i=1}^{N_s} n_i \bar{u}_i(T) - U = 0
$$

where $\mu_i$ is the chemical potential of species $i$ and $\bar{c}_{v,i}(T) = d\bar{u}_i/dT$ is the molar heat capacity at constant volume.

### Newton-Raphson Iteration

The CEA algorithm uses a **reduced Newton-Raphson method** that iterates on temperature and a subset of independent composition variables [1,2]:

1. **Variable transformation**: Use $\ln n_i$ and $\ln T$ to enforce positivity
2. **Reduced system**: Eliminate dependent variables using element constraints
3. **Jacobian construction**: Compute derivatives with respect to $\{\ln n_i, \ln T\}$
4. **Linear solve**: Solve for Newton correction: $\mathbf{J} \, \Delta \mathbf{z} = -\mathbf{r}$
5. **Line search**: Apply damping factor $\alpha \in (0, 1]$ to ensure descent and feasibility
6. **Convergence check**: Stop when residuals satisfy $\|\mathbf{r}\| < \epsilon$ (typically $\epsilon = 10^{-10}$)

---

## References

[1] Gordon, S., and McBride, B. J. "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications. Part 1: Analysis." NASA Reference Publication 1311, 1994.

[2] McBride, B. J., and Gordon, S. "Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications. Part 2: Users Manual and Program Description." NASA Reference Publication 1311, 1996.