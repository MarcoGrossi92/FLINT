# FLINT Algorithms and Computational Methods

## 1. Chemical Equilibrium Analysis (CEA)

### Overview
Calculates equilibrium composition at constant Temperature-Pressure (TP) or Temperature-Volume (TV).

### Mathematical Formulation
Minimize Gibbs free energy subject to element conservation:

$$\min_n G = \sum_i n_i(\hat{g}^°_i(T) + RT \ln(x_i P/P°))$$

Subject to:
- Element conservation: $\sum_i a_{ij} n_i = b_j$ for each element j
- Non-negativity: $n_i ≥ 0$ for all species i

### Solution Method
1. **Lagrangian**: Introduce multipliers λⱼ for element constraints
2. **Optimization**: Newton-Raphson iteration on reduced problem
3. **Convergence**: Check element balance to tolerance (e.g., 1e-10)
4. **Stability**: Handle zero/negative mole fractions, singular Jacobian

### Algorithm Steps
```
function cea_tp(T, P, X_initial):
  n ← X_initial * sum(n)              // Initial guess
  for iter = 1 to max_iter:
    compute G(n), ∂G/∂n
    compute constraints b_j(n)
    build Jacobian J = ∂b/∂n
    solve: J Δn = -b(n)               // Newton step
    if |b(n)| < tol: break
    n ← n + α Δn                       // Line search
  return n
end
```

### Numerical Considerations
- **Singular Jacobian**: Use SVD or regularization near constraints
- **Negative Mole Fractions**: Clamp to zero and remove species
- **Pressure Dependence**: Only appears in entropy term (small effect)

### References
- NASA CEA code (Gordon & McBride)
- Equilibrium Gibbs energy minimization (classical thermodynamics)

---

## 2. Finite-Rate Kinetics

### Overview
Solves ODEs for species concentrations including chemical reactions.

### Mathematical Formulation
$$\frac{dc_i}{dt} = \sum_{r} \left(\nu'_{i,r} - \nu_{i,r}\right) R_r$$

Where:
- $c_i$ = concentration of species i
- $\nu_{i,r}$, $\nu'_{i,r}$ = stoichiometric coefficients (reactant, product)
- $R_r$ = reaction rate for reaction r

### Reaction Rate (Arrhenius)
$$R_r = A_r T^{b_r} \exp\left(-\frac{E_{a,r}}{RT}\right) \prod_j c_j^{m_j}$$

Where:
- $A_r$ = pre-exponential factor
- $b_r$ = temperature exponent
- $E_{a,r}$ = activation energy
- $m_j$ = concentration order in reaction

### Pressure-Dependent Reactions (Troe)
For three-body and falloff reactions:

$$k(P) = k_∞ \frac{P_r}{1 + P_r} F^{(\log_{10} P_r, N)}$$

Where:
- $k_∞$ = high-pressure limit
- $P_r$ = reduced pressure
- $F$ = Troe correction factor

### ODE Integration
**Solver**: DASSL (implicit Runge-Kutta, variable order)

```
function kinetics_integrate(c_init, T, dt, tol):
  y ← [c_init]
  for time_step:
    call DASSL(y, f, df/dc, dt, tol)
    store output
  return concentrations vs time
end
```

**Parameters**:
- Relative tolerance: 1e-6
- Absolute tolerance: 1e-20 mol/cm³
- Max substeps: 1000

### Stiffness
- **Stiff**: Large Ea → fast equilibration
- **Mitigation**: Implicit solver (DASSL) + adaptive stepping

---

## 3. Transport Properties

### Viscosity (Chapman-Enskog)
$$\mu_i = \frac{5k_B T \sqrt{\pi m_i k_B T}}{16\pi\sigma_i^2 \Omega^{(2,2)}}$$

For mixtures:
$$\mu_{mix} = \sum_i \frac{\phi_i \mu_i}{\sum_j \phi_j \psi_{ij}}$$

Where $\phi_i$ = mole fraction, $\psi_{ij}$ = Wilke correction factor.

### Thermal Conductivity
$$\lambda = \lambda_{\text{heavy}} + \lambda_{\text{reactive}} + \lambda_{\text{Eucken}}$$

**Heavy particle contribution** (Chapman-Enskog):
$$\lambda_h = \frac{25k_B^2 T \sqrt{\pi m_i k_B T}}{32\pi\sigma_i^2 \Omega^{(2,2)}}$$

**Eucken correction** (for polyatomic molecules):
$$\lambda_E = \frac{f C_v}{M} \mu$$

### Binary Diffusion (Chapman-Enskog)
$$D_{ij} = \frac{3k_B T \sqrt{\pi m_i m_j / (m_i + m_j)}}{8\pi\sigma_{ij}^2 \Omega^{(1,1)}}$$

For multicomponent mixtures: Use mixture-averaged or full matrix solution.

### Temperature/Pressure Dependence
All transport properties scale as:
$$\mu, \lambda, D \propto T^{0.5} \text{ (approximately)}$$
Pressure independence (in Chapman-Enskog limit).

---

## 4. Thermodynamic Properties

### Enthalpy (NASA7 Polynomial)
$$\frac{H^°}{RT} = a_0 + \frac{a_1}{2}T + \frac{a_2}{3}T^2 + \frac{a_3}{4}T^3 + \frac{a_4}{5}T^4 + \frac{a_5}{T}$$

Heat capacity:
$$\frac{C_p^°}{R} = a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4$$

### Entropy (NASA7 Polynomial)
$$\frac{S^°}{R} = a_0 \ln T + a_1 T + \frac{a_2}{2}T^2 + \frac{a_3}{3}T^3 + \frac{a_4}{4}T^4 + a_6$$

### Gibbs Free Energy
$$G^° = H^° - TS^°$$

Used in CEA module for equilibrium calculations.

### Mixture Properties
$$\bar{H} = \sum_i x_i H_i^°, \quad \bar{S} = \sum_i x_i S_i^° - R \sum_i x_i \ln x_i$$

---

## Key Parameters & Tolerances

| Parameter | Default | Description |
|-----------|---------|-------------|
| CEA max iterations | 100 | Gibbs minimization |
| CEA tolerance | 1e-10 | Element conservation |
| Kinetics relative tol | 1e-6 | ODE solution |
| Kinetics absolute tol | 1e-20 | Species concentration floor |
| Transport T range | 200–5000 K | Valid for NASA polynomials |

---

## Validation Against Theory

- **CEA**: Compare to NIST CEA code output
- **Transport**: Compare to Cantera or NIST viscosity database
- **Kinetics**: Validate against shock tube or flame structure

---

## References

[Add your references here]
- NASA Glenn Report SP-273, Glenn Research Center (1994)
- Chapman-Enskog kinetic theory references
- Troe falloff formulation papers
```

### How to Customize
1. Replace `[Add your references here]` with real citations
2. Add figures (flowcharts for solver loops)
3. Include numerical examples from test cases
4. Add error analysis section if needed

---