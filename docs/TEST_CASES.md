# FLINT Test Cases

Complete validation test suite for FLINT. Each test case validates specific capabilities.

## Test Suite Overview

| Case | Physical System | Purpose | Files |
|------|-----------------|---------|-------|
| **CORIA** | CH₄/O₂ laminar flame | Finite-rate kinetics validation | test/CORIA/ |
| **Smooke** | Flame structure | Transport properties validation | test/Smooke/ |
| **TSR-CDF-13** | Shock tube reactor | CEA + kinetics stiffness | test/TSR-CDF-13/ |
| **WD** | Detonation wave | High-speed reaction validation | test/WD/ |
| **Cross** | Crossover reactions | Complex chemistry | test/Cross/ |

---

## CORIA Test Case

### Physical Description
**System**: 1D laminar premixed methane-air flame
**Fuel**: CH₄
**Oxidizer**: O₂ (21%) + N₂ (79%)
**Equivalence ratio**: φ = 1.0 (stoichiometric)
**Inlet temperature**: 300 K
**Outlet temperature**: ~2050 K (adiabatic flame temperature)

### Purpose
Validates:
- Finite-rate kinetics solver
- Transport property computation
- Multispecies diffusion
- Heat transfer calculations

### Reference
F. Monnier and G. Ribert, "A rapid compression machine study of autoignition of syngas/air mixtures," *Combustion and Flame*, 2021.

### Input File
Location: `test/CORIA/INPUT/coria.yaml`
- 21 species (H₂, H, O, OH, H₂O, HO₂, CH₄, CO, CO₂, etc.)
- 53 elementary reactions
- GRI-Mech compatible

### How to Run
```bash
cd FLINT
./install.sh build --compiler=gnu
cd test/CORIA
../../bin/test-FLINT INPUT/coria.yaml
```

**Expected runtime**: ~30 seconds (depending on system)

### Expected Output
Output files in `test/CORIA/OUTPUT/`:
- `temperature.dat` – Temperature profile [x (cm), T (K)]
- `species.dat` – Species concentrations [x, c₁, c₂, ...]
- `heatrelease.dat` – Volumetric heat release [x, Q̇]
- `flint.log` – Solver diagnostics

### Key Variables & Expected Ranges

| Variable | Description | Units | Expected Value | Tolerance |
|----------|-------------|-------|---|---|
| **T_inlet** | Inlet temperature | K | 300 | Fixed |
| **T_peak** | Peak flame temperature | K | 2030–2080 | ±50 |
| **T_out** | Exit temperature | K | 1950–2050 | ±100 |
| **x_H2O_max** | Peak water mole fraction | — | 0.11–0.13 | ±0.01 |
| **x_CO_max** | Peak CO mole fraction | — | 0.02–0.04 | ±0.01 |
| **x_OH_max** | Peak OH mole fraction | — | 5e-3–8e-3 | ±1e-3 |
| **Flame speed** | Laminar burning velocity | cm/s | 38–42 | ±2 |

### Plots
Use Tecplot with `test/CORIA/plot.lay`:
```
tecplot plot.lay
```

Shows:
- Temperature profile (smooth, single peak)
- Major species (H₂O increasing, fuel decreasing)
- Radicals (OH, H, O peak at flame front)
- Heat release (sharp peak at flame)

### Success Criteria

✅ **Acceptable results if:**
1. Temperature rises smoothly from 300 K to ~2050 K
2. Species profiles are monotonic (no oscillations)
3. Radical (OH) peaks at flame front, drops in products
4. Flame speed within ±5% of literature (38–42 cm/s)
5. No negative concentrations
6. Solver converges in <2000 iterations

❌ **Red flags:**
- Negative concentrations → Solver tolerance too loose
- Oscillating temperature → Mesh/time-stepping issue
- Flat temperature profile → Chemistry not converged
- Multiple peaks → Spurious numerical effects

### Troubleshooting

**Flame too cold (T_peak < 1900 K):**
- Check reaction mechanism (Ea values reasonable?)
- Verify thermodynamic data (H°, S° reasonable?)
- Reduce solver tolerance (tol_rel = 1e-8 instead of 1e-6)

**Flame too hot (T_peak > 2150 K):**
- Check heat loss model (radiative cooling implemented?)
- Verify boundary conditions (outlet pressure set correctly?)

**Divergence (| Error | > 100%):**
- Initial guess poor; try running from previous solution
- Time step too large; reduce dt

### Data Extraction
```fortran
! Typical post-processing
read(*, *) x, T, c_H2O, c_OH, c_CO, ...
flame_speed = integrate(Q_dot) / rho_inlet
peak_OH = max(c_OH)
```

---

## Smooke Test Case

### Physical Description
**System**: Stagnation-point diffusion flame (1D)
**Fuel stream**: CH₄
**Air stream**: O₂/N₂
**Configuration**: Counterflow
**Strain rate**: a = 100–500 s⁻¹

### Purpose
Validates:
- Multi-species diffusion (full matrix)
- Thermal diffusion effects
- Sensitivity to strain rate
- Transport property accuracy at high temperature

### Reference
[Add reference]

### Expected Results
[Add expected temperature, species profiles, extinction strain rate]

---

## TSR-CDF-13 Test Case

### Physical Description
**System**: Shock Tube Reactor (0D, uniform mixture)
**Fuel**: H₂/O₂
**Initial conditions**:
- T₀ = 1000 K
- P₀ = 1 atm
- φ = 1.0

### Purpose
Validates:
- CEA at high temperature
- Fast kinetics (stiff ODE system)
- Pressure rise during ignition
- Induction time predictions

### Expected Results
[Add timeline, temperature rise, species evolution]

---

## Additional Test Cases

### Cross [Minimal description]
### WD [Minimal description]
### ZK [Minimal description]
### Others [...]

---

## Validation Checklist

Use this to verify your installation:

```
Run CORIA:
  [ ] Temperature profile smooth?
  [ ] Peak T ≈ 2050 K?
  [ ] H₂O production ≈ 0.12?
  [ ] No negative species?

Run Smooke:
  [ ] Flame anchored at stagnation point?
  [ ] Temperature profile S-shaped?

Run TSR:
  [ ] Ignition delay predictions reasonable?
  [ ] Pressure rise physical?

Compare to literature:
  [ ] Flame speeds within 5%?
  [ ] Ignition delay correlation matches?
```

---

## Extending Test Suite

To add your own test case:

1. **Create directory**: `test/MY_CASE/`
2. **Add input file**: `test/MY_CASE/INPUT/mycase.yaml`
3. **Document**: Add section to TEST_CASES.md
4. **Run validation**: Compare against experimental or CFD data
5. **Add plot file**: `test/MY_CASE/plot.lay` for Tecplot visualization

---

## References

[Add your experimental/simulation references]
```

### Customize For Your Cases
1. Add actual reference papers
2. Include real expected values from your runs
3. Add comparison plots (experimental vs. simulation)
4. Include troubleshooting tips from your experience

---