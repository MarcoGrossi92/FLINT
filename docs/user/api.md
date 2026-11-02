# API Reference

This page lists the procedures and the module variables that a solver uses, module by module, with their arguments
and return codes. [Using the Library](using.md) shows how they fit together.

Conventions:

- Values are `real(8)` in SI units: K, Pa, kg/m³, J/kg, J/(kg·K), kg/kmol.
- Arrays over species run over `1:ns`, in the order of `phase.txt`. A mixture state is given by its partial
  densities `rhoi(ns)` (kg/m³), with `rho = sum(rhoi)`.
- The loaders read from the folder given, or from `INPUT/` when `folder` is absent, and prepend the module variable
  `FLINT_phase_prefix` (empty by default) to every file name. They return an `integer`: 0 on success.
- A table with a temperature axis has one row per kelvin: row `T` holds the value at `T` K, from `Tmin` to `Tmax`.
  Between rows the value is interpolated linearly. Out of the range, temperatures are clamped to it.

---

## Loading data

Module `FLINT_Load_ThermoTransport`, except `read_chemistry` (module `FLINT_Load_Chemistry`). The input files are
described in [Native Format](input/native.md). Each table is read from `<name>.dat` (Tecplot ASCII) or `<name>.szplt`
(Tecplot binary, which needs FLINT built with TecIO).

### Ideal gas

`err = read_idealgas_thermo([folder])` reads `phase.txt` (species names and molecular weights), `thermo.dat`
($c_p$, $h$, $s$ per species) and `composition.txt` (elemental composition, used by the equilibrium solver).

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | `phase.txt` not found |
| 2 | `phase.txt` cannot be read |
| 3 | Neither `thermo.dat` nor `thermo.szplt` found |
| 4 | The thermo table cannot be read, or a zone is not on the 1 K grid |
| 5 | `composition.txt` not found: the thermo tables are loaded, the equilibrium solver cannot be used |
| 6 | `composition.txt` cannot be read: as for 5 |

`err = read_idealgas_transport([folder])` reads `transport.dat`: viscosity and thermal conductivity per species. Call
it after `read_idealgas_thermo`.

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Neither `transport.dat` nor `transport.szplt` found |
| 2 | The table cannot be read |
| 3 | A zone is not on the 1 K grid, or the table does not cover the thermo range (it must end at the last thermo temperature and start at or below the first) |

`err = read_idealgas_diffusion([folder])` reads `diffusion.dat`: the binary diffusion coefficients $D_{ij}$, one zone
per species pair in the order (1,2), (1,3), …, (1,ns), (2,3), …, at the reference pressure given in its title
(`Pref=<value> Pa`; 101325 Pa without it). Call it after `read_idealgas_thermo`.

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Neither `diffusion.dat` nor `diffusion.szplt` found, or a single species (no pairs) |
| 2 | The table cannot be read |
| 3 | A zone is not on the 1 K grid, or the table does not cover the thermo range |
| 4 | The number of zones is not ns(ns−1)/2 |

`err = read_idealgas_composition([folder])` reads `composition.txt` alone (0 success, 1 not found, 2 cannot be
read); `read_idealgas_thermo` calls it.

### Real fluid

`err = read_realfluid_thermo([folder])` reads the species name from `phase.txt` and the $(p, h)$ table `thermo.dat`.

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | `phase.txt` not found |
| 2 | `phase.txt` cannot be read |
| 3 | Neither `thermo.dat` nor `thermo.szplt` found |
| 4 | The table cannot be read |
| 5 | The table has more than one block |

`err = read_realfluid_transport([folder])` reads the viscosity and thermal conductivity on the same $(p, h)$ grid.

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Neither `transport.dat` nor `transport.szplt` found |
| 2 | The table cannot be read |
| 3 | The table has more than one block |
| 4 | The grid is not that of `thermo.dat` |

`err = ph2pT()` builds the $(p, T)$ table of the enthalpy, `h_tab2D`, that `pT2h` reads, by inverting `T_tab`. Call it
once, after `read_realfluid_thermo`: 0 success, 1 already built, 2 grid too small, 3 real-fluid thermo not loaded, 4 a
row of `T_tab` cannot be inverted.

### Chemistry

`err = read_chemistry([folder], [mech_name])` reads `chemistry-info.txt` (the reactions, their third-body efficiencies
and, optionally, their reaction orders) and the rate tables `chemistry-Arrhenius.dat`, `chemistry-Troe.dat` and
`chemistry-Lindemann.dat`. `mech_name` (`character(*)`, out) returns the mechanism name, the first line of
`chemistry-info.txt`. Call it after `read_idealgas_thermo`.

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | `chemistry-info.txt` not found |
| 2 | `chemistry-info.txt` cannot be read, or gives reaction orders for a falloff reaction |
| 3 | A rate table cannot be read |
| 4 | A rate table has fewer zones than its reactions |
| 5 | A falloff table is not admissible: $k_\infty$ or $k_0$ negative, or $F_\text{cent}$ not finite |
| 6 | A zone is not on the 1 K grid, or the table does not cover the thermo range |

`call free_chemistry_data()` (module `FLINT_Lib_Chemistry_data`) deallocates the chemistry tables, before loading
another mechanism.

---

## Module variables

Module `FLINT_Lib_Thermodynamic`. The loaders set them; a solver reads them.

| Variable | Description | Unit |
|---|---|---|
| `ns` | Number of species | |
| `species_names(ns)` | Species names | |
| `wm_tab(ns)` | Molecular weights | kg/kmol |
| `Ri_tab(ns)` | Specific gas constants, `Runiv/wm_tab` | J/(kg·K) |
| `Runiv` | Universal gas constant, 8314.46261815324 (exact SI value) | J/(kmol·K) |
| `ne`, `elements_names(ne)`, `species_composition(ne, ns)` | Elements and elemental composition | |
| `Tmin`, `Tmax` | First and last row of the temperature tables | K |
| `cp_tab(Tmin:Tmax, ns)` | Specific heat at constant pressure | J/(kg·K) |
| `h_tab(Tmin:Tmax, ns)` | Specific enthalpy | J/kg |
| `s_tab(Tmin:Tmax, ns)` | Specific entropy | J/(kg·K) |
| `dcpi_tab(Tmin:Tmax, ns)` | The `dCp` column of `thermo.dat` | |
| `mi_tab(Tmin:Tmax, ns)` | Dynamic viscosity | Pa·s |
| `k_tab(Tmin:Tmax, ns)` | Thermal conductivity | W/(m·K) |
| `dij_tab(Tmin:Tmax, ndij)` | Binary diffusion coefficients at `dij_pref`, `ndij = ns(ns−1)/2` pairs | m²/s |
| `dij_pref` | Reference pressure of `dij_tab` | Pa |

Real fluid, on the $(p, h)$ grid `pmin + i*deltap`, `hmin + j*deltah`, with `i` and `j` from 0:

| Table | Property | Unit |
|---|---|---|
| `rho_tab` | Density | kg/m³ |
| `T_tab` | Temperature | K |
| `dT_tab` | $(\partial\rho/\partial T)_p$ | kg/(m³·K) |
| `rh_tab` | $(\partial\rho/\partial h)_p$ | kg²/(m³·J) |
| `hT_tab` | $c_p = (\partial h/\partial T)_p$ | J/(kg·K) |
| `s_tab2D` | Specific entropy | J/(kg·K) |
| `rp_tab` | $(\partial\rho/\partial p)_h$ | kg/(m³·Pa) |
| `sound_tab` | Speed of sound | m/s |
| `mi_tab2D` | Dynamic viscosity | Pa·s |
| `k_tab2D` | Thermal conductivity | W/(m·K) |

---

## Ideal-gas properties

Module `FLINT_Lib_Thermodynamic`. Every procedure is `pure`.

| Procedure | Returns |
|---|---|
| `f_Rtot(rhoi)` | Gas constant of the mixture, J/(kg·K) |
| `f_molecularWeight(rhoi)` | Molecular weight of the mixture, kg/kmol |
| `call co_rotot_Rtot(rhoi, rho, R)` | `rho = sum(rhoi)` and `R = f_Rtot(rhoi)` |
| `EOS(p, rho, R)` | Temperature of the perfect gas, `p/(rho*R)`, K |
| `f_cp(rhoi, T, rho)` | $c_p$ of the mixture, J/(kg·K) |
| `f_gamma(rhoi, p, rho, Rtot)` | Ratio of specific heats $\gamma$ |
| `f_ss(rhoi, p, rho, Rtot)` | Speed of sound, m/s |
| `H(p, rhoi)`, `H0(p, rhoi, u)` | Static and total ($+u^2/2$) specific enthalpy of the mixture, J/kg |
| `E(p, rhoi)`, `E0(p, rhoi, u)` | Static ($H - p/\rho$) and total specific internal energy, J/kg |
| `cpi(i, T)`, `H0i(i, T)`, `S0i(i, T)` | $c_p$, $h$ and $s$ of species `i` |
| `f_tabT(T, i, tab)` | Any `(Tmin:Tmax, ns)` table of species `i` at `T` |
| `mass2molar(Y)`, `molar2mass(X)` | Mole fractions from mass fractions, and the converse |
| `f_laminarViscosity(rhoi, p, rho, Rtot)` | Viscosity of the mixture (Wilke), Pa·s |
| `call co_k_mi_lam_Wilke(rhoi, rho, T, mu, k)` | Viscosity and thermal conductivity of the mixture (Wilke) |
| `call co_DS(rhoi, p, rho, Rtot, Dm)` | Mixture-averaged diffusion coefficients `Dm(ns)`, m²/s: $(1-X_s)/\sum_{j\neq s} X_j/D_{sj}$, scaled from `dij_pref` to `p` |
| `prim2cons(prim)`, `cons2prim(cons, T)` | Conservative from primitive variables, and the converse (see below) |
| `call check_gas_state(rhoi, p, error)` | Raises every `rhoi` to at least 1e-20; `error` = 1 for `p < 0`, 2 for a temperature out of `Tmin..Tmax`, 0 otherwise |

The primitive variables are `(rhoi(1:ns), u(1:nd), p)` and the conservative ones `(rhoi(1:ns), rho*u(1:nd), rho*E0)`,
for `nd` velocity components. `cons2prim` finds the temperature by Newton iterations on the energy from the guess
`T` (relative tolerance 1e-6, at most 20 iterations).

### Real fluid

| Procedure | Returns |
|---|---|
| `ph2vars(p, h, tab2D)` | The value of a real-fluid table (`rho_tab`, `T_tab`, …) at $(p, h)$, by bilinear interpolation; $(p, h)$ must lie on the grid |
| `pT2h(p, T)` | The specific enthalpy at $(p, T)$, J/kg; needs `ph2pT()` |

---

## Finite-rate chemistry

Module `FLINT_Lib_Chemistry_wdot`.

`call Assign_Mechanism(name)` selects the procedure that computes the source terms. A name with a compiled routine
selects that routine, and its analytical Jacobian when it has one:

`WD`, `WD-Andersen`, `JLR-Nasuti`, `Frassoldati`, `CKJLR-10sp`, `Smooke`, `CORIA-CNRS`, `TSR-CDF-13`, `TSR-GP-24`,
`TSR-Rich-31`, `ZK`, `CoronettiC4H6`, `Singh`, `Singh-WC32`, `Frolov`, `Frolov_nopressure`, `Nassini`, `Troyes`,
`Ecker`, `Cross`, `Pelucchi`, `OSK`, `ONERA-7`, `SanDiego`, `FFCMy-12`, `Gerlinger-9`.

For these names the loaded species and rate tables must match the routine: `Assign_Mechanism` checks them and stops
with the expected and the loaded lists when they differ. Called before the data are loaded, it checks them at the
first call of `chemistry_source`. Any other name selects the general procedure, which evaluates the mechanism of
`chemistry-info.txt`, and prints a WARNING; in strict mode it stops instead. Strict mode is set by the module
variable `FLINT_strict_mechanism` (`logical`, default `.false.`) or by the environment variable
`FLINT_STRICT_MECHANISM` set to `1`, `true`, `yes` or `on`, in any letter case.

| Procedure pointer | Arguments |
|---|---|
| `call chemistry_source(roi, temp, omegadot)` | `roi(ns)`, partial densities, kg/m³ (inout); `temp`, K (in); `omegadot(ns)`, net mass production rates, kg/(m³·s) (out) |
| `call chemistry_jacobian(roi, temp, dwdr, dwdT)` | `dwdr(i,j)` = $\partial\dot\omega_i/\partial\rho_j$, `dwdT(i)` = $\partial\dot\omega_i/\partial T$ (out). Null when the mechanism has no analytical Jacobian |

### Constant-volume reactor

Module `FLINT_Lib_Chemistry_rhs`. The state of an adiabatic constant-volume reactor is `Z = (rhoi(1:ns), T)`, with
`nz = ns + 1`, and the signatures are those of the ODE solvers of OSLO.

| Procedure | Description |
|---|---|
| `call rhs_native(nz, time, Z, F)` | `F = (omegadot(1:ns), dT/dt)`. Every `F` is −1 when `T` is NaN or out of the thermo or rate tables |
| `call jac_native(nz, time, Z, DFY, LDFY, RPAR, IPAR)` | Analytical Jacobian of `rhs_native`, from `chemistry_jacobian` |
| `call set_analytical_jacobian(flag)` | Sets `IJAC_chem` to 1, so that the ODE solver calls `jac_native`, when `flag` is true and the mechanism has an analytical Jacobian; otherwise 0 (finite differences, with a message if `flag` is true). Call it after `Assign_Mechanism` |
| `call rhs_cantera(nz, time, Z, F)` | As `rhs_native`, with the rates of Cantera's phase `gas` (module variable; `USE_CANTERA` only) |

The rate coefficients of the tables are read by `f_kf`, `f_kb` (module `FLINT_Lib_Chemistry_data`), `f_k_troe` and
`f_k_lindemann` (module `FLINT_Lib_Chemistry_falloff`), which the mechanism routines call; see
[Generating New Mechanisms](chemistry_routines.md).

---

## Chemical equilibrium

Module `FLINT_CEA_solver`. It needs the thermo tables and `composition.txt` (`read_idealgas_thermo` returning 0).

`call CEA_solve(T_in, rhoi_in, T_out, rhoi_out, [status])` computes the equilibrium at constant internal energy and
volume of the state `(T_in, rhoi_in)`: `T_out`, K, and `rhoi_out(ns)`, the equilibrium **mass fractions**. It
initializes the solver (`CEA_initialize_global`, module `FLINT_CEA_setup`) at its first call in each thread.

`status` (`integer`, out, optional) is 0 when the iteration converged and 1 when it did not, for instance for a state
whose energy lies beyond the thermo tables. Without convergence `T_out = T_in`, and `rhoi_out` holds the mass
fractions of the input state, `rhoi_in/sum(rhoi_in)`. Each such call also writes one line, to standard output and to
standard error, with the input state:

```text
[WARNING] FLINT CEA_solve: no equilibrium for T_in =  3.0000E+04 K, rho =  3.2500E+00 kg/m3: returning T_in and the input mass fractions (status 1)
```

---

## Cantera

Module `FLINT_cantera_load`, built with `USE_CANTERA` only.

`call load_phase(phase, filename, [phase_name])` imports the phase `phase_name` of the Cantera file `filename` into
`phase` (`type(phase_t)`, from Cantera's Fortran interface).
