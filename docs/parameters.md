# Parameters

Values are taken directly from the scripts in `code/`. They are model assumptions in normalised units; this repository does not include calibration against experimental data.

## Script 04: 84-day model (22^3 grid, 1000 depots)

### Grid, time and depots

| Parameter | Value | Description |
|---|---|---|
| `Ng` | 22 | Grid points per spatial dimension |
| `L` | 1.0 | Domain length (normalised) |
| `dx` | `L/(Ng-1)` | Grid spacing |
| `t_final` / `nDays` | 84 | Simulated days (12 weeks), one image per day |
| `substeps_per_day` | 10 | Euler sub-steps per day |
| `dt_day` | 0.1 | Time step (days) |
| `numDepots` | 1000 | Gaussian depot sources |
| `sigma` | `0.02 L` = 0.02 | Standard deviation of each depot |

### Diffusion coefficient ranges and depot amplitudes

| Field | `D_min` | `D_max` | Amplitude variable | Amplitude range |
|---|---|---|---|---|
| TIMP (T) | 0.004 | 0.02 | `amps_T` | 0.5 - 2.0 |
| MMP (M) | 0.003 | 0.015 | `amps_M` | 0.7 - 2.3 |
| Scaffold (S) | 0.001 | 0.006 | `amps_S` | 0.3 - 1.5 |
| ECM (E) | 0.002 | 0.01 | `amps_E` | 0.4 - 1.4 |
| NPC (P) | 0.002 | 0.008 | `amps_P` | 0.4 - 1.4 |
| Astrocytes (AST) | 0.0008 | 0.006 | `amps_AST` | 0.8 - 2.0 |
| Endothelial (ECL) | 0.001 | 0.007 | `amps_ECL` | 0.8 - 2.0 |

### Kinetic parameters

| Parameter | Value | Meaning |
|---|---|---|
| `R0` | 1.0 | TIMP source coefficient |
| `mu` | 0.25 | TIMP release decay rate |
| `lambdaT` | 0.10 | TIMP removal rate |
| `alpha_P` | 0.3 | NPC-driven MMP production |
| `alpha_A` | 0.25 | Basal MMP production |
| `beta_M` | 0.20 | MMP decay |
| `gamma_T` | 0.40 | TIMP-MMP inhibition |
| `k_d` | 0.05 | Scaffold degradation (reduced) |
| `k_E` | 0.10 | ECM production by NPCs |
| `k_mE` | 0.05 | ECM degradation by MMP |
| `r_P` | 0.50 | NPC proliferation |
| `P_max` | 0.8 | Maximum NPC density at a point |
| `k_diff` | 1.0 | Neuron formation (increased) |
| `delta_N` | 0.01 | Neuron decay |
| `r_ast` | 0.20 | Astrocyte formation rate |
| `r_ecl` | 0.25 | Endothelial formation rate |
| `k_ast` | 0.10 | Astrocyte loss by MMP |
| `k_ecl` | 0.12 | Endothelial loss by MMP |
| `k_w_p` | 0.08 | Wnt-driven NPC proliferation |
| `BMP_base` | 0.02 | BMP value before the ramp |
| `BMP_final` | 0.18 | BMP value at day 84 |
| `bmp_ramp_start_day` / `bmp_ramp_end_day` | 60 / 84 | BMP ramp window |

### Cell seeding

| Cell type | Depots / points in the code | Design ratio | Function |
|---|---|---|---|
| NPC | present throughout the scaffold (initial central cluster) | 50 | Progenitor / neurogenesis |
| Endothelial | 164 depots (`round(1000 * 10/61)`) | 10 | Vascularisation |
| Astrocyte | 16 depots (`round(1000 * 1/61)`) | 1 | Glial support / ECM |

## Scripts 02 and 03: 12-week spatial model (25^3 grid, 25 depots)

| Parameter | Value | | Parameter | Value |
|---|---|---|---|---|
| `Ng` | 25 | | `alpha_P` | 0.3 |
| `L` | 1.0 | | `alpha_A` | 0.25 |
| `t_final` | 12 weeks | | `beta_M` | 0.20 |
| `dt` | 0.01 | | `gamma_T` | 0.40 |
| `numDepots` | 25 | | `D_s` | 0.002 |
| `sigma` | `0.04 L` | | `k_d` | 0.15 |
| `D_T` | 0.01 | | `k_E` | 0.10 |
| `R0` | 1.0 | | `k_mE` | 0.05 |
| `mu` | 0.25 | | `r_P` | 0.50 |
| `lambdaT` | 0.10 | | `P_max` | 0.8 |
| `D_m_background` | 0.003 | | `k_diff` | 0.05 |
| `D_m_depot_max` | 0.015 | | `eta` | 0.25 |
| | | | `delta_N` | 0.01 |

Initial conditions: `T = 0`, `M = 0.8 + noise (SD 0.05)`, `S = 1`, `E = 0`, `P = 0.05 + 0.20 exp(-r^2/0.08)`.

## Script 01: lumped 12-week dashboard model

| Group | Parameters |
|---|---|
| NPC and neurons | `rP = 0.45`, `deltaN = 0.18`, carrying capacity 1, `P(0) = 0.05`, `N(0) = 0` |
| Scaffold | `S0 = 1`, `kd = 0.08` |
| MMP / TIMP | `M0 = 0.2`, `T0 = 0.5`, `alphaP = 0.25`, `alphaA = 0.15`, `betaM = 0.12`, `gammaT = 0.5`, `R0 = 0.3`, `lambdaT = 0.05` |
| Morphogens (initial values) | Wnt 0.9 (CHIR99021), BMP 0.4, BDNF 0.3, Laminin-511 0.8 |

`Ds = 0.01` is defined in the script but not used by the lumped model.
