# Model equations

The model couples six fields on a 3D cube of normalised side length `L = 1`:

| Symbol | Field | Role |
|---|---|---|
| S | Hyaluronic-acid (HA) scaffold fraction | Material that is degraded; 1 = intact |
| M | MMP-9 concentration | Enzyme that degrades the scaffold and native ECM |
| T | TIMP-1 concentration | Inhibitor released from depots to neutralise MMP |
| E | New extracellular matrix (ECM) density | Deposited by progenitor cells to replace the scaffold |
| P | Neural progenitor cell (NPC) density | Proliferates and differentiates |
| N | Mature neurons (one global number) | Functional output of the model |

The 84-day script adds two more fields, astrocytes (AST) and endothelial cells (ECL), and a BMP schedule (see below).

All quantities are normalised or in arbitrary units. Time is in **weeks** in scripts 01-03 and in **days** in script 04.

## Core equations (as presented in the final presentation)

**1. Scaffold degradation (S)**

$$\frac{\partial S}{\partial t} = D_s \nabla^2 S - k_d\, M\, S$$

**2. MMP dynamics (M)**

$$\frac{\partial M}{\partial t} = D_m \nabla^2 M + \alpha_P P + \alpha_A - \beta_M M - \gamma_T\, T\, M$$

**3. TIMP release (T)**

$$\frac{\partial T}{\partial t} = D_T \nabla^2 T + R_0 e^{-\mu t}\,\mathrm{Depot}(\mathbf{x}) - \lambda_T T$$

**4. ECM deposition (E)**

$$\frac{\partial E}{\partial t} = k_E P - k_{mE}\, M\, E$$

**5. NPC proliferation and differentiation (P)**

$$\frac{\partial P}{\partial t} = r_P P \left(1 - \frac{P}{P_{max}}\right) - k_{diff}\, S\, P + \eta\, W(\mathbf{x})$$

**6. Neuronal maturation (N)**

$$\frac{dN}{dt} = k_{diff} \int_{\Omega} S(\mathbf{x},t)\,P(\mathbf{x},t)\,dV - \delta_N N$$

`Depot(x)` is a sum of Gaussian sources, normalised to a maximum of 1. `W(x)` is a static Gaussian Wnt field centred in the domain, normalised to a maximum of 1.

## Spatially varying diffusion

Each field has its own diffusion range. The local coefficient depends on the local depot density:

$$D_{field}(\mathbf{x}) = D_{min} + (D_{max} - D_{min})\cdot \mathrm{Depot}_{field}(\mathbf{x})$$

- Script 02/03 (25 depots): only MMP uses a spatially varying `D_m`; the other diffusion coefficients are constants.
- Script 04 (1000 depots): all seven fields (T, M, S, E, P, AST, ECL) use spatially varying coefficients, and the divergence form is expanded as `D ∇²u + ∇D · ∇u`.

## Changes in the 84-day model (script 04)

| Equation | Core version | 84-day version |
|---|---|---|
| S | `- k_d M S` | `- k_d min(M,1) S`, with reduced `k_d = 0.05` |
| P | `+ η W` | `+ k_w_p W - 0.5 B(t) P` (Wnt-driven growth, BMP suppresses proliferation) |
| N | `δ_N` decay, integral of `S·P` | integral of `S·P` scaled by `Ng^3`, `k_diff = 1.0` |
| AST | not present | `∂AST/∂t = ∇·(D_ast ∇AST) + r_ast · Depot_AST - k_ast · M · AST` |
| ECL | not present | `∂ECL/∂t = ∇·(D_ecl ∇ECL) + r_ecl · Depot_ECL - k_ecl · M · ECL` |

BMP schedule `B(t)`: constant `BMP_base = 0.02` until day 60, then a linear ramp to `BMP_final = 0.18` at day 84.

Initial conditions in script 04: `S = 1`, `T = 0`, `E = 0`, `M = 0.4 + noise (SD 0.05)`, `P = 0.05 + 0.2·exp(-r²/0.08)` (central NPC cluster, capped at `P_max`), `AST = P/50`, `ECL = P/5`, both floored at `1e-6`.

## Lumped dashboard model (script 01)

A week-by-week Euler model of the same system without space. NPCs grow logistically and are modulated by `(Wnt - BMP)`. Neurons form from NPCs at a rate proportional to `(BDNF + 0.2)`. Scaffold loss is `-k_d · S · M / (T + 0.1)`. MMP and TIMP follow simple production/decay terms, and the morphogens (Wnt, BMP, BDNF, laminin-511) follow prescribed week-by-week update rules (Wnt and laminin decline slowly, BMP and BDNF rise).

## Numerical notes

- Time stepping is explicit Euler. Scripts 02/03 use `dt = 0.01`, script 04 uses `dt = 0.1 day` (10 substeps per day).
- Spatial derivatives use MATLAB's `del2`. For a 3D array, `del2` returns **one sixth** of the Laplacian (`∇²/6`). The effective diffusivities in the simulations are therefore the listed `D` values divided by 6. Multiply `del2(...)` by 6 if the nominal `D` should apply exactly. This was left unchanged so that the figures in `results/` can be reproduced.
- Scripts 02/03 clamp `S` to [0, 1] and `P` to [0, `P_max`]; script 04 does the same and floors the other fields at 0.
