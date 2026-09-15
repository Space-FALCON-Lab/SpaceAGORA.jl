---
id: gnc.switch_window_solver_switch_calculation
label: switch_calculation
kind: function
source:
  file: src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl
  symbol: switch_calculation
  lines:
  - 64
  - 64
inputs:
- id: ip
  type: Any
  units: n/a
  required: true
  description: Positional argument `ip`.
- id: m
  type: Any
  units: n/a
  required: true
  description: Positional argument `m`.
- id: position
  type: Any
  units: n/a
  required: true
  description: Positional argument `position`.
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: t
  type: Any
  units: n/a
  required: true
  description: Positional argument `t`.
- id: heat_rate_control
  type: Any
  units: n/a
  required: true
  description: Positional argument `heat_rate_control`.
- id: reevaluation_mode
  type: Any
  units: n/a
  required: true
  description: Positional argument `reevaluation_mode`.
- id: current_position
  type: Any
  units: n/a
  required: false
  description: Positional argument `current_position` (default `0`).
- id: cnf
  type: Any
  units: n/a
  required: false
  description: Keyword argument `cnf` (default `nothing`).
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: result
  type: AbstractArray
  units: n/a
  description: Return value of `switch_calculation`. Returns `[0.0, 0.0]` or `[0.0,
    t_cf[end]/2]` or `t_switch`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# switch_calculation

## Purpose
Closed-form solver for the two angle-of-attack switch times of an aerobraking drag passage. It avoids the full trajectory integration used by the integration-based variant by propagating an analytic approximation of the pass, then locating the arc over which the costate favours the high-drag attitude.

## Theory & Math
Molecular speed ratio $S = v/\sqrt{2 R T}$ with $v$ the relative velocity (m/s), $R$ the planetary specific gas constant (J/kg/K) and $T$ the atmospheric temperature (K). The drag coefficient is linearised in angle of attack $\alpha$ as $C_D(\alpha) \approx C_{D,0} + \alpha\,(C_{D,90} - C_{D,0})/(\pi/2)$, with $C_{D,90}$ and $C_{D,0}$ the free-molecular values at $\alpha = \pi/2$ and $\alpha = 0$. The switch is placed where the switching function $\lambda_v(t)$ crosses the threshold $\lambda_{switch}$, the bang-bang condition of the underlying optimal-control problem.

## Design & Implementation
It first calls `closed_form(args, m, position, T, true, m.aerodynamics.α)` with the planet temperature `T = m.planet.T`, obtaining approximate time, altitude, flight-path-angle and velocity histories `t_cf, h_cf, γ_cf, v_cf`. The molecular speed ratio is formed as `S = v_cf / sqrt(2*R*T)`, and free-molecular coefficients are evaluated at incidence `pi/2` and `0` via `aerodynamic_coefficient_fM`, yielding the linearised triple `coeff = (CD_slope, CL_0, CD_0)` with `CD_slope = (CD_90 - CD_0)/(pi/2)`. Bracket residuals `delta_Q_max` and `delta_Q_min` come from `func` at gains `0.1` and `0.0`; when they straddle zero, `fzero(..., [0.0, 0.1], Roots.Brent())` finds the gain. If the upper residual is negative it returns `[0.0, 0.0]`; if the lower residual is positive it returns `[0.0, t_cf[end]/2]`. With the root gain it re-evaluates the profile, computes `lambda_switch` and `lambdav` from `lambdas`, takes the mask `lambdav .< lambda_switch`, multiplies it into `t_cf`, filters out zeros, and takes the first and last surviving times. The window is then shrunk: if the end is within 5 s of the pass end the opening time drops by 4% of the pass duration, within 60 s it drops by 5 s, and the closing time is always reduced by 10%.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `ip` | Any | n/a | yes | Positional argument `ip`. |
| in | `m` | Any | n/a | yes | Positional argument `m`. |
| in | `position` | Any | n/a | yes | Positional argument `position`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `heat_rate_control` | Any | n/a | yes | Positional argument `heat_rate_control`. |
| in | `reevaluation_mode` | Any | n/a | yes | Positional argument `reevaluation_mode`. |
| in | `current_position` | Any | n/a | no | Positional argument `current_position` (default `0`). |
| in | `cnf` | Any | n/a | no | Keyword argument `cnf` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractArray | n/a | — | Return value of `switch_calculation`. Returns `[0.0, 0.0]` or `[0.0, t_cf[end]/2]` or `t_switch`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.e_edg_strategy_compute_e_edg_guidance_window_bang|compute_e_edg_guidance_window!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/e_edg/e_edg_strategy.jl:32-32`

**Downstream**

- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:75-75`
- `callees` → [[gnc.heat_rate_models_aoa|aoa]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:83-83`
- `callees` → [[gnc.heat_rate_models_func|func]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:87-87`
- `callees` → [[gnc.second_switch_solver_func|func]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:87-87`
- `callees` → [[gnc.switch_window_solver_func|func]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:87-87`
- `callees` → [[gncx.closed_form_solution_closed_form|closed_form]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:70-70`
- `callees` → [[gncx.heat_rate_models_lambdas|lambdas]] · `callers` · call · `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl:100-100`
<!-- vulcan:connections:end -->

## Limitations
`filter(!iszero, temp)` silently drops any genuine sample at `t_cf == 0`, and if the mask selects nothing the subsequent `temp[1]` raises a `BoundsError`. The 4%, 5 s and 10% trims are unexplained magic numbers with no physical derivation and can push `t_switch[1]` negative on a short pass. The Brent bracket `[0.0, 0.1]` is fixed, so a problem whose root lies above 0.1 is never found, and the whole result inherits the accuracy of the closed-form pass approximation rather than the integrated dynamics.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/e_edg/switch_window_solver.jl` line 64.
