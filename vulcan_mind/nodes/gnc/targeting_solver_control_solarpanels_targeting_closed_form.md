---
id: gnc.targeting_solver_control_solarpanels_targeting_closed_form
label: control_solarpanels_targeting_closed_form
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl
  symbol: control_solarpanels_targeting_closed_form
  lines:
  - 290
  - 290
inputs:
- id: energy_target
  type: Any
  units: n/a
  required: true
  description: Positional argument `energy_target`.
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
  type: Any
  units: n/a
  description: Return value of `control_solarpanels_targeting_closed_form`. Returns
    `t_switch`.
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

# control_solarpanels_targeting_closed_form

## Purpose
Estimates the pair of switch times `[t_on, t_off]` for solar-panel drag modulation using the closed-form aerobraking approximation and a Lagrange-multiplier switching law, instead of repeated numerical integration. It is intended as a fast initial guess or replanning path for the T-EDG aerobraking guidance.

## Design & Implementation
Calls `closed_form(args, m, position, T, true, m.aerodynamics.α)` with `T = m.planet.T` to obtain `(t_cf, h_cf, γ_cf, v_cf)`, computes the speed ratio `S = v_cf / sqrt(2 R T)` and evaluates free-molecular coefficients at α = π/2 and α = 0 via `aerodynamic_coefficient_fM`, forming `coeff = (CD_slope, CL_0, CD_0)` with `CD_slope = (CD_90 - CD_0)/(π/2)`. The multiplier `nu_E_root` is found with `fzero(nu_E -> func_e(nu_E, ...), [1, 100], Roots.Brent())`. A second `func_e` call with `approx_calc=true` returns the converged closed-form trajectory, then `lambdas(...)` yields `lambda_switch` and the multiplier vector `lambdav`; times where `lambdav .< lambda_switch` are collected and the first and last non-zero entries form `t_switch`. Heuristic adjustments follow: if the last switch is within 5 s of `t_cf[end]` the first switch is pulled earlier by 4 % of the pass duration, else if within 60 s by 5 s; the second switch is always reduced by 10 %.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `energy_target` | Any | n/a | yes | Positional argument `energy_target`. |
| in | `ip` | Any | n/a | yes | Positional argument `ip`. |
| in | `m` | Any | n/a | yes | Positional argument `m`. |
| in | `position` | Any | n/a | yes | Positional argument `position`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `heat_rate_control` | Any | n/a | yes | Positional argument `heat_rate_control`. |
| in | `reevaluation_mode` | Any | n/a | yes | Positional argument `reevaluation_mode`. |
| in | `current_position` | Any | n/a | no | Positional argument `current_position` (default `0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `control_solarpanels_targeting_closed_form`. Returns `t_switch`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl`

**Downstream**

- `callees` → [[dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_fm|aerodynamic_coefficient_fM]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:300-300`
- `callees` → [[gnc.targeting_solver_func_e|func_e]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:318-318`
- `callees` → [[gncx.closed_form_solution_closed_form|closed_form]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:295-295`
- `callees` → [[gncx.heat_rate_models_lambdas|lambdas]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:327-327`
<!-- vulcan:connections:end -->

## Limitations
The call to `lambdas(m, aoa_cf, ...)` references `aoa_cf`, which is only assigned in a commented-out line above, so as written this function raises `UndefVarError` when reached. The arguments `ip`, `t`, `reevaluation_mode` and `current_position` are accepted but never used. `k_cf = 1.0`, the bracket `[1, 100]`, and the 5 s / 60 s / 4 % / 10 % adjustments are unexplained magic numbers. `filter(!iszero, temp)` discards a genuine switch at t = 0 and throws `BoundsError` if no time satisfies the switching condition.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl` line 290.
