---
id: gnc.targeting_solver_func_e
label: func_e
kind: function
source:
  file: src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl
  symbol: func_e
  lines:
  - 348
  - 348
inputs:
- id: nu_E_var
  type: Any
  units: n/a
  required: true
  description: Positional argument `nu_E_var`.
- id: m
  type: Any
  units: n/a
  required: true
  description: Positional argument `m`.
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: coeff
  type: Any
  units: n/a
  required: true
  description: Positional argument `coeff`.
- id: position
  type: Any
  units: n/a
  required: true
  description: Positional argument `position`.
- id: heat_rate_control
  type: Any
  units: n/a
  required: true
  description: Positional argument `heat_rate_control`.
- id: approx_sol
  type: Any
  units: n/a
  required: true
  description: Positional argument `approx_sol`.
- id: energy_target
  type: Any
  units: n/a
  required: true
  description: Positional argument `energy_target`.
- id: initial_guess
  type: Any
  units: n/a
  required: false
  description: Positional argument `initial_guess` (default `false`).
- id: approx_calc
  type: Any
  units: n/a
  required: false
  description: Positional argument `approx_calc` (default `false`).
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
  description: Return value of `func_e`. Returns `in_cond_lambda` or `t_cf, v_cf,
    γ_cf, h_cf` or `delta_E`.
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

# func_e

## Purpose
Iteratively evaluates the closed-form aerobraking pass under a candidate energy multiplier `nu_E_var`, converging the angle-of-attack profile and trajectory to a self-consistent fixed point, and returns the energy miss `E_fin - energy_target`. It is the residual function that Brent's method drives to zero in `control_solarpanels_targeting_closed_form`, and doubles as a trajectory/initial-guess provider via its flag arguments.

## Theory & Math
Fixed-point iteration on the pass energy: $\varepsilon_k = \tfrac{1}{2} v_{f,k}^2 - \mu/(R_{p} + h_{f,k})$, stopping when $\varepsilon_k - \varepsilon_{k-1} \le 10^{-3}$ J/kg, where $v_{f,k}$ and $h_{f,k}$ are the final closed-form speed and altitude at iteration $k$ and $R_p$ is `m.planet.Rp_e`. Speed ratio used by the free-molecular heating model: $S = \sqrt{\gamma/2}\,M$, $M = v/\sqrt{\gamma R T}$, with $\gamma$ = `m.planet.γ`, $R$ = `m.planet.R`, $T$ = `m.planet.T`.

## Design & Implementation
Reads `multiplicative_factor_heatload` and `max_heat_rate` from `args` via `_bridge_required_field`, unpacks `approx_sol` into `(t_cf, h_cf, γ_cf, v_cf)`, and loops `while temp_E_diff > 1e-3`. Each iteration calls `aoa(m, k_cf, t_cf, h_cf, γ_cf, v_cf, coeff, nu_E_var, aoa_cf)` to update the angle-of-attack profile and multiplier initial condition `in_cond_lambda`, re-runs `closed_form(args, m, position, T, true, m.aerodynamics.α, aoa_cf)`, computes the speed ratio `S = sqrt(γ/2)·M` with `M = v_cf / sqrt(γ R T)`, exponential density `density_exp(h_cf, m.planet)[1]`, and `heat_rate_calc(...)`, clipping to `max_heat_rate` when `heat_rate_control` is true. The final energy `E_fin = |v_end|^2/2 - μ/(Rp_e + h_end)` is compared with the previous iterate to update `temp_E_diff`. When `initial_guess` is true the function returns `in_cond_lambda`; when `approx_calc` is true it returns `(t_cf, v_cf, γ_cf, h_cf)`; otherwise it returns `delta_E`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `nu_E_var` | Any | n/a | yes | Positional argument `nu_E_var`. |
| in | `m` | Any | n/a | yes | Positional argument `m`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `coeff` | Any | n/a | yes | Positional argument `coeff`. |
| in | `position` | Any | n/a | yes | Positional argument `position`. |
| in | `heat_rate_control` | Any | n/a | yes | Positional argument `heat_rate_control`. |
| in | `approx_sol` | Any | n/a | yes | Positional argument `approx_sol`. |
| in | `energy_target` | Any | n/a | yes | Positional argument `energy_target`. |
| in | `initial_guess` | Any | n/a | no | Positional argument `initial_guess` (default `false`). |
| in | `approx_calc` | Any | n/a | no | Positional argument `approx_calc` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `func_e`. Returns `in_cond_lambda` or `t_cf, v_cf, γ_cf, h_cf` or `delta_E`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_solver_control_solarpanels_targeting_closed_form|control_solarpanels_targeting_closed_form]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:318-318`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl`

**Downstream**

- `callees` → [[gnc.bridge_helpers__bridge_required_field|_bridge_required_field]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:349-349`
- `callees` → [[gnc.heat_rate_models_aoa|aoa]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:367-367`
- `callees` → [[gnc.tracking_executor_heat_rate_calc|heat_rate_calc]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:375-375`
- `callees` → [[gncx.closed_form_solution_closed_form|closed_form]] · `callers` · call · `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl:368-368`
<!-- vulcan:connections:end -->

## Limitations
The convergence test is signed (`E_fin - E_prev > 1e-3`), so a sequence whose energy decreases between iterations terminates after one step regardless of magnitude, and an oscillating sequence can loop forever since there is no iteration cap (`count` is incremented but never checked). `heat_rate` is computed and clipped but not used afterwards, so `heat_rate_control` has no effect on the returned values. `k_cf = 1.0` is hard-coded. The multiple return shapes (scalar, vector, 4-tuple) selected by boolean flags make the function type-unstable and easy to misuse.

## Provenance
Mapped from `src/gnc/guidance/aerobraking/t_edg/targeting_solver.jl` line 348.
