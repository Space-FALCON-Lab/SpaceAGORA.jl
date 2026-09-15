---
id: gnc.heat_load_control__edg_heat_load_profile_for_k
label: _edg_heat_load_profile_for_k
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_heat_load_profile_for_k
  lines:
  - 542
  - 542
inputs:
- id: config
  type: Any
  units: n/a
  required: true
  description: Positional argument `config`.
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: spacecraft
  type: Any
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
- id: pos
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos`.
- id: vel
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel`.
- id: mass
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mass`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: env
  type: Any
  units: n/a
  required: true
  description: Positional argument `env`.
- id: coeffs
  type: Any
  units: n/a
  required: true
  description: Positional argument `coeffs`.
- id: k
  type: Float64
  units: n/a
  required: true
  description: Positional argument `k`.
- id: controlled_panel_links
  type: Tuple{Vararg{Int}}
  units: n/a
  required: true
  description: Positional argument `controlled_panel_links`.
- id: heat_rate_control
  type: Bool
  units: n/a
  required: true
  description: Positional argument `heat_rate_control`.
- id: structural_control
  type: Bool
  units: n/a
  required: true
  description: Positional argument `structural_control`.
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
  description: Return value of `_edg_heat_load_profile_for_k`. Returns `track, alpha_profile`.
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

# _edg_heat_load_profile_for_k

## Purpose
For one candidate heating weight `k`, produces the predicted track and the constrained bang-bang alpha profile, iterating between costate solution, constraint enforcement, and re-propagation when the TPBVP solver is enabled.

## Design & Implementation
Signature `(config, p, spacecraft, pos, vel, mass, t, env, coeffs, k, controlled_panel_links, heat_rate_control, structural_control)`. It computes `area`, `scale_height`, and the closed-form `base_track`, seeding `alpha_profile` at `max_alpha_rad`. In `:tpbvp_integration` mode it computes `altitude_offset_m = env.altitude_m - (norm(pos) - Rp_e)`, constrains the seed, integrates the track with RK4, then performs two iterations of (costate profile via `_edg_heat_load_alpha_profile`, constrain, re-integrate) followed by one final costate-plus-constrain step, returning `(track, alpha_profile)`. In the closed-form mode it keeps `base_track` and applies the costate map three times without constraints.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | Any | n/a | yes | Positional argument `config`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `pos` | SVector{3, Float64} | n/a | yes | Positional argument `pos`. |
| in | `vel` | SVector{3, Float64} | n/a | yes | Positional argument `vel`. |
| in | `mass` | Float64 | n/a | yes | Positional argument `mass`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `env` | Any | n/a | yes | Positional argument `env`. |
| in | `coeffs` | Any | n/a | yes | Positional argument `coeffs`. |
| in | `k` | Float64 | n/a | yes | Positional argument `k`. |
| in | `controlled_panel_links` | Tuple{Vararg{Int}} | n/a | yes | Positional argument `controlled_panel_links`. |
| in | `heat_rate_control` | Bool | n/a | yes | Positional argument `heat_rate_control`. |
| in | `structural_control` | Bool | n/a | yes | Positional argument `structural_control`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_heat_load_profile_for_k`. Returns `track, alpha_profile`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control_residual|residual]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:649-649`
- [[gncx.heat_load_control__edg_solve_heat_load_switches|_edg_solve_heat_load_switches]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:649-649`

**Downstream**

- `callees` → [[gnc.heat_load_control__edg_closed_form_heat_load_trajectory|_edg_closed_form_heat_load_trajectory]] · `callers` · call · `src/gnc/control/heat_load_control.jl:560-560`
- `callees` → [[gnc.heat_load_control__edg_constrained_heat_load_alpha_profile|_edg_constrained_heat_load_alpha_profile]] · `callers` · call · `src/gnc/control/heat_load_control.jl:565-565`
- `callees` → [[gnc.heat_load_control__edg_heat_load_alpha_profile|_edg_heat_load_alpha_profile]] · `callers` · call · `src/gnc/control/heat_load_control.jl:579-579`
- `callees` → [[gnc.heat_load_control__edg_heat_load_scale_height|_edg_heat_load_scale_height]] · `callers` · call · `src/gnc/control/heat_load_control.jl:559-559`
- `callees` → [[gnc.heat_load_control__edg_integrated_heat_load_trajectory|_edg_integrated_heat_load_trajectory]] · `callers` · call · `src/gnc/control/heat_load_control.jl:575-575`
- `callees` → [[gnc.heat_load_control__edg_total_ref_area|_edg_total_ref_area]] · `callers` · call · `src/gnc/control/heat_load_control.jl:558-558`
<!-- vulcan:connections:end -->

## Limitations
Iteration counts (2 and 3) are hard-coded with no convergence test, and each TPBVP iteration re-runs a full RK4 propagation with panel aerodynamics at every stage, so a single `k` evaluation is expensive. The closed-form branch returns an unconstrained profile; constraints are applied later by the caller. `altitude_offset_m` is treated as constant for the whole pass.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 542.
