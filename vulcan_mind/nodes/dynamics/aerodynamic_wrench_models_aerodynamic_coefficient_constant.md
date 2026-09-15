---
id: dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_constant
label: aerodynamic_coefficient_constant
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: aerodynamic_coefficient_constant
  lines:
  - 1005
  - 1005
inputs:
- id: alpha
  type: Any
  units: n/a
  required: true
  description: Positional argument `α`.
- id: body
  type: Any
  units: n/a
  required: true
  description: Positional argument `body`.
- id: T
  type: Any
  units: n/a
  required: true
  description: Positional argument `T`.
- id: S
  type: Any
  units: n/a
  required: true
  description: Positional argument `S`.
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: montecarlo
  type: Any
  units: n/a
  required: false
  description: Positional argument `montecarlo` (default `false`).
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
  description: Return value of `aerodynamic_coefficient_constant`. Returns `CL_body,
    CD_body`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# aerodynamic_coefficient_constant

## Purpose
Legacy constant-model coefficient function returning zero lift and the linear CD law, with an optional Monte Carlo perturbation hook.

## Design & Implementation
Signature `(α, body, T, S, args, montecarlo=false)`. Computes `CL_body = 0.0` and `CD_body = 2*(2.2-0.8)/π * args.α + 0.8`, then if `montecarlo == true` calls `monte_carlo_aerodynamics(CL_body, CD_body, args)`. Returns `(CL_body, CD_body)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `alpha` | Any | n/a | yes | Positional argument `α`. |
| in | `body` | Any | n/a | yes | Positional argument `body`. |
| in | `T` | Any | n/a | yes | Positional argument `T`. |
| in | `S` | Any | n/a | yes | Positional argument `S`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `montecarlo` | Any | n/a | no | Positional argument `montecarlo` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `aerodynamic_coefficient_constant`. Returns `CL_body, CD_body`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.trajectory_predictor_f_ctrl_rf_bang|f_ctrl_rf!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:316-316`
- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:316-316`
- [[grp.src_gnc_guidance|gnc/guidance/]] · `members_out` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:316-316`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Uses `args.α` instead of the `α` positional argument, so the passed incidence is ignored and the function throws unless `args` has an `α` field. `body`, `T`, and `S` are unused. Not called by the active `wrench` path, which uses `_constant_drag_coefficient` instead. The docstring is an empty string placed inside the function body.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 1005.
