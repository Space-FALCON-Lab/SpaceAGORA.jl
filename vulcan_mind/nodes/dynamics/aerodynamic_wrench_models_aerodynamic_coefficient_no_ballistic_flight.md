---
id: dynamics.aerodynamic_wrench_models_aerodynamic_coefficient_no_ballistic_flight
label: aerodynamic_coefficient_no_ballistic_flight
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: aerodynamic_coefficient_no_ballistic_flight
  lines:
  - 1091
  - 1091
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
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
- id: T
  type: Any
  units: n/a
  required: false
  description: Positional argument `T` (default `0`).
- id: S
  type: Any
  units: n/a
  required: false
  description: Positional argument `S` (default `0`).
- id: a
  type: Any
  units: n/a
  required: false
  description: Positional argument `a` (default `0`).
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
  description: Return value of `aerodynamic_coefficient_no_ballistic_flight`. Returns
    `CL_body, CD_body`.
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

# aerodynamic_coefficient_no_ballistic_flight

## Purpose
Modified-Newtonian lift and drag coefficients for a sphere-cone blunt body, parameterised by nose-to-base radius ratio and cone half-angle.

## Theory & Math
Modified Newtonian sphere-cone: $C_A = (1-\sin^4\delta)k^2 + (2\sin^2\delta\cos^2\alpha + \cos^2\delta\sin^2\alpha)(1 - k^2\cos^2\delta)$, $C_N = (1 - k^2\cos^2\delta)\cos^2\delta\,\sin 2\alpha$, with $k = R_n/R_b$ and $\delta$ the cone half-angle.

## Design & Implementation
Signature `(α, body, args, T=0, S=0, a=0, montecarlo=false)`. Reads `k = body.nose_radius / body.base_radius` and half-angle `δ = body.δ`. Computes `CA = (1 - sin⁴δ) k² + (2 sin²δ cos²α + cos²δ sin²α)(1 - (k cosδ)²)`, `CN = (1 - (k cosδ)²) cos(δ²) sin(2α)`, then `CD = CA cosα + CN sinα - 0.15` and `CL = CN cosα - CA sinα`, with an optional `monte_carlo_aerodynamics` perturbation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `alpha` | Any | n/a | yes | Positional argument `α`. |
| in | `body` | Any | n/a | yes | Positional argument `body`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `T` | Any | n/a | no | Positional argument `T` (default `0`). |
| in | `S` | Any | n/a | no | Positional argument `S` (default `0`). |
| in | `a` | Any | n/a | no | Positional argument `a` (default `0`). |
| in | `montecarlo` | Any | n/a | no | Positional argument `montecarlo` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `aerodynamic_coefficient_no_ballistic_flight`. Returns `CL_body, CD_body`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.trajectory_predictor_f_ctrl_rf_bang|f_ctrl_rf!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:320-320`
- [[gncy.trajectory_predictor_asim_ctrl_targeting|asim_ctrl_targeting]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:320-320`
- [[grp.src_gnc_guidance|gnc/guidance/]] · `members_out` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/trajectory_predictor.jl:320-320`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The code writes `cos(δ^2)` where the standard form is `cos(δ)^2`, so `CN` is wrong for any non-trivial half-angle. The `- 0.15` base-drag correction and `Cp_max = 2` (assigned but unused) are hard-coded. Not wired into any `wrench` path, so it is currently unreachable in simulation. Unused parameters `T`, `S`, `a`.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 1091.
