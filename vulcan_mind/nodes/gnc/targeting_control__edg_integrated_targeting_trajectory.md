---
id: gnc.targeting_control__edg_integrated_targeting_trajectory
label: _edg_integrated_targeting_trajectory
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_integrated_targeting_trajectory
  lines:
  - 457
  - 457
inputs:
- id: config
  type: AerobrakingEnergyDepletionConfig
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
- id: pos0
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos0`.
- id: vel0
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel0`.
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
- id: times
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `times`.
- id: switch_time_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `switch_time_s`.
- id: heat_rate_control
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `heat_rate_control`.
- id: structural_control
  type: Bool
  units: n/a
  required: true
  description: Keyword argument `structural_control`.
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
  description: Return value of `_edg_integrated_targeting_trajectory`. Returns `gravity
    + aero` or `(`.
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

# _edg_integrated_targeting_trajectory

## Purpose
Integrates the vehicle through a drag passage under the targeting profile — maximum drag until a switch time, then minimum — returning the full predicted track with its angle profile.

## Design & Implementation
Preallocates position, velocity and angle arrays over `times`. At each step it picks the base angle by comparing `t + tau` with `switch_time_s`, constrains it through `_edg_targeting_constrained_alpha` warm-started from the previous step, and advances with a hand-written RK4 holding `alpha` fixed across the stages. After the loop it evaluates the final angle and then makes a second pass over every sample to record altitude, flight-path angle, speed, density, temperature and speed ratio. Returns a named tuple.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `config` | AerobrakingEnergyDepletionConfig | n/a | yes | Positional argument `config`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `pos0` | SVector{3, Float64} | n/a | yes | Positional argument `pos0`. |
| in | `vel0` | SVector{3, Float64} | n/a | yes | Positional argument `vel0`. |
| in | `mass` | Float64 | n/a | yes | Positional argument `mass`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `times` | Vector{Float64} | n/a | yes | Positional argument `times`. |
| in | `switch_time_s` | Float64 | n/a | yes | Positional argument `switch_time_s`. |
| in | `heat_rate_control` | Bool | n/a | yes | Keyword argument `heat_rate_control`. |
| in | `structural_control` | Bool | n/a | yes | Keyword argument `structural_control`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_integrated_targeting_trajectory`. Returns `gravity + aero` or `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_predict_targeting_outcome|_edg_predict_targeting_outcome]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:685-685`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Holding the angle fixed across the four RK4 stages is a first-order treatment of the switch; the second pass re-samples the environment at every point, doubling the density-model cost of the prediction.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 457.
