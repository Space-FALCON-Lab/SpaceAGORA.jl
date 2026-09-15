---
id: vehicle.model_guidancemodel
label: GuidanceModel
kind: struct
source:
  file: src/vehicle/spacecraft/model.jl
  symbol: GuidanceModel
  lines:
  - 426
  - 426
inputs:
- id: guidance_effectors
  type: T_Effectors
  units: n/a
  required: true
  description: Field `guidance_effectors`.
- id: guidance_rates
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `guidance_rates`.
- id: n_effectors
  type: Any
  units: n/a
  required: false
  description: Field `n_effectors` (default `length(guidance_effectors)`).
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
  type: GuidanceModel
  units: n/a
  description: Constructed `GuidanceModel` (keyword constructor via @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# GuidanceModel

## Purpose
Groups guidance effectors (maneuver planners, apoapsis targeting) with their call rates, enforcing at construction that the two vectors line up and that every rate is a finite positive number of seconds.

## Design & Implementation
`@kwdef struct GuidanceModel{T_Effectors<:Tuple}` holding `guidance_effectors::T_Effectors` and `guidance_rates::Vector{Float64}`. The inner constructor throws `ArgumentError("GuidanceModel guidance_rates length (...) must match guidance_effectors length (...)")` on mismatch and `ArgumentError("GuidanceModel guidance_rate at index ... must be finite and > 0.0, got ...")` for any bad rate, iterating with `@inbounds pairs(guidance_rates)`. Structure and validation are identical to `NavigationModel` and `ControlModel`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `guidance_effectors` | T_Effectors | n/a | yes | Field `guidance_effectors`. |
| in | `guidance_rates` | Vector{Float64} | n/a | yes | Field `guidance_rates`. |
| in | `n_effectors` | Any | n/a | no | Field `n_effectors` (default `length(guidance_effectors)`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | GuidanceModel | n/a | — | Constructed `GuidanceModel` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__with_campaign_maneuvers|_with_campaign_maneuvers]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:439-439`
- [[envana.ana_example_support_make_example_config|make_example_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:161-161`
- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/model.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The three GNC model structs duplicate the same validation code rather than sharing it. Rates are stored in a mutable vector so post-construction edits bypass validation. Nothing relates a guidance rate to the corresponding control rate even though guidance outputs are consumed by controllers.

## Provenance
Mapped from `src/vehicle/spacecraft/model.jl` line 426.
