---
id: vehicle.model_controlmodel
label: ControlModel
kind: struct
source:
  file: src/vehicle/spacecraft/model.jl
  symbol: ControlModel
  lines:
  - 460
  - 460
inputs:
- id: control_effectors
  type: T_Effectors
  units: n/a
  required: true
  description: Field `control_effectors`.
- id: control_rates
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `control_rates`.
- id: n_effectors
  type: Any
  units: n/a
  required: false
  description: Field `n_effectors` (default `length(control_effectors)`).
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
  type: ControlModel
  units: n/a
  description: Constructed `ControlModel` (keyword constructor via @kwdef).
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

# ControlModel

## Purpose
Bundles the tuple of control effectors (reaction wheels, thrusters, solar-panel controllers) together with the cadence at which each is invoked, validated at construction so the simulation loop can trust the rate vector.

## Design & Implementation
`@kwdef struct ControlModel{T_Effectors<:Tuple}` with fields `control_effectors::T_Effectors` and `control_rates::Vector{Float64}` (seconds between calls). The inner constructor checks `length(control_rates) == length(control_effectors)` and that every rate `isfinite` and `> 0.0`, throwing `ArgumentError` with the index and offending value otherwise. The tuple type parameter keeps effector dispatch static.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `control_effectors` | T_Effectors | n/a | yes | Field `control_effectors`. |
| in | `control_rates` | Vector{Float64} | n/a | yes | Field `control_rates`. |
| in | `n_effectors` | Any | n/a | no | Field `n_effectors` (default `length(control_effectors)`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ControlModel | n/a | — | Constructed `ControlModel` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__with_campaign_maneuvers|_with_campaign_maneuvers]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:444-444`
- [[envana.ana_example_support_make_example_config|make_example_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:163-163`
- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/model.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Rates are validated for positivity but not against the integration step; a rate smaller than the solver step simply results in calls every step. Effectors are heterogeneous by design, so a large tuple increases compile time and specialisation. `control_rates` is a mutable `Vector`, so it can be altered after validation.

## Provenance
Mapped from `src/vehicle/spacecraft/model.jl` line 460.
