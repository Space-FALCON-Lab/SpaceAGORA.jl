---
id: vehicle.model_navigationmodel
label: NavigationModel
kind: struct
source:
  file: src/vehicle/spacecraft/model.jl
  symbol: NavigationModel
  lines:
  - 443
  - 443
inputs:
- id: navigation_effectors
  type: T_Effectors
  units: n/a
  required: true
  description: Field `navigation_effectors`.
- id: navigation_rates
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `navigation_rates`.
- id: n_effectors
  type: Any
  units: n/a
  required: false
  description: Field `n_effectors` (default `length(navigation_effectors)`).
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
  type: NavigationModel
  units: n/a
  description: Constructed `NavigationModel` (keyword constructor via @kwdef).
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

# NavigationModel

## Purpose
Pairs navigation effectors (sensor and estimator models) with the interval at which each is executed, validating the pairing at construction so the callback scheduler can rely on the rate vector.

## Design & Implementation
`@kwdef struct NavigationModel{T_Effectors<:Tuple}` with `navigation_effectors::T_Effectors` and `navigation_rates::Vector{Float64}` (s). The inner constructor throws `ArgumentError` when the lengths differ or when any rate is non-finite or `<= 0.0`, reporting the index and value. It mirrors `GuidanceModel` and `ControlModel` exactly apart from field names.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `navigation_effectors` | T_Effectors | n/a | yes | Field `navigation_effectors`. |
| in | `navigation_rates` | Vector{Float64} | n/a | yes | Field `navigation_rates`. |
| in | `n_effectors` | Any | n/a | no | Field `n_effectors` (default `length(navigation_effectors)`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | NavigationModel | n/a | — | Constructed `NavigationModel` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[envana.ana_example_support_make_example_config|make_example_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:162-162`
- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/model.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Validation code is triplicated across the three GNC models. The rates vector is mutable after construction. There is no ordering guarantee or dependency declaration between navigation outputs and the guidance/control effectors that consume them; ordering is fixed by callback registration elsewhere.

## Provenance
Mapped from `src/vehicle/spacecraft/model.jl` line 443.
