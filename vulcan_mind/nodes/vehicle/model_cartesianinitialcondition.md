---
id: vehicle.model_cartesianinitialcondition
label: CartesianInitialCondition
kind: struct
source:
  file: src/vehicle/spacecraft/model.jl
  symbol: CartesianInitialCondition
  lines:
  - 226
  - 226
inputs:
- id: pos
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `pos`.
- id: vel
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `vel`.
- id: q
  type: SVector{4, Float64}
  units: n/a
  required: true
  description: Field `q`.
- id: ang_vel
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `ang_vel`.
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
  type: CartesianInitialCondition
  units: n/a
  description: Constructed `CartesianInitialCondition`.
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

# CartesianInitialCondition

## Purpose
Initial condition expressed directly as an inertial position and velocity, for scenarios (telemetry replays, RPO setups) where orbital elements are less natural than a state vector.

## Design & Implementation
`struct CartesianInitialCondition <: AbstractInitialCondition` with fields `pos::SVector{3,Float64}` (ECI, m), `vel::SVector{3,Float64}` (ECI, m/s), `q::SVector{4,Float64}` (orientation quaternion in `(x, y, z, w)` order) and `ang_vel::SVector{3,Float64}` (rad/s). A convenience constructor `CartesianInitialCondition(pos, vel; q=DEFAULT_INITIAL_CONDITION_Q, ang_vel=DEFAULT_INITIAL_CONDITION_ANG_VEL)` converts any 3-element iterables to `SVector{3,Float64}`; the defaults are the identity quaternion `(0,0,0,1)` and zero rates.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos` | SVector{3, Float64} | n/a | yes | Field `pos`. |
| in | `vel` | SVector{3, Float64} | n/a | yes | Field `vel`. |
| in | `q` | SVector{4, Float64} | n/a | yes | Field `q`. |
| in | `ang_vel` | SVector{3, Float64} | n/a | yes | Field `ang_vel`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | CartesianInitialCondition | n/a | — | Constructed `CartesianInitialCondition`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__initial_condition_from_time_aligned_telemetry|_initial_condition_from_time_aligned_telemetry]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:108-108`
- [[analysis.scenario_builders__initial_condition_in_j2000|_initial_condition_in_j2000]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:515-515`
- [[analysis.scenario_builders__scenario_initial_condition|_scenario_initial_condition]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:532-532`
- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/model.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No validation of the inputs: a zero position, a quaternion that is not unit norm, or non-finite values are accepted silently. The frame is documented as ECI only in comments; nothing records the epoch, so the caller must supply a consistent `initial_time` elsewhere.

## Provenance
Mapped from `src/vehicle/spacecraft/model.jl` line 226.
