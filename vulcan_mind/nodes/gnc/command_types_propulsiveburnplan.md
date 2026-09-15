---
id: gnc.command_types_propulsiveburnplan
label: PropulsiveBurnPlan
kind: struct
source:
  file: src/gnc/command_types.jl
  symbol: PropulsiveBurnPlan
  lines:
  - 12
  - 12
inputs:
- id: valid
  type: Bool
  units: n/a
  required: false
  description: Field `valid` (default `false`).
- id: delta_v_mps
  type: Float64
  units: n/a
  required: false
  description: Field `delta_v_mps` (default `0.0`).
- id: direction_rad
  type: Float64
  units: n/a
  required: false
  description: Field `direction_rad` (default `0.0`).
- id: source_orbit
  type: Int64
  units: n/a
  required: false
  description: Field `source_orbit` (default `0`).
- id: start_burn_s
  type: Float64
  units: n/a
  required: false
  description: Field `start_burn_s` (default `-1.0`).
- id: stop_burn_s
  type: Float64
  units: n/a
  required: false
  description: Field `stop_burn_s` (default `-1.0`).
- id: thrust_n
  type: Float64
  units: n/a
  required: false
  description: Field `thrust_n` (default `0.0`).
- id: isp_s
  type: Float64
  units: n/a
  required: false
  description: Field `isp_s` (default `0.0`).
- id: commanded_impulse_n_s
  type: Float64
  units: n/a
  required: false
  description: Field `commanded_impulse_n_s` (default `0.0`).
- id: propellant_required_kg
  type: Float64
  units: n/a
  required: false
  description: Field `propellant_required_kg` (default `0.0`).
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
  type: PropulsiveBurnPlan
  units: n/a
  description: Constructed `PropulsiveBurnPlan` (keyword constructor via @kwdef).
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

# PropulsiveBurnPlan

## Purpose
`PropulsiveBurnPlan` is the fully-scheduled form of a propulsive manoeuvre: it carries not only the commanded delta-v and its in-plane direction, but the absolute simulation times at which the thruster must open and close, plus the engine parameters and propellant budget needed to realise the burn. Controllers consume it to gate thrust on and off during integration.

## Design & Implementation
Declared with `Base.@kwdef` as an immutable struct with ten concretely-typed fields: `valid::Bool`, `delta_v_mps::Float64`, `direction_rad::Float64`, `source_orbit::Int64`, `start_burn_s`/`stop_burn_s::Float64`, `thrust_n::Float64`, `isp_s::Float64`, `commanded_impulse_n_s::Float64` and `propellant_required_kg::Float64`. Times are seconds of simulation time, angles radians, thrust newtons, impulse newton-seconds. `start_burn_s` and `stop_burn_s` default to `-1.0` so an unplanned burn can never fall inside a non-negative integration window.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `valid` | Bool | n/a | no | Field `valid` (default `false`). |
| in | `delta_v_mps` | Float64 | n/a | no | Field `delta_v_mps` (default `0.0`). |
| in | `direction_rad` | Float64 | n/a | no | Field `direction_rad` (default `0.0`). |
| in | `source_orbit` | Int64 | n/a | no | Field `source_orbit` (default `0`). |
| in | `start_burn_s` | Float64 | n/a | no | Field `start_burn_s` (default `-1.0`). |
| in | `stop_burn_s` | Float64 | n/a | no | Field `stop_burn_s` (default `-1.0`). |
| in | `thrust_n` | Float64 | n/a | no | Field `thrust_n` (default `0.0`). |
| in | `isp_s` | Float64 | n/a | no | Field `isp_s` (default `0.0`). |
| in | `commanded_impulse_n_s` | Float64 | n/a | no | Field `commanded_impulse_n_s` (default `0.0`). |
| in | `propellant_required_kg` | Float64 | n/a | no | Field `propellant_required_kg` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | PropulsiveBurnPlan | n/a | — | Constructed `PropulsiveBurnPlan` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.runtime_types_sharedbuffers|SharedBuffers]] · `callees` → `callers` · call · `src/core/types/runtime_types.jl:736-736`
- [[gnc.propulsive_maneuvers__clear_burn_plan_bang|_clear_burn_plan!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:96-96`
- [[gnc.propulsive_maneuvers__validated_burn_plan|_validated_burn_plan]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:237-237`
- [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:546-546`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The struct is a record only: it never validates that `commanded_impulse_n_s` equals `thrust_n * (stop_burn_s - start_burn_s)`, nor that `propellant_required_kg` follows from `isp_s` by the rocket equation, so an inconsistent plan constructed by a caller propagates unchecked into the effector. `direction_rad` is not wrapped to any interval, and `source_orbit` carries no meaning if the trajectory has no orbit counter.

## Provenance
Mapped from `src/gnc/command_types.jl` line 12.
