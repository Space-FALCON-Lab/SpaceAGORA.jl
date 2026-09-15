---
id: vehicle.thruster_hooks_schmitt_trigger
label: schmitt_trigger
kind: function
source:
  file: src/vehicle/actuators/thruster/thruster_hooks.jl
  symbol: schmitt_trigger
  lines:
  - 76
  - 76
inputs:
- id: input
  type: Float64
  units: n/a
  required: true
  description: Positional argument `input`.
- id: level_on
  type: Float64
  units: n/a
  required: true
  description: Positional argument `level_on`.
- id: level_off
  type: Float64
  units: n/a
  required: true
  description: Positional argument `level_off`.
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
  description: Return value of `schmitt_trigger`. Returns `state`.
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

# schmitt_trigger

## Purpose

Implements the hysteresis decision that says whether a requested on-time is large enough to command a minimum-width thruster pulse. Given `input`, an upper threshold `level_on` and a lower threshold `level_off`, it returns `1.0` when the input exceeds `level_on` and `0.0` when it falls below `level_off`.

## Design & Implementation

The function initialises `state = 0.0`, sets it to `1.0` if `input > level_on`, and leaves it at `0.0` if `input < level_off`. The returned `Float64` is used multiplicatively by `thrust_calculation_schmitt_trigger!` to scale `thruster.min_firing_time`, so a `0.0` result suppresses the pulse entirely and a `1.0` result stretches it to the minimum firing width.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `input` | Float64 | n/a | yes | Positional argument `input`. |
| in | `level_on` | Float64 | n/a | yes | Positional argument `level_on`. |
| in | `level_off` | Float64 | n/a | yes | Positional argument `level_off`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `schmitt_trigger`. Returns `state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl`
- [[vehicle.thruster_hooks_thrust_calculation_schmitt_trigger_bang|thrust_calculation_schmitt_trigger!]] · `callees` → `callers` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl:65-65`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

Despite the name, the implementation is memoryless: no previous state is carried across calls, so the deadband between `level_off` and `level_on` resolves to `0.0` rather than holding the prior output. That makes it a two-threshold comparator, not a true latching Schmitt trigger, and it will chatter differently from a stateful one near the thresholds. All three arguments must be `Float64`.

## Provenance
Mapped from `src/vehicle/actuators/thruster/thruster_hooks.jl` line 76.
