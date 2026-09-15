---
id: vehicle.thruster_hooks_thrust_calculation_schmitt_trigger_bang
label: thrust_calculation_schmitt_trigger!
kind: function
source:
  file: src/vehicle/actuators/thruster/thruster_hooks.jl
  symbol: thrust_calculation_schmitt_trigger!
  lines:
  - 55
  - 55
inputs:
- id: link
  type: Link
  units: n/a
  required: true
  description: Positional argument `link`.
- id: thruster
  type: Thruster
  units: n/a
  required: true
  description: Positional argument `thruster`.
- id: thrust
  type: Float64
  units: n/a
  required: true
  description: Positional argument `thrust`.
- id: time
  type: Float64
  units: n/a
  required: true
  description: Positional argument `time`.
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
  type: Nothing
  units: n/a
  description: Return value of `thrust_calculation_schmitt_trigger!`; mutates `link`
    in place. Returns `nothing`.
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

# thrust_calculation_schmitt_trigger!

## Purpose

Converts one thruster's continuous thrust demand into a pulse-width command, integrates the resulting impulse and writes the period-averaged thrust back into `thruster.thrust`. It is called once per thruster per attitude control step by `update_thrusters!`.

## Design & Implementation

If `thruster.max_thrust` is non-finite or non-positive the thrust is zeroed and the function returns immediately. Otherwise the on-time is `ti = min(thrust / max_thrust * link.attitude_control_rate, link.attitude_control_rate)`, capping duty cycle at 100%. When `ti` is below `thruster.min_firing_time` it is replaced by `schmitt_trigger(ti, thruster.level_on, thruster.level_off) * thruster.min_firing_time`, so short demands either snap up to the minimum pulse or drop to zero. With debug logging enabled the tuple `(time, on_time_request, thrust_req)` is appended to `thruster_debug.csv`. Finally `integrate_impulse!` returns the total impulse, which is divided by `link.attitude_control_rate` to give the average thrust.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `link` | Link | n/a | yes | Positional argument `link`. |
| in | `thruster` | Thruster | n/a | yes | Positional argument `thruster`. |
| in | `thrust` | Float64 | n/a | yes | Positional argument `thrust`. |
| in | `time` | Float64 | n/a | yes | Positional argument `time`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `thrust_calculation_schmitt_trigger!`; mutates `link` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[vehx.actuators_thruster_hooks_update_thrusters_bang|update_thrusters!]] · `callees` → `callers` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl:51-51`

**Downstream**

- `callees` → [[vehicle.thruster_hooks_integrate_impulse_bang|integrate_impulse!]] · `callers` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl:71-71`
- `callees` → [[vehicle.thruster_hooks_schmitt_trigger|schmitt_trigger]] · `callers` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl:65-65`
- `callees` → [[vehicle.thruster_hooks_thruster_debug_enabled|thruster_debug_enabled]] · `callers` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl:68-68`
<!-- vulcan:connections:end -->

## Limitations

Pulse-width modulation is single-pulse per control period and the quantisation to `min_firing_time` introduces a torque error the caller does not compensate. The debug write re-opens and appends to `thruster_debug.csv` in the working directory on every call, which is slow and not safe across concurrently simulated vehicles. The averaged thrust reported back is a period mean, so instantaneous force during the pulse is not visible to downstream consumers.

## Provenance
Mapped from `src/vehicle/actuators/thruster/thruster_hooks.jl` line 55.
