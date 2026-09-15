---
id: gnc.propulsive_maneuvers__active_burn_plan
label: _active_burn_plan
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _active_burn_plan
  lines:
  - 73
  - 73
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: i
  type: Int64
  units: n/a
  required: true
  description: Positional argument `i`.
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
  description: 'Return value of `_active_burn_plan`. Returns `nothing` or `plan.valid
    ? plan : nothing`.'
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

# _active_burn_plan

## Purpose
Returns the currently valid burn plan for spacecraft `i`, or `nothing` if there is no buffer, the index is out of range, or the stored plan is not marked valid.

## Design & Implementation
Calls `_burn_plan_buffer(p)` and returns `nothing` when it yields `nothing` or when `i < 1 || i > length(plans)`. Otherwise it reads `plan = plans[i]` and returns it only when `plan.valid`. It is the override source consulted by `_effective_burn_window`, `_effective_direction_rad` and `_effective_thrust_isp`, letting a planned burn supply thrust, specific impulse, direction and window even when the thruster model's own arrays are unset.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_active_burn_plan`. Returns `nothing` or `plan.valid ? plan : nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers__effective_burn_window|_effective_burn_window]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:154-154`
- [[gnc.propulsive_maneuvers__effective_direction_rad|_effective_direction_rad]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:162-162`
- [[gnc.propulsive_maneuvers__effective_thrust_isp|_effective_thrust_isp]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:173-173`
- [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:445-445`

**Downstream**

- `callees` → [[gnc.propulsive_maneuvers__burn_plan_buffer|_burn_plan_buffer]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:74-74`
<!-- vulcan:connections:end -->

## Limitations
Validity is a flag on the record with no expiry, so a plan whose window has already passed is still returned as active until `_clear_burn_plan!` overwrites it. The three collapsed failure modes — no buffer, bad index, invalid plan — are indistinguishable to callers, all surfacing as `nothing`.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 73.
