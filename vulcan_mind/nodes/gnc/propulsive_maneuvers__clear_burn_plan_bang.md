---
id: gnc.propulsive_maneuvers__clear_burn_plan_bang
label: _clear_burn_plan!
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _clear_burn_plan!
  lines:
  - 91
  - 91
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
  description: Return value of `_clear_burn_plan!`; mutates `p` in place. Returns
    `nothing`.
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

# _clear_burn_plan!

## Purpose
Retires the burn plan for spacecraft `i` once its burn has completed, so a later campaign maneuver can be planned.

## Design & Implementation
Mirrors `_set_burn_plan!`: it resolves the buffer, bounds-checks `i`, and assigns `plans[i] = PropulsiveBurnPlan()` — a default-constructed record whose `valid` field is false, which makes `_active_burn_plan` return `nothing` from then on. It returns `nothing` unconditionally. `calcControlEffect!` invokes it on the `schedule_cleared` branch, alongside resetting `start_burn_time[i]` and `stop_burn_time[i]` to `-1.0`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_clear_burn_plan!`; mutates `p` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:487-487`

**Downstream**

- `callees` → [[gnc.command_types_propulsiveburnplan|PropulsiveBurnPlan]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:96-96`
- `callees` → [[gnc.propulsive_maneuvers__burn_plan_buffer|_burn_plan_buffer]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:92-92`
<!-- vulcan:connections:end -->

## Limitations
Clearing replaces the record rather than archiving it, so the executed plan's commanded impulse and propellant figures are lost unless they were traced. Like its sibling it fails silently when the buffer is missing, and it performs an unlocked write to shared state.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 91.
