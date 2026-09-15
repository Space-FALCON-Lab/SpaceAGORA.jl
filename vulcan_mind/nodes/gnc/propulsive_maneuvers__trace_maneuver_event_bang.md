---
id: gnc.propulsive_maneuvers__trace_maneuver_event_bang
label: _trace_maneuver_event!
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _trace_maneuver_event!
  lines:
  - 277
  - 277
inputs:
- id: event
  type: String
  units: n/a
  required: true
  description: Positional argument `event`.
- id: controlModel
  type: BaseThrusterModel
  units: n/a
  required: true
  description: Positional argument `controlModel`.
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
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: start_burn_s
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `start_burn_s` (default `NaN`).
- id: stop_burn_s
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `stop_burn_s` (default `NaN`).
- id: alt_m
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `alt_m` (default `NaN`).
- id: e
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `e` (default `NaN`).
- id: nu_rad
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `nu_rad` (default `NaN`).
- id: a_m
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `a_m` (default `NaN`).
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
  description: Return value of `_trace_maneuver_event!`; mutates `event` in place.
    Returns `nothing`.
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

# _trace_maneuver_event!

## Purpose
Appends one row to the maneuver-event CSV trace, recording a burn start, burn end, schedule set, schedule update or schedule clear together with the state that produced it.

## Design & Implementation
Returns immediately unless `_maneuver_trace_enabled()`. It resolves the destination with `_maneuver_trace_path()`, reads the orbit counter through `_safe_orbit_counter`, and index-guards `controlModel.Δv`, `controlModel.direction` and `controlModel.thrust` individually, substituting `NaN` where `i` exceeds a length. Burn duration is derived as `stop_burn_s - start_burn_s` when both are finite, else `NaN`. All file work happens inside `lock(_MANEUVER_TRACE_LOCK)`: it calls `mkpath(dirname(path))`, records whether the file already existed, opens it in append mode, writes the fourteen-column header only for a new file, and then writes the row. Optional keywords `alt_m`, `e`, `nu_rad` and `a_m` all default to `NaN` and carry the orbital context for schedule events.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `event` | String | n/a | yes | Positional argument `event`. |
| in | `controlModel` | BaseThrusterModel | n/a | yes | Positional argument `controlModel`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `start_burn_s` | Float64 | n/a | no | Keyword argument `start_burn_s` (default `NaN`). |
| in | `stop_burn_s` | Float64 | n/a | no | Keyword argument `stop_burn_s` (default `NaN`). |
| in | `alt_m` | Float64 | n/a | no | Keyword argument `alt_m` (default `NaN`). |
| in | `e` | Float64 | n/a | no | Keyword argument `e` (default `NaN`). |
| in | `nu_rad` | Float64 | n/a | no | Keyword argument `nu_rad` (default `NaN`). |
| in | `a_m` | Float64 | n/a | no | Keyword argument `a_m` (default `NaN`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_trace_maneuver_event!`; mutates `event` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:475-475`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:293-293`
- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:303-303`
- `callees` → [[gnc.propulsive_maneuvers__maneuver_trace_enabled|_maneuver_trace_enabled]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:290-290`
- `callees` → [[gnc.propulsive_maneuvers__maneuver_trace_path|_maneuver_trace_path]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:291-291`
- `callees` → [[gnc.propulsive_maneuvers__safe_orbit_counter|_safe_orbit_counter]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:292-292`
<!-- vulcan:connections:end -->

## Limitations
The file is opened, written and closed under the global lock on every event, so tracing serialises all threads in the control callback and can dominate step time in a large constellation. Values are written with `string(...)`, giving full Float64 repr including literal `NaN`, which many CSV readers will not coerce to a missing value. Append mode with a header written only for a new file means reruns silently concatenate onto stale data.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 277.
