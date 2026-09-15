---
id: gnc.propulsive_maneuvers__maneuver_trace_enabled
label: _maneuver_trace_enabled
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _maneuver_trace_enabled
  lines:
  - 30
  - 30
inputs:
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
  type: Bool
  units: n/a
  description: Return value of `_maneuver_trace_enabled`.
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

# _maneuver_trace_enabled

## Purpose
Master switch for the maneuver-event CSV trace, enabled either explicitly or implicitly by naming an output file.

## Design & Implementation
Returns the logical OR of two conditions: `_trace_bool_enabled(get(ENV, "SPACEAGORA_TRACE_MANEUVERS", "0"))`, and `!isempty(strip(get(ENV, "SPACEAGORA_MANEUVER_TRACE_CSV", "")))`. The second arm means that setting a destination path is itself sufficient to turn tracing on. `_trace_maneuver_event!` calls this first and returns `nothing` immediately when it is false.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_maneuver_trace_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers__trace_maneuver_event_bang|_trace_maneuver_event!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:290-290`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`

**Downstream**

- `callees` → [[gnc.propulsive_maneuvers__trace_bool_enabled|_trace_bool_enabled]] · `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:31-31`
<!-- vulcan:connections:end -->

## Limitations
Because a non-empty CSV path forces tracing on, a user who sets the path once and later wants tracing off must unset the variable entirely; setting `SPACEAGORA_TRACE_MANEUVERS=0` will not disable it. Both variables are read on every traced event rather than captured once.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 30.
