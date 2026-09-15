---
id: gnc.propulsive_maneuvers__safe_orbit_counter
label: _safe_orbit_counter
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _safe_orbit_counter
  lines:
  - 44
  - 44
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
  type: Int64
  units: n/a
  description: Return value of `_safe_orbit_counter`.
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

# _safe_orbit_counter

## Purpose
Reads the current orbit-revolution counter for one spacecraft for trace annotation, without letting a missing or malformed counter abort the traced event.

## Design & Implementation
Wraps `Int64(p.orbit_counter[i])` in a `try` block and returns `Int64(-1)` from the `catch`. The sentinel `-1` therefore appears in the `orbit_counter` column of the CSV whenever the parameter object has no counter, the index is out of range, or the stored value is not integer-convertible.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int64 | n/a | — | Return value of `_safe_orbit_counter`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers__trace_maneuver_event_bang|_trace_maneuver_event!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:292-292`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/propulsive_maneuvers.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The bare `catch` cannot distinguish an absent field from an out-of-range index or a `NaN` counter, so all three appear identically as `-1` in the trace. A legitimate stored value of `-1` is indistinguishable from the failure sentinel.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 44.
