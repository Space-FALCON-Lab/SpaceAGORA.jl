---
id: gnc.propulsive_maneuvers__maneuver_trace_key
label: _maneuver_trace_key
kind: function
source:
  file: src/gnc/control/propulsive_maneuvers.jl
  symbol: _maneuver_trace_key
  lines:
  - 40
  - 40
inputs:
- id: controlModel
  type: BaseThrusterModel
  units: n/a
  required: true
  description: Positional argument `controlModel`.
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
  type: Tuple{UInt64,
  units: n/a
  description: Return value of `_maneuver_trace_key`.
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

# _maneuver_trace_key

## Purpose
Builds the dictionary key that identifies one spacecraft's burn schedule within the process-global maneuver-trace state.

## Design & Implementation
Returns the tuple `(UInt64(objectid(controlModel.start_burn_time)), i)`. Using `objectid` of the `start_burn_time` array rather than of the model itself keys on the identity of the mutable schedule buffer, so two thruster models sharing one buffer collapse to the same key and a model whose buffer is replaced gets a fresh key. The tuple type `Tuple{UInt64, Int64}` matches `_MANEUVER_TRACE_LAST_WINDOW` and `_MANEUVER_TRACE_BURN_ACTIVE`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `controlModel` | BaseThrusterModel | n/a | yes | Positional argument `controlModel`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{UInt64, | n/a | — | Return value of `_maneuver_trace_key`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.propulsive_maneuvers_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:443-443`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`objectid` is a hash of object identity and is not guaranteed collision-free, and it is not stable across processes or sessions, so trace state cannot be checkpointed and restored. Because the arrays are garbage-collectable, entries in the two global dictionaries keyed by a freed array's id are never reclaimed except by the explicit `pop!` on schedule clear.

## Provenance
Mapped from `src/gnc/control/propulsive_maneuvers.jl` line 40.
