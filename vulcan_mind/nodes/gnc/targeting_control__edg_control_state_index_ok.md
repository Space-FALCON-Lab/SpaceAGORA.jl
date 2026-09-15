---
id: gnc.targeting_control__edg_control_state_index_ok
label: _edg_control_state_index_ok
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_control_state_index_ok
  lines:
  - 49
  - 49
inputs:
- id: state
  type: AerobrakingEnergyDepletionState
  units: n/a
  required: true
  description: Positional argument `state`.
- id: i
  type: Int
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
  type: Bool
  units: n/a
  description: Return value of `_edg_control_state_index_ok`.
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

# _edg_control_state_index_ok

## Purpose
Guards every per-satellite access to the energy-depletion state against an index outside the allocated vectors.

## Design & Implementation
Returns whether `1 <= i <= length(state.selected_mode)`, using the mode vector as the canonical length. `@inline` with a `::Bool` return. `calcControlEffect!` returns early on false.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `state` | AerobrakingEnergyDepletionState | n/a | yes | Positional argument `state`. |
| in | `i` | Int | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_edg_control_state_index_ok`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:254-254`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the mode vector's length is checked; the other per-satellite vectors in the state are assumed to have been allocated to the same length.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 49.
