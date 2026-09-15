---
id: gnc.target_energy_bracketing__edg_state_index_ok
label: _edg_state_index_ok
kind: function
source:
  file: src/gnc/guidance/target_energy_bracketing.jl
  symbol: _edg_state_index_ok
  lines:
  - 176
  - 176
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
  description: Return value of `_edg_state_index_ok`.
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

# _edg_state_index_ok

## Purpose
`_edg_state_index_ok` is the bounds guard that `calcGuidanceEffect!` applies before touching any per-spacecraft vector in `AerobrakingEnergyDepletionState`, so a guidance call for a spacecraft index outside the allocated state is ignored rather than raising a `BoundsError` inside the ODE right-hand side.

## Design & Implementation
Declared `@inline` with signature `(state::AerobrakingEnergyDepletionState, i::Int)::Bool`. It returns `1 <= i <= length(state.selected_mode)`, using `selected_mode` as the canonical length because every vector in the state is constructed with the same `num_sats`. There are no side effects and no allocation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `state` | AerobrakingEnergyDepletionState | n/a | yes | Positional argument `state`. |
| in | `i` | Int | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_edg_state_index_ok`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.target_energy_bracketing_calcguidanceeffect_bang|calcGuidanceEffect!]] · `callees` → `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:199-199`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/target_energy_bracketing.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the length of `selected_mode` is checked; if any other state vector were resized independently the guard would not detect the mismatch. A `false` result causes `calcGuidanceEffect!` to return silently, so a mis-sized state produces a spacecraft with no guidance and no diagnostic. The check is repeated on every guidance evaluation, although its cost is negligible.

## Provenance
Mapped from `src/gnc/guidance/target_energy_bracketing.jl` line 176.
