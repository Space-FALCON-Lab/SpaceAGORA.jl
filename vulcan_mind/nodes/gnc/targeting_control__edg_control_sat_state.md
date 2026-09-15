---
id: gnc.targeting_control__edg_control_sat_state
label: _edg_control_sat_state
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_control_sat_state
  lines:
  - 53
  - 53
inputs:
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
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
  type: Any
  units: n/a
  description: 'Return value of `_edg_control_sat_state`. Returns `hasproperty(u,
    :sc) ? u.sc[i] : u`.'
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

# _edg_control_sat_state

## Purpose
Selects satellite `i`'s state from the integrator's component-tree state, or passes a single-satellite state through.

## Design & Implementation
Returns `u.sc[i]` when `u` has an `sc` property and `u` itself otherwise. `@inline`. This lets the same control code run under both the multi-satellite tree layout and a flat single-vehicle test harness.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `i` | Int | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_control_sat_state`. Returns `hasproperty(u, :sc) ? u.sc[i] : u`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_environment_state|_edg_environment_state]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:65-65`
- [[gnc.targeting_control_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:259-259`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
For a flat state the index `i` is ignored, so calling with `i > 1` silently returns the same vehicle.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 53.
