---
id: gnc.heat_load_control__edg_heat_load_low_drag_active
label: _edg_heat_load_low_drag_active
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_heat_load_low_drag_active
  lines:
  - 754
  - 754
inputs:
- id: model
  type: Any
  units: n/a
  required: true
  description: Positional argument `model`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  description: Return value of `_edg_heat_load_low_drag_active`.
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

# _edg_heat_load_low_drag_active

## Purpose
Fast predicate used in the guidance loop to test whether the current time falls inside a satellite's solved low-drag heat-load window.

## Design & Implementation
`@inline` function `(model, t::Float64, i::Int)::Bool`. Reads `switches = model.state.heat_load_switches_s[i]` and returns true only when `model.state.selected_mode[i] == :max_energy_depletion`, `:heat_load in model.config.max_energy_submodes`, both switch times are finite, and `switches[1] <= t <= switches[2]`. The `(Inf, Inf)` sentinel produced by the solver therefore always yields `false`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | Any | n/a | yes | Positional argument `model`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `i` | Int | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_edg_heat_load_low_drag_active`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_base_alpha|_edg_base_alpha]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:166-166`
- [[gnc.targeting_control_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:265-265`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_load_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No bounds check on `i`; an invalid satellite index throws `BoundsError`. The membership test `:heat_load in config.max_energy_submodes` is a linear scan on every RHS call, though the collection is tiny. The predicate does not verify `switches[1] <= switches[2]`.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 754.
