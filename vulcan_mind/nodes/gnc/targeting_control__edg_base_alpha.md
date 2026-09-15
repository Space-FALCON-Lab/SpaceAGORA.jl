---
id: gnc.targeting_control__edg_base_alpha
label: _edg_base_alpha
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_base_alpha
  lines:
  - 157
  - 157
inputs:
- id: model
  type: AerobrakingEnergyDepletionControlModel
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
  type: Float64
  units: n/a
  description: Return value of `_edg_base_alpha`.
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

# _edg_base_alpha

## Purpose
Determines the unconstrained angle-of-attack command from the guidance mode and the solved switch times, before any thermal or structural limits are applied.

## Design & Implementation
Returns `min_alpha_rad` for safe-low-drag mode or when the safe flag is set. In targeting mode it returns the maximum before `targeting_switch_s` and the minimum after. In max-energy-depletion it returns the minimum while the heat-load low-drag window is active. Every other case returns `max_alpha_rad`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | AerobrakingEnergyDepletionControlModel | n/a | yes | Positional argument `model`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `i` | Int | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_edg_base_alpha`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:266-266`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[gnc.heat_load_control__edg_heat_load_low_drag_active|_edg_heat_load_low_drag_active]] · `callers` · call · `src/gnc/control/targeting_control.jl:166-166`
<!-- vulcan:connections:end -->

## Limitations
The bang-bang switch produces a step in commanded angle; there is no slew-rate limiting here, so the panel effector is asked to jump instantaneously.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 157.
