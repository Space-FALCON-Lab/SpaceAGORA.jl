---
id: gnc.targeting_control__apply_solar_panel_aoa_bang
label: _apply_solar_panel_aoa!
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _apply_solar_panel_aoa!
  lines:
  - 229
  - 229
inputs:
- id: effector
  type: SolarPanelAngleOfAttackControlModel
  units: n/a
  required: true
  description: Positional argument `effector`.
- id: spacecraft
  type: Any
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
- id: alpha
  type: Float64
  units: n/a
  required: true
  description: Positional argument `alpha`.
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
  description: Return value of `_apply_solar_panel_aoa!`; mutates `effector` in place.
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

# _apply_solar_panel_aoa!

## Purpose
Rotates each controlled panel link to the commanded angle and records it on the link, so the aerodynamic model sees the new incidence.

## Design & Implementation
For each index in `controlled_panel_links`, validates that the index is in range and the link is not the root, raising `ArgumentError` otherwise. It forms the rotation axis as the absolute value of the link's mounting vector `r`, rotates by `π/2 - alpha` through `rotate_link`, and stores `alpha` in `link.α`. Mutates the spacecraft's links in place.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effector` | SolarPanelAngleOfAttackControlModel | n/a | yes | Positional argument `effector`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `alpha` | Float64 | n/a | yes | Positional argument `alpha`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_apply_solar_panel_aoa!`; mutates `effector` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:268-268`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[vehicle.kinematics_rotate_link|rotate_link]] · `callers` · call · `src/gnc/control/targeting_control.jl:240-240`
<!-- vulcan:connections:end -->

## Limitations
`rotate_link` sets an absolute orientation from the axis and angle each tick rather than applying an increment, which is what makes repeated calls idempotent, but it also means any other rotation applied to the panel is overwritten.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 229.
