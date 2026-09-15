---
id: gnc.targeting_control_solarpanelangleofattackcontrolmodel
label: SolarPanelAngleOfAttackControlModel
kind: struct
source:
  file: src/gnc/control/targeting_control.jl
  symbol: SolarPanelAngleOfAttackControlModel
  lines:
  - 10
  - 10
inputs:
- id: controlled_panel_links
  type: Tuple{Vararg{Int}}
  units: n/a
  required: true
  description: Field `controlled_panel_links`.
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
  type: SolarPanelAngleOfAttackControlModel
  units: n/a
  description: Constructed `SolarPanelAngleOfAttackControlModel`.
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

# SolarPanelAngleOfAttackControlModel

## Purpose
The control effector that physically articulates the solar-panel links to realise a commanded angle of attack during an aerobraking pass.

## Design & Implementation
An immutable subtype of `AbstractControlEffectorModel` holding `controlled_panel_links` as a tuple of positive integers. The keyword constructor defaults to links 2 and 3 — the two panels of the standard three-body vehicle — and validates through `_edg_panel_link_tuple`. Its `calcControlForceTorque` returns zero because articulation changes geometry, not applied force; the aerodynamic effectors then see the new link angle.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `controlled_panel_links` | Tuple{Vararg{Int}} | n/a | yes | Field `controlled_panel_links`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SolarPanelAngleOfAttackControlModel | n/a | — | Constructed `SolarPanelAngleOfAttackControlModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_panel_link_tuple|_edg_panel_link_tuple]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:21-21`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Rotation is applied about the absolute value of each link's mounting vector, which is only correct for panels mounted along a principal axis; an obliquely mounted panel rotates about the wrong axis.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 10.
