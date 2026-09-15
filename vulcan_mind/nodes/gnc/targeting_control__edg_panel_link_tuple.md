---
id: gnc.targeting_control__edg_panel_link_tuple
label: _edg_panel_link_tuple
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_panel_link_tuple
  lines:
  - 14
  - 14
inputs:
- id: controlled_panel_links
  type: Any
  units: n/a
  required: true
  description: Positional argument `controlled_panel_links`.
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
  description: Return value of `_edg_panel_link_tuple`. Returns `links`.
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

# _edg_panel_link_tuple

## Purpose
Normalises and validates the controlled-link specification for the panel effector at construction time.

## Design & Implementation
Converts each entry to `Int` into a tuple, raises `ArgumentError` if the tuple is empty or any index is non-positive. `@inline`. The upper bound is not known until a spacecraft is attached, so it is checked later by `_apply_solar_panel_aoa!`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `controlled_panel_links` | Any | n/a | yes | Positional argument `controlled_panel_links`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_panel_link_tuple`. Returns `links`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[gnc.targeting_control_solarpanelangleofattackcontrolmodel|SolarPanelAngleOfAttackControlModel]] · `callers` · call · `src/gnc/control/targeting_control.jl:21-21`
<!-- vulcan:connections:end -->

## Limitations
Accepts duplicate indices, which would rotate the same link twice per tick, doubling its deflection.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 14.
