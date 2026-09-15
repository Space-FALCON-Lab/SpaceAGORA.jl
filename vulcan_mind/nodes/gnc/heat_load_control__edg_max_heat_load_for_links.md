---
id: gnc.heat_load_control__edg_max_heat_load_for_links
label: _edg_max_heat_load_for_links
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_max_heat_load_for_links
  lines:
  - 35
  - 35
inputs:
- id: sc
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc`.
- id: links
  type: Tuple{Vararg{Int}}
  units: n/a
  required: true
  description: Positional argument `links`.
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
  description: Return value of `_edg_max_heat_load_for_links`.
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

# _edg_max_heat_load_for_links

## Purpose
Returns the largest accumulated heat load among a selected tuple of link indices, used to compare against the configured heat-load limit for the controlled panels.

## Design & Implementation
Takes `sc` (a spacecraft state record) and `links::Tuple{Vararg{Int}}`. Returns `0.0` immediately if `sc` lacks a `heat_loads` property. Otherwise it iterates the indices, skips any outside `1:length(heat_loads)`, converts each to `Float64`, ignores non-finite entries, and accumulates `max(value, max(0.0, candidate))`. Units follow `heat_loads`, which the caller treats as J/cm^2.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `sc` | Any | n/a | yes | Positional argument `sc`. |
| in | `links` | Tuple{Vararg{Int}} | n/a | yes | Positional argument `links`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_edg_max_heat_load_for_links`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:263-263`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_load_control.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/heat_load_control.jl:41-41`
<!-- vulcan:connections:end -->

## Limitations
Out-of-range and NaN indices are silently ignored, so a typo in `controlled_panel_links` produces a zero heat load rather than an error. Negative loads are clamped to zero. The tuple type forces compile-time specialisation per link-count.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 35.
