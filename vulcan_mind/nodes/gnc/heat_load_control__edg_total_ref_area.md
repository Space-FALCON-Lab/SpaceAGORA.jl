---
id: gnc.heat_load_control__edg_total_ref_area
label: _edg_total_ref_area
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_total_ref_area
  lines:
  - 12
  - 12
inputs:
- id: spacecraft
  type: Any
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
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
  description: Return value of `_edg_total_ref_area`.
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

# _edg_total_ref_area

## Purpose
Sums the reference areas of every spacecraft link to obtain the total aerodynamic reference area used to normalise weighted coefficients and drag acceleration in the heat-load predictor.

## Design & Implementation
`@inline` loop over `spacecraft.links` converting each `link.ref_area` to `Float64` and adding it only when finite and strictly positive. Returns `max(area, eps(Float64))` so the result is always usable as a divisor. Called by `_edg_weighted_aero_coefficients`, `_edg_integrated_heat_load_trajectory`, and `_edg_heat_load_profile_for_k`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_edg_total_ref_area`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_heat_load_profile_for_k|_edg_heat_load_profile_for_k]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:558-558`
- [[gnc.heat_load_control__edg_integrated_heat_load_trajectory|_edg_integrated_heat_load_trajectory]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:186-186`
- [[gnc.heat_load_control__edg_weighted_aero_coefficients|_edg_weighted_aero_coefficients]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:51-51`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_load_control.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/heat_load_control.jl:15-15`
<!-- vulcan:connections:end -->

## Limitations
Areas are summed regardless of orientation, so links that shadow one another are double counted. A spacecraft with no valid areas silently yields `eps(Float64)` (about 2.2e-16 m^2), which makes area-normalised coefficients blow up rather than raising an error.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 12.
