---
id: gnc.heat_load_control__edg_integrate_series
label: _edg_integrate_series
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_integrate_series
  lines:
  - 414
  - 414
inputs:
- id: time
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `time`.
- id: values
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `values`.
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
  description: Return value of `_edg_integrate_series`. Returns `total`.
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

# _edg_integrate_series

## Purpose
Trapezoidal integration of a sampled time series, used to convert heat-rate histories (W/cm^2) into accumulated heat load (J/cm^2).

## Theory & Math
$$Q = \sum_{j=1}^{n-1} \tfrac{1}{2}\,(q_j + q_{j+1})\,(t_{j+1} - t_j)$$ where $q_j$ are the sampled values and $t_j$ the sample times (s).

## Design & Implementation
Takes `time::Vector{Float64}` and `values::Vector{Float64}` of equal length and accumulates `0.5 (values[j] + values[j+1]) (time[j+1] - time[j])` for `j in 1:length(time)-1`. Returns the scalar total. Handles non-uniform spacing correctly because each interval width is taken from the time vector.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `time` | Vector{Float64} | n/a | yes | Positional argument `time`. |
| in | `values` | Vector{Float64} | n/a | yes | Positional argument `values`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_integrate_series`. Returns `total`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.heat_load_control__edg_balanced_tpbvp_heat_load_window|_edg_balanced_tpbvp_heat_load_window]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:488-488`
- [[gnc.heat_load_control__edg_profile_heat_load|_edg_profile_heat_load]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:424-424`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/heat_load_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No length check between `time` and `values`; a shorter `values` vector raises `BoundsError`. Trapezoidal accuracy is second order, so sharp heating peaks between 1 s samples are underestimated. An empty or single-element series returns 0.0.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 414.
