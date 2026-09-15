---
id: gnc.heat_load_control__edg_predict_mass
label: _edg_predict_mass
kind: function
source:
  file: src/gnc/control/heat_load_control.jl
  symbol: _edg_predict_mass
  lines:
  - 21
  - 21
inputs:
- id: spacecraft
  type: Any
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
- id: mass
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mass`.
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
  description: Return value of `_edg_predict_mass`.
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

# _edg_predict_mass

## Purpose
Chooses the spacecraft mass in kilograms for the heat-load prediction, preferring the live integrated mass state and falling back to the structural plus propellant mass when the state is invalid.

## Design & Implementation
`@inline` function `(spacecraft, mass::Float64)::Float64`. If `mass` is finite and positive it is returned unchanged. Otherwise it sums `max(0, link.m)` over `spacecraft.links`, adds `max(0, spacecraft.prop_mass)` when that property exists, and returns `max(total, 1.0)` so a degenerate model still yields at least 1 kg.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `mass` | Float64 | n/a | yes | Positional argument `mass`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_edg_predict_mass`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control__edg_predict_max_energy_depletion_outcome|_edg_predict_max_energy_depletion_outcome]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:724-724`
- [[gnc.targeting_control__edg_predict_targeting_outcome|_edg_predict_targeting_outcome]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:682-682`
- [[gnc.targeting_control__edg_solve_targeting_switch|_edg_solve_targeting_switch]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:966-966`
- [[gncx.heat_load_control__edg_solve_heat_load_switches|_edg_solve_heat_load_switches]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:632-632`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/heat_load_control.jl:27-27`
<!-- vulcan:connections:end -->

## Limitations
The 1 kg floor masks a misconfigured spacecraft rather than failing. Negative link masses are clamped to zero silently. The function does not distinguish between a NaN mass caused by integrator failure and a deliberately zeroed mass.

## Provenance
Mapped from `src/gnc/control/heat_load_control.jl` line 21.
