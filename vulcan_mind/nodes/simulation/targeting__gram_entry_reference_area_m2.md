---
id: simulation.targeting__gram_entry_reference_area_m2
label: _gram_entry_reference_area_m2
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/targeting.jl
  symbol: _gram_entry_reference_area_m2
  lines:
  - 185
  - 185
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: sat_index
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_index`.
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
  description: Return value of `_gram_entry_reference_area_m2`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _gram_entry_reference_area_m2

## Purpose
Determines the aerodynamic reference area (m²) for the Allen-Eggers entry-target estimate by summing the `ref_area` of every link on the satellite, with graceful fallbacks to the root body area or 1 m².

## Design & Implementation
Inside a `try`, reads `spacecraft = p.args.dynamics_model.spacecraft[sat_index]` and accumulates `Float64(link.ref_area)` over `spacecraft.links` for entries that are finite and positive. If the sum is positive it is returned; otherwise `spacecraft.root.ref_area` is used when finite and positive, else `1.0`. Any exception yields `1.0`. The function is `@inline` and annotated `::Float64`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_index` | Int | n/a | yes | Positional argument `sat_index`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_gram_entry_reference_area_m2`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:107-107`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:190-190`
<!-- vulcan:connections:end -->

## Limitations
Summing all link areas assumes every link's reference area is simultaneously flow-facing, which overestimates area for a bus-plus-panel geometry unless panels are edge-on. Negative or NaN link areas are silently skipped rather than reported. The 1 m² fallback is a hard-coded guess. The bare `catch` masks structural errors in the configuration.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/targeting.jl` line 185.
