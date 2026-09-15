---
id: simulation.targeting__gram_entry_mass_kg
label: _gram_entry_mass_kg
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/targeting.jl
  symbol: _gram_entry_mass_kg
  lines:
  - 205
  - 205
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
- id: current_mass_kg
  type: Float64
  units: n/a
  required: true
  description: Positional argument `current_mass_kg`.
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
  description: Return value of `_gram_entry_mass_kg`.
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

# _gram_entry_mass_kg

## Purpose
Supplies a positive spacecraft mass for the Allen-Eggers entry-target predictor, preferring the live propagated mass and falling back to the configured dry-plus-propellant mass or a 100 kg default.

## Design & Implementation
If `current_mass_kg` is finite and positive it is returned immediately. Otherwise, inside a `try`, the function reads `p.args.dynamics_model.spacecraft[sat_index]` and computes `Float64(spacecraft.dry_mass + spacecraft.prop_mass)`, returning it when finite and positive and `100.0` otherwise. Any exception (missing satellite index, missing field) is caught and `100.0` is returned. The function is `@inline` with a `Float64` return annotation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_index` | Int | n/a | yes | Positional argument `sat_index`. |
| in | `current_mass_kg` | Float64 | n/a | yes | Positional argument `current_mass_kg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_gram_entry_mass_kg`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:106-106`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:211-211`
<!-- vulcan:connections:end -->

## Limitations
The 100 kg fallback is an arbitrary hard-coded constant that can badly mis-predict ballistic coefficient for large or small vehicles when configuration data is unavailable. The bare `catch` hides programming errors such as a wrong field name. Units are assumed kilograms throughout with no verification. Because the live mass is trusted whenever positive, a corrupted small value (for example 1e-300) would pass the check.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/targeting.jl` line 205.
