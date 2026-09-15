---
id: simulation.targeting__gram_entry_reference_density
label: _gram_entry_reference_density
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/targeting.jl
  symbol: _gram_entry_reference_density
  lines:
  - 218
  - 218
inputs:
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: h
  type: Float64
  units: n/a
  required: true
  description: Positional argument `h`.
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
  description: Return value of `_gram_entry_reference_density`.
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

# _gram_entry_reference_density

## Purpose
Cheap exponential-atmosphere density estimate used inside the Allen-Eggers entry integrator so the target predictor never calls the expensive GRAM model it is trying to pre-position the cache for.

## Theory & Math
$\rho(h) = \rho_{ref}\,\exp\!\left(\dfrac{h_{ref} - h}{H}\right)$ with $H = \max(1, H_{planet})$ metres.

## Design & Implementation
Computes `H = max(1.0, planet.H)` (scale height in metres, floored to avoid division by zero) and `ρ = planet.ρ_ref * exp((planet.h_ref - h) / H)`. Returns `ρ` when finite and positive, otherwise `0.0`. The function is `@inline` with `Float64` arguments and return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `h` | Float64 | n/a | yes | Positional argument `h`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_gram_entry_reference_density`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.targeting_gram_entry_target_allen_eggers|_gram_entry_target_allen_eggers]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/targeting.jl:277-277`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A single-scale-height exponential is a crude representation of any real atmosphere; errors of a factor of several are typical across a 100 km altitude range, so the predicted track endpoint is approximate. Negative altitudes far below `h_ref` can overflow `exp` to `Inf`, which the guard converts to `0.0`, an unphysical result that silently removes drag. Depends on `planet.ρ_ref`, `planet.h_ref`, and `planet.H` being consistently in kg/m³ and metres.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/targeting.jl` line 218.
