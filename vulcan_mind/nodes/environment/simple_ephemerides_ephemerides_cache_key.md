---
id: environment.simple_ephemerides_ephemerides_cache_key
label: ephemerides_cache_key
kind: function
source:
  file: src/environment/ephemerides/simple_ephemerides.jl
  symbol: ephemerides_cache_key
  lines:
  - 122
  - 122
inputs:
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
  type: Tuple
  units: n/a
  description: Return value of `ephemerides_cache_key`. Returns `(:spice,)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- environment
origin: agent
---

# ephemerides_cache_key

## Purpose
Produces a hashable, value-based key that identifies an ephemeris model configuration, so that caches keyed on ephemeris state (planet frame caches, n-body ephemeris caches) are reused only across simulations with equivalent ephemeris settings.

## Design & Implementation
For `SpiceEphemeridesModel` the key is the constant tuple `(:spice,)` because the model has no fields. For `SimpleEphemeridesModel` it is `(:simple, round(Int64, reference_epoch_seconds * 1e6), pm_key)` where `pm_key` is `typemin(Int64)` when `prime_meridian_at_reference_rad` is `NaN` (the planet-true default sentinel, since `round` would throw on NaN) and otherwise `round(Int64, pm * 1e12)`. Quantising to microseconds and picoradians makes floating-point noise in equal configurations hash identically. Both methods are `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `ephemerides_cache_key`. Returns `(:spice,)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__harmonics_lpi_cache_key|_harmonics_lpi_cache_key]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:563-563`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/simple_ephemerides.jl`
- [[simulation.setup__ephemerides_model_reuse_key|_ephemerides_model_reuse_key]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:263-263`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`round(Int64, x * 1e6)` overflows for `reference_epoch_seconds` beyond about 9.2e12 s (roughly 290,000 years), throwing `InexactError`. Two prime meridians differing by less than 1e-12 rad collide, which is harmless. The key does not include the SPICE kernel set loaded, so two `SpiceEphemeridesModel` runs with different kernels share `(:spice,)` and may reuse stale cached frames.

## Provenance
Mapped from `src/environment/ephemerides/simple_ephemerides.jl` line 122.
