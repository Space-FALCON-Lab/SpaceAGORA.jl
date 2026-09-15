---
id: simulation.interpolation__gram_track_cache_enabled
label: _gram_track_cache_enabled
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/interpolation.jl
  symbol: _gram_track_cache_enabled
  lines:
  - 1
  - 1
inputs:
- id: cfg
  type: GramTrackCacheConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: density_model
  type: Any
  units: n/a
  required: true
  description: Positional argument `density_model`.
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
  type: Bool
  units: n/a
  description: Return value of `_gram_track_cache_enabled`.
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

# _gram_track_cache_enabled

## Purpose
Decides whether the GRAM ground-track cache should be consulted at all for the current configuration and density model, gating the entire cached-atmosphere fast path.

## Design & Implementation
Takes `cfg::GramTrackCacheConfig` and an untyped `density_model`, returning `Bool`. It short-circuits to `false` when `cfg.mode == :off`, so the cache can be disabled purely by configuration without touching the model. Otherwise the answer is `_is_gram_density_model(density_model)`: the cache stores GRAM-specific outputs — density, temperature, and a three-component wind vector — so it is only valid when the underlying atmosphere really is a GRAM model.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | GramTrackCacheConfig | n/a | yes | Positional argument `cfg`. |
| in | `density_model` | Any | n/a | yes | Positional argument `density_model`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_gram_track_cache_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.runtime__density_state_from_kinematics_bang|_density_state_from_kinematics!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:118-118`
- [[simulation.runtime_affect_bang|affect!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:262-262`
- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:262-262`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the `:off` mode is handled explicitly; any other symbol value of `cfg.mode`, including a typo, is treated as enabled, so misconfiguration silently activates the cache rather than erroring. The function re-runs the model type test on every call rather than caching the answer, and the result depends on `_is_gram_density_model` correctly recognising every GRAM wrapper type — a newly added GRAM variant that this predicate does not know about disables the cache without any diagnostic.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/interpolation.jl` line 1.
