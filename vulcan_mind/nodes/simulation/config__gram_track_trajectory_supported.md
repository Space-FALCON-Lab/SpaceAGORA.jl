---
id: simulation.config__gram_track_trajectory_supported
label: _gram_track_trajectory_supported
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _gram_track_trajectory_supported
  lines:
  - 154
  - 154
inputs:
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
  description: Return value of `_gram_track_trajectory_supported`.
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

# _gram_track_trajectory_supported

## Purpose
Capability probe answering whether the configured density model can generate a GRAM trajectory, which is the precondition for using the track-cache path instead of point-by-point atmospheric queries.

## Design & Implementation
Returns `false` immediately unless `_is_gram_density_model(density_model)` holds, i.e. the model is a `GRAMAtmosphereModel` or a `GRAMAtmosphereModelSurrogate`. It then requires both `hasproperty(density_model, :gram)` and `hasproperty(density_model, :gram_atmosphere)`. The `:gram` property is fetched inside a `try`, returning `false` on any error. Because `hasproperty` on a `Module` only sees exported names, a `Module` driver is probed with `isdefined(gram_driver, :generate_trajectory)`; a non-module driver is probed with `hasproperty`. The source comment notes this matches the capability check GRAMSuite itself uses for entry points such as `get_winds_state`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `density_model` | Any | n/a | yes | Positional argument `density_model`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_gram_track_trajectory_supported`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:195-195`

**Downstream**

- `callees` → [[simulation.config__callback_env_config|_callback_env_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:175-175`
- `callees` → [[simulation_a.density_config_snapshot_callback_env_config|_snapshot_callback_env_config]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:170-170`
<!-- vulcan:connections:end -->

## Limitations
The probe confirms only that the name `generate_trajectory` exists, not that it has a usable signature or that the underlying native library is loaded, so a mismatched GRAMSuite version still fails later at call time. The bare `catch` around the property fetch swallows every error class, including genuine initialisation faults, and reports them as an absent capability.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 154.
