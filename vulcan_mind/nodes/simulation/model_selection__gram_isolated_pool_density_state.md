---
id: simulation.model_selection__gram_isolated_pool_density_state
label: _gram_isolated_pool_density_state
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/model_selection.jl
  symbol: _gram_isolated_pool_density_state
  lines:
  - 86
  - 86
inputs:
- id: model
  type: EnvironmentModels.GRAMAtmosphereModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: h
  type: Float64
  units: n/a
  required: true
  description: Positional argument `h`.
- id: lat
  type: Float64
  units: n/a
  required: true
  description: Positional argument `lat`.
- id: lon
  type: Float64
  units: n/a
  required: true
  description: Positional argument `lon`.
- id: el_time
  type: Float64
  units: n/a
  required: true
  description: Positional argument `el_time`.
- id: wind
  type: Bool
  units: n/a
  required: true
  description: Positional argument `wind`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: model_lock
  type: ReentrantLock
  units: n/a
  required: true
  description: Positional argument `model_lock`.
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
  type: Tuple{Float64,
  units: n/a
  description: Return value of `_gram_isolated_pool_density_state`.
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

# _gram_isolated_pool_density_state

## Purpose
Evaluates atmospheric density, temperature and wind for one satellite using a worker-private GRAM model instance, applying the same altitude regime switching as the serial path so results are identical regardless of pool use.

## Design & Implementation
Arguments are the worker's `model::GRAMAtmosphereModel`, geodetic altitude `h` (m), `lat` and `lon` (rad), `el_time` (s), a `wind::Bool` flag, the parameter object `p`, and the worker's `model_lock`. The entry interface altitude `EI = p.args.environment_model.EI * 1e3` (km to m) defines `drag_state = h - EI <= 0`. Above 2000 km the function returns `(0.0, planet.T_ref, zero wind)`; outside the drag regime in a non-Keplerian mission it returns `EnvironmentModels.density_polyfit(h, p)`; otherwise it calls `EnvironmentModels._gram_core_density_state(model.core, h, lat, lon, el_time, wind, model_lock, T_ref)`. The return is a `Tuple{Float64, Float64, SVector{3,Float64}}`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | EnvironmentModels.GRAMAtmosphereModel | n/a | yes | Positional argument `model`. |
| in | `h` | Float64 | n/a | yes | Positional argument `h`. |
| in | `lat` | Float64 | n/a | yes | Positional argument `lat`. |
| in | `lon` | Float64 | n/a | yes | Positional argument `lon`. |
| in | `el_time` | Float64 | n/a | yes | Positional argument `el_time`. |
| in | `wind` | Bool | n/a | yes | Positional argument `wind`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `model_lock` | ReentrantLock | n/a | yes | Positional argument `model_lock`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{Float64, | n/a | — | Return value of `_gram_isolated_pool_density_state`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.model_selection_gram_isolated_pool_batch_eval__gram_isolated_pool_batch_eval_bang|_gram_isolated_pool_batch_eval!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:208-208`

**Downstream**

- `callees` → [[environment.density_models__gram_core_density_state|_gram_core_density_state]] · `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:103-103`
- `callees` → [[environment.density_models_density_polyfit|density_polyfit]] · `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:101-101`
- `callees` → [[ext.gram_core_density_state|_gram_core_density_state]] · `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:103-103`
<!-- vulcan:connections:end -->

## Limitations
The 2000 km cutoff is a hard-coded constant with no configuration hook. `EI` is assumed to be in kilometres and is multiplied by 1e3 every call. When `mission_configuration.keplerian` is true, GRAM is consulted even above the entry interface, which may be slow. The lock is passed through but whether it is acquired depends on `_gram_core_density_state`.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/model_selection.jl` line 86.
