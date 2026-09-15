---
id: environment.density_models_gramatmospheremodel
label: GRAMAtmosphereModel
kind: struct
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: GRAMAtmosphereModel
  lines:
  - 163
  - 163
inputs:
- id: core
  type: Any
  units: n/a
  required: true
  description: Field `core`.
- id: instance_lock
  type: ReentrantLock
  units: n/a
  required: true
  description: Field `instance_lock`.
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
  type: GRAMAtmosphereModel
  units: n/a
  description: Constructed `GRAMAtmosphereModel`.
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

# GRAMAtmosphereModel

## Purpose
The core package's wrapper around a native GRAMSuite atmosphere model, preserving `AbstractDensityModel` dispatch while the real implementation lives in the package extension.

## Design & Implementation
An immutable struct with an untyped `core` holding the GRAMSuite object and a `ReentrantLock` `instance_lock` that serialises native calls on this instance when `SPACEAGORA_GRAM_LOCK_SCOPE=model`. A single-argument constructor creates a fresh lock, which `deepcopy_internal` and deserialization in the extension also rely on. `getproperty` and `propertynames` forward unknown names to `core`, so `model.planet_name` works transparently. Constructors, density evaluation, precompute and serialization are all provided by `SpaceAGORAGRAMSuiteExt`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `core` | Any | n/a | yes | Field `core`. |
| in | `instance_lock` | ReentrantLock | n/a | yes | Field `instance_lock`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | GRAMAtmosphereModel | n/a | — | Constructed `GRAMAtmosphereModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[ext.spaceagoragramsuiteext__gram_spice_ephemeris_state|_gram_spice_ephemeris_state]] · `callees` → `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:122-122`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`
- [[parallel.worker_pool__warm_gram_wrapper_bang|_warm_gram_wrapper!]] · `callees` → `callers` · call · `src/parallel/process/worker_pool.jl:153-153`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`core` is untyped, so every field access through the forwarding `getproperty` is dynamic; without the extension loaded, constructing one is impossible and evaluating one raises the not-loaded error.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 163.
