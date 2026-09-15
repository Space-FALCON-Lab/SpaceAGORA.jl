---
id: ext.spaceagoragramsuiteext__gram_call_lock
label: _gram_call_lock
kind: function
source:
  file: ext/SpaceAGORAGRAMSuiteExt.jl
  symbol: _gram_call_lock
  lines:
  - 17
  - 17
inputs:
- id: model
  type: EM.GRAMAtmosphereModel
  units: n/a
  required: true
  description: Positional argument `model`.
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
  type: ReentrantLock
  units: n/a
  description: Return value of `_gram_call_lock`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- ext
charts:
- ext
origin: agent
---

# _gram_call_lock

## Purpose
Chooses which lock serialises native GRAM calls for a given model, trading process-wide safety against per-instance concurrency.

## Design & Implementation
Reads `EM._gram_lock_scope()` and returns `model.instance_lock` when it is `:model`, otherwise the process-wide `GRAM_LOCK`. Under the default global scope every GRAM call in the process contends on one lock; under model scope only calls on the same wrapper instance serialise, which is what lets per-sample or per-worker model copies evaluate density concurrently. Declared `@inline` with a `::ReentrantLock` return because it is consulted on every density evaluation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | EM.GRAMAtmosphereModel | n/a | yes | Positional argument `model`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | ReentrantLock | n/a | — | Return value of `_gram_call_lock`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[ext.gram_core_density_state|_gram_core_density_state]] · `callees` → `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:264-264`
- [[ext.spaceagoragramsuiteext__gram_spice_ephemeris_state|_gram_spice_ephemeris_state]] · `callees` → `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:264-264`

**Downstream**

- `callees` → [[environment.density_models__gram_lock_scope|_gram_lock_scope]] · `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:18-18`
<!-- vulcan:connections:end -->

## Limitations
Model scope is only sound if each copy genuinely owns an isolated native instance, which this function assumes rather than verifies; selecting it for models that share native state would allow concurrent entry into the same library.

## Provenance
Mapped from `ext/SpaceAGORAGRAMSuiteExt.jl` line 17.
