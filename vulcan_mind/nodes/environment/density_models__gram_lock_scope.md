---
id: environment.density_models__gram_lock_scope
label: _gram_lock_scope
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _gram_lock_scope
  lines:
  - 416
  - 416
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
  type: Symbol
  units: n/a
  description: Return value of `_gram_lock_scope`.
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

# _gram_lock_scope

## Purpose
Reads whether native GRAM calls serialise on the process-wide lock or on each model instance's own lock.

## Design & Implementation
Parses `SPACEAGORA_GRAM_LOCK_SCOPE`, returning `:global` for empty or `global` and `:model` for `model`, `per_model`, `per-model` or `instance`, raising `ArgumentError` otherwise. `@inline`. Model scope relies on the same premise as the isolated-pool batch path: independent instances may run concurrently as long as each is serialised.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_gram_lock_scope`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__nrlmsise_flux_value|_nrlmsise_flux_value]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:404-404`
- [[ext.spaceagoragramsuiteext__gram_call_lock|_gram_call_lock]] · `callees` → `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:18-18`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It reads the environment on every call, and model scope is unsafe if two wrappers share native state — the function cannot detect that.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 416.
