---
id: environment.density_models__gram_default_surrogate_file
label: _gram_default_surrogate_file
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _gram_default_surrogate_file
  lines:
  - 434
  - 434
inputs:
- id: planet
  type: String
  units: n/a
  required: true
  description: Positional argument `planet`.
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
  type: String
  units: n/a
  description: Return value of `_gram_default_surrogate_file`.
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

# _gram_default_surrogate_file

## Purpose
Resolves the default offline surrogate file for a planet by calling through the slot the GRAMSuite extension fills at load.

## Design & Implementation
Invokes `_GRAM_DEFAULT_SURROGATE_FILE_FN[](planet)`, whose default returns an empty string until the extension installs the real resolver. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | String | n/a | yes | Positional argument `planet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_gram_default_surrogate_file`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Without the extension it silently returns `""`, and callers must treat the empty string as absence rather than a path.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 434.
