---
id: environment.density_models__nrlmsise_space_indices_lookup
label: _nrlmsise_space_indices_lookup
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _nrlmsise_space_indices_lookup
  lines:
  - 514
  - 514
inputs:
- id: index
  type: Val
  units: n/a
  required: true
  description: Positional argument `index`.
- id: instant
  type: DateTime
  units: n/a
  required: true
  description: Positional argument `instant`.
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
  type: SpaceIndices.space_index
  units: n/a
  description: Return value of `_nrlmsise_space_indices_lookup`. Returns `SpaceIndices.space_index(index,
    instant)`.
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

# _nrlmsise_space_indices_lookup

## Purpose
The single indirection through which the built-in provider reads a space index, so tests can substitute a fake lookup.

## Design & Implementation
An `@inline` wrapper calling `SpaceIndices.space_index(index, instant)` with a `Val`-typed index selector. Every `_nrlmsise_space_indices_*` function takes the lookup as an argument rather than calling `SpaceIndices` directly.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `index` | Val | n/a | yes | Positional argument `index`. |
| in | `instant` | DateTime | n/a | yes | Positional argument `instant`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SpaceIndices.space_index | n/a | — | Return value of `_nrlmsise_space_indices_lookup`. Returns `SpaceIndices.space_index(index, instant)`. |
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
It is only as thread-safe as `SpaceIndices` itself, which the caller does not guard.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 514.
