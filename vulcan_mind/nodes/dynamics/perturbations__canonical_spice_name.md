---
id: dynamics.perturbations__canonical_spice_name
label: _canonical_spice_name
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _canonical_spice_name
  lines:
  - 56
  - 56
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Positional argument `name`.
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
  type: Any
  units: n/a
  description: Return value of `_canonical_spice_name`. Returns `replace(lowercase(strip(name)),
    ' ' => '_')`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# _canonical_spice_name

## Purpose
Normalises a body name to the lowercase, underscore-joined form used as the key for SPICE queries and the GM table.

## Design & Implementation
Strips, lowercases and replaces spaces with underscores. `@inline`. `_mu_lookup_name` further strips a `_barycenter` suffix for GM lookup.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_canonical_spice_name`. Returns `replace(lowercase(strip(name)), ' ' => '_')`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__mu_lookup_name|_mu_lookup_name]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:57-57`
- [[dynamics.perturbations__spice_query_name|_spice_query_name]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:59-59`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No validation that the result is a body SPICE knows; that surfaces as a SPICE error at query time.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 56.
