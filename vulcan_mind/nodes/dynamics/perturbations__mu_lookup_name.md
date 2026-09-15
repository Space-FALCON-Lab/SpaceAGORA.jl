---
id: dynamics.perturbations__mu_lookup_name
label: _mu_lookup_name
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _mu_lookup_name
  lines:
  - 57
  - 57
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
  description: Return value of `_mu_lookup_name`. Returns `replace(_canonical_spice_name(name),
    "_barycenter" => "")`.
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

# _mu_lookup_name

## Purpose
Derives the key for the third-body GM table from a body name by stripping any barycentre suffix.

## Design & Implementation
Applies `_canonical_spice_name` then removes `_barycenter`. `@inline`. The GM table is keyed by planet, while SPICE queries for outer planets use the barycentre.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_mu_lookup_name`. Returns `replace(_canonical_spice_name(name), "_barycenter" => "")`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__resolve_third_body_mu|_resolve_third_body_mu]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:66-66`

**Downstream**

- `callees` → [[dynamics.perturbations__canonical_spice_name|_canonical_spice_name]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:57-57`
<!-- vulcan:connections:end -->

## Limitations
Assumes the barycentre GM equals the planet GM, which ignores the moons' mass — acceptable at the perturbation level.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 57.
