---
id: dynamics.perturbations__resolve_third_body_mu
label: _resolve_third_body_mu
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _resolve_third_body_mu
  lines:
  - 65
  - 65
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
  type: Float64
  units: n/a
  description: Return value of `_resolve_third_body_mu`.
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

# _resolve_third_body_mu

## Purpose
Looks up a third body's gravitational parameter from the built-in table at model construction, failing clearly for unknown bodies.

## Design & Implementation
Derives the lookup key, reads `_THIRD_BODY_MU` with a `NaN` default, and raises `ArgumentError` directing the user to add the GM if it is not finite. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Positional argument `name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_resolve_third_body_mu`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations_solarradiationpressuremodel|SolarRadiationPressureModel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:635-635`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[dynamics.perturbations__mu_lookup_name|_mu_lookup_name]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:66-66`
<!-- vulcan:connections:end -->

## Limitations
The table is a hard-coded constant in this file; there is no way to supply a GM from configuration.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 65.
