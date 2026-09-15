---
id: environment.density_models__nrlmsise_provider_indices
label: _nrlmsise_provider_indices
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _nrlmsise_provider_indices
  lines:
  - 565
  - 565
inputs:
- id: result
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `result`.
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
  description: Return value of `_nrlmsise_provider_indices`. Returns `(`.
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

# _nrlmsise_provider_indices

## Purpose
Normalises whatever a user-supplied index provider returned into validated `(f107a, f107, ap)`, accepting tuple or named-tuple forms with either lowercase or STK-style uppercase names.

## Design & Implementation
Three methods. A `Tuple` must have length three and is validated positionally. A `NamedTuple` may use `f107a` or `F10A`, `f107` or `F10`, and `ap` or `Ap`, with a specific error for each missing key. Any other type raises an error describing the accepted forms.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `result` | Tuple | n/a | yes | Positional argument `result`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_nrlmsise_provider_indices`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models__nrlmsise_resolved_indices|_nrlmsise_resolved_indices]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:635-635`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/atmosphere/density_models.jl`

**Downstream**

- `callees` → [[environment.density_models__nrlmsise_ap_value|_nrlmsise_ap_value]] · `callers` · call · `src/environment/atmosphere/density_models.jl:572-572`
- `callees` → [[environment.density_models__nrlmsise_flux_value|_nrlmsise_flux_value]] · `callers` · call · `src/environment/atmosphere/density_models.jl:570-570`
- `callees` → [[environment.density_models__nrlmsise_space_indices_indices|_nrlmsise_space_indices_indices]] · `callers` · call · `src/environment/atmosphere/density_models.jl:612-612`
- `callees` → [[environment.density_models_init_nrlmsise_space_indices_bang|init_nrlmsise_space_indices!]] · `callers` · call · `src/environment/atmosphere/density_models.jl:610-610`
<!-- vulcan:connections:end -->

## Limitations
Mixed-case variants such as `F107a` are not recognised; the provider contract is documented only in the error messages.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 565.
