---
id: environment.planets__furnsh_first_existing_if_available
label: _furnsh_first_existing_if_available
kind: function
source:
  file: src/environment/ephemerides/planets.jl
  symbol: _furnsh_first_existing_if_available
  lines:
  - 223
  - 223
inputs:
- id: spice_path
  type: String
  units: n/a
  required: true
  description: Positional argument `spice_path`.
- id: relpaths
  type: NTuple{N, String}
  units: n/a
  required: true
  description: Positional argument `relpaths`.
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
  type: Nothing
  units: n/a
  description: 'Return value of `_furnsh_first_existing_if_available`. Returns `kernel_path`
    or `nothing`. Type parameters: `{N}`.'
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

# _furnsh_first_existing_if_available

## Purpose
Optional variant of the first-existing kernel loader that returns `nothing` instead of throwing when none of the candidates are present, used for high-precision Earth orientation and ITRF93 frame kernels that the in-repo starter bundle may omit.

## Design & Implementation
Same loop as `_furnsh_first_existing`: for each `relpath` in the `NTuple{N, String}`, build `joinpath(spice_path, relpath)`, and on the first `isfile` match call `_furnsh_once(kernel_path)` and return the path. The only difference is the fall-through `return nothing`. The Earth constructor calls it twice, once for `earth_latest_high_prec.bpc`/`earth_200101_990628_predict.bpc` and once for `tf/earth_assoc_itrf93.tf`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `spice_path` | String | n/a | yes | Positional argument `spice_path`. |
| in | `relpaths` | NTuple{N, String} | n/a | yes | Positional argument `relpaths`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_furnsh_first_existing_if_available`. Returns `kernel_path` or `nothing`. Type parameters: `{N}`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.planets__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:382-382`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- `callees` → [[environment.planets__furnsh_once|_furnsh_once]] · `callers` · call · `src/environment/ephemerides/planets.jl:227-227`
<!-- vulcan:connections:end -->

## Limitations
Silent fallback means a misconfigured `spice_path` degrades Earth frame accuracy (to the generic IAU_EARTH model) with no log message. Return type is `Union{Nothing, String}`, and no caller currently inspects it.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 223.
