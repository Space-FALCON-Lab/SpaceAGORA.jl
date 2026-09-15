---
id: environment.planets__furnsh_first_existing
label: _furnsh_first_existing
kind: function
source:
  file: src/environment/ephemerides/planets.jl
  symbol: _furnsh_first_existing
  lines:
  - 212
  - 212
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
  type: Any
  units: n/a
  description: 'Return value of `_furnsh_first_existing`. Returns `kernel_path`. Type
    parameters: `{N}`.'
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

# _furnsh_first_existing

## Purpose
Loads the first SPICE kernel from an ordered list of candidate relative paths that exists on disk, failing loudly when none is present, so planet constructors can prefer newer ephemeris files while tolerating older bundles.

## Design & Implementation
Signature `_furnsh_first_existing(spice_path::String, relpaths::NTuple{N, String})`. It iterates `relpaths` in order, forms `joinpath(spice_path, relpath)`, and on the first `isfile` hit calls `_furnsh_once` and returns that absolute-ish kernel path. If the loop completes it throws `ArgumentError` listing `spice_path` and every path tried joined by `", "`. Used by `_furnsh_planetary_kernel`, `_furnsh_mars_system_kernel`, and the Titan constructor.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `spice_path` | String | n/a | yes | Positional argument `spice_path`. |
| in | `relpaths` | NTuple{N, String} | n/a | yes | Positional argument `relpaths`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_furnsh_first_existing`. Returns `kernel_path`. Type parameters: `{N}`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.planets__furnsh_mars_system_kernel|_furnsh_mars_system_kernel]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:312-312`
- [[environment.planets__furnsh_planetary_kernel|_furnsh_planetary_kernel]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:298-298`
- [[environment.planets__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:440-440`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- `callees` → [[environment.planets__furnsh_once|_furnsh_once]] · `callers` · call · `src/environment/ephemerides/planets.jl:216-216`
<!-- vulcan:connections:end -->

## Limitations
Only file existence is checked, not readability or kernel validity; a truncated kernel file passes this check and fails later inside CSPICE. Order matters and is fixed at the call site, so a caller cannot express a preference by date without rewriting the tuple. The `NTuple` type restriction rejects a `Vector{String}` argument.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 212.
