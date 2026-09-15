---
id: environment.planets__furnsh_required
label: _furnsh_required
kind: function
source:
  file: src/environment/ephemerides/planets.jl
  symbol: _furnsh_required
  lines:
  - 205
  - 205
inputs:
- id: spice_path
  type: String
  units: n/a
  required: true
  description: Positional argument `spice_path`.
- id: relpath
  type: String
  units: n/a
  required: true
  description: Positional argument `relpath`.
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
  description: Return value of `_furnsh_required`. Returns `kernel_path`.
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

# _furnsh_required

## Purpose
Loads a SPICE kernel that the caller cannot proceed without, converting a missing file into an early `ArgumentError` with the full path instead of an opaque CSPICE failure later.

## Design & Implementation
`@inline` function taking `spice_path::String` and `relpath::String`. It computes `kernel_path = joinpath(spice_path, relpath)`, throws `ArgumentError("Required SPICE kernel not found: ...")` unless `isfile` is true, then calls `_furnsh_once(kernel_path)` and returns the joined path. Used for the leapseconds kernel `lsk/naif0012.tls`, the generic PCKs, and the lunar frame kernels.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `spice_path` | String | n/a | yes | Positional argument `spice_path`. |
| in | `relpath` | String | n/a | yes | Positional argument `relpath`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_furnsh_required`. Returns `kernel_path`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.planets__furnsh_mars_pck|_furnsh_mars_pck]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:334-334`
- [[environment.planets__furnsh_planetary_kernel|_furnsh_planetary_kernel]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:296-296`
- [[environment.planets__spice_body_fixed_frame|_spice_body_fixed_frame]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:375-375`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- `callees` → [[environment.planets__furnsh_once|_furnsh_once]] · `callers` · call · `src/environment/ephemerides/planets.jl:208-208`
<!-- vulcan:connections:end -->

## Limitations
`isfile` follows symlinks and does not distinguish a zero-byte file from a valid one. The returned path is the joined (not absolute) form even though `_furnsh_once` records the absolute path internally, so callers comparing paths must normalise them.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 205.
