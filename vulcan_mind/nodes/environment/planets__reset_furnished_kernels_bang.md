---
id: environment.planets__reset_furnished_kernels_bang
label: _reset_furnished_kernels!
kind: function
source:
  file: src/environment/ephemerides/planets.jl
  symbol: _reset_furnished_kernels!
  lines:
  - 30
  - 30
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
  type: Nothing
  units: n/a
  description: Return value of `_reset_furnished_kernels!`. Returns `nothing`.
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

# _reset_furnished_kernels!

## Purpose
Clears the module's memory of which kernels have been furnished and drops all cached planet instances, which must be called after any `kclear()` so subsequent constructors reload kernels into the now-empty CSPICE pool.

## Design & Implementation
Under `lock(SPICE_LOCK)` it calls `empty!` on `_FURNISHED_KERNELS` and on each of `_EARTH_CACHE`, `_MARS_CACHE`, `_VENUS_CACHE`, `_TITAN_CACHE`, and `_MOON_CACHE` (all `Dict{Tuple{String,String}, T}`). Returns `nothing`. Mutates only these six module-level containers; it does not itself call `kclear()`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_reset_furnished_kernels!`. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
There is no coupling that forces this to run with `kclear()`; forgetting it produces confusing failures such as `utc2et` erroring for lack of a leapseconds kernel. Planet instances previously handed out to callers remain valid Julia objects but reference frame data that may no longer be loaded.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 30.
