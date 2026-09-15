---
id: environment.planets__furnsh_once
label: _furnsh_once
kind: function
source:
  file: src/environment/ephemerides/planets.jl
  symbol: _furnsh_once
  lines:
  - 42
  - 42
inputs:
- id: kernel_path
  type: String
  units: n/a
  required: true
  description: Positional argument `kernel_path`.
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
  description: Return value of `_furnsh_once`. Returns `nothing`.
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

# _furnsh_once

## Purpose
Idempotent wrapper around SPICE `furnsh` that guarantees each kernel path is loaded at most once per process, preventing CSPICE kernel-table exhaustion when planet constructors are called thousands of times in Monte Carlo loops.

## Design & Implementation
Takes `kernel_path::String`, resolves it with `abspath`, then under `lock(SPICE_LOCK)` checks membership in the module-global `Set{String}` `_FURNISHED_KERNELS`. On a miss it calls `furnsh(resolved)` and pushes the resolved path into the set. Both branches return `nothing`. The lock is the shared `RuntimeServices.SPICE_LOCK` reentrant lock, so callers already holding it (the planet constructors) do not deadlock.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `kernel_path` | String | n/a | yes | Positional argument `kernel_path`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_furnsh_once`. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.planets__furnsh_first_existing|_furnsh_first_existing]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:216-216`
- [[environment.planets__furnsh_first_existing_if_available|_furnsh_first_existing_if_available]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:227-227`
- [[environment.planets__furnsh_mars_pck|_furnsh_mars_pck]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:327-327`
- [[environment.planets__furnsh_required|_furnsh_required]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:208-208`
- [[environment.planets__gravity_constants_kernel_if_available|_gravity_constants_kernel_if_available]] · `callees` → `callers` · call · `src/environment/ephemerides/planets.jl:263-263`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/planets.jl`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/environment/ephemerides/planets.jl:47-47`
<!-- vulcan:connections:end -->

## Limitations
Deduplication is by absolute path string, so the same kernel reached via a symlink or a differently normalised path is loaded twice. The set is not synchronised with `kclear()`; anyone clearing the pool must also call `_reset_furnished_kernels!` or later loads will be skipped against an empty pool. Errors thrown by `furnsh` propagate before the path is recorded, so a retry will attempt the load again.

## Provenance
Mapped from `src/environment/ephemerides/planets.jl` line 42.
