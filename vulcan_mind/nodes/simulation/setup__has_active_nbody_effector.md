---
id: simulation.setup__has_active_nbody_effector
label: _has_active_nbody_effector
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _has_active_nbody_effector
  lines:
  - 1464
  - 1464
inputs:
- id: dynamic_effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
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
  type: Bool
  units: n/a
  description: Return value of `_has_active_nbody_effector`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _has_active_nbody_effector

## Purpose
Tests whether the configuration includes at least one N-body effector with a non-empty body list, gating ephemeris cache construction.

## Design & Implementation
Loops the effector tuple returning true on the first `_is_nbody_effector_like` effector whose `body_names` is non-empty. `@inline` with a `::Bool` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_has_active_nbody_effector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_nbody_ephemeris_cache_bang|_initialize_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1803-1803`
- [[simulation.setup__prewarm_nbody_ephemeris_cache|_prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1688-1688`
- [[simulation.setup__validate_ephemerides_support_bang|_validate_ephemerides_support!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:111-111`

**Downstream**

- `callees` → [[simulation.setup__is_nbody_effector_like|_is_nbody_effector_like]] · `callers` · call · `src/simulation/engine/setup.jl:1466-1466`
<!-- vulcan:connections:end -->

## Limitations
An N-body effector configured with no bodies is treated as absent, which is correct for caching but means a misconfigured effector is silently inert.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1464.
