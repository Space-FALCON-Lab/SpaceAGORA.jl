---
id: simulation.setup__collect_nbody_query_names
label: _collect_nbody_query_names
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _collect_nbody_query_names
  lines:
  - 1473
  - 1473
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
  type: Vector{String}
  units: n/a
  description: Return value of `_collect_nbody_query_names`.
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

# _collect_nbody_query_names

## Purpose
Gathers the distinct SPICE query names of every third body any N-body effector references, in first-seen order, to size and label the ephemeris cache.

## Design & Implementation
Iterates the effector tuple, skips non-N-body effectors via `_is_nbody_effector_like`, maps each `body_names` entry through `_spice_query_name`, and appends it if not already in a `Set`. Returns a `Vector{String}`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{String} | n/a | — | Return value of `_collect_nbody_query_names`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_nbody_ephemeris_cache_bang|_initialize_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1808-1808`
- [[simulation.setup__prewarm_nbody_ephemeris_cache|_prewarm_nbody_ephemeris_cache]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1706-1706`

**Downstream**

- `callees` → [[dynamics.perturbations__spice_query_name|_spice_query_name]] · `callers` · call · `src/simulation/engine/setup.jl:1481-1481`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/simulation/engine/setup.jl:1483-1483`
- `callees` → [[simulation.setup__is_nbody_effector_like|_is_nbody_effector_like]] · `callers` · call · `src/simulation/engine/setup.jl:1477-1477`
<!-- vulcan:connections:end -->

## Limitations
Order depends on effector and body order, so two configurations with the same bodies in different order produce different cache keys and do not share a reuse entry.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1473.
