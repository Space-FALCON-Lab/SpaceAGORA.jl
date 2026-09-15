---
id: simulation.setup__is_nbody_effector_like
label: _is_nbody_effector_like
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _is_nbody_effector_like
  lines:
  - 358
  - 358
inputs:
- id: effector
  type: Any
  units: n/a
  required: true
  description: Positional argument `effector`.
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
  description: Return value of `_is_nbody_effector_like`.
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

# _is_nbody_effector_like

## Purpose
Duck-types an effector as an N-body gravity source, accepting both the built-in `NBodyGravityModel` and any user type exposing `body_names` and `primary_body_name` properties, so ephemeris cache setup covers custom third-body models.

## Design & Implementation
Returns `effector isa SimulationModel.NBodyGravityModel || (hasproperty(effector, :body_names) && hasproperty(effector, :primary_body_name))`. The `hasproperty` checks are evaluated at runtime on the concrete value, so structs, named tuples, and objects with custom `propertynames` all qualify. Pure and `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effector` | Any | n/a | yes | Positional argument `effector`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_is_nbody_effector_like`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__collect_nbody_query_names|_collect_nbody_query_names]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1477-1477`
- [[simulation.setup__has_active_nbody_effector|_has_active_nbody_effector]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1466-1466`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Property presence is the only test; the types of `body_names` (expected `Vector{String}`) and `primary_body_name` (expected `String`) are not verified, so a false positive fails later inside `_collect_nbody_query_names` with a `MethodError`. It also does not check whether the effector is active (non-empty body list).

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 358.
