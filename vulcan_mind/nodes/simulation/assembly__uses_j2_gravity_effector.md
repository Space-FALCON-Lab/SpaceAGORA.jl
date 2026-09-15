---
id: simulation.assembly__uses_j2_gravity_effector
label: _uses_j2_gravity_effector
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/assembly.jl
  symbol: _uses_j2_gravity_effector
  lines:
  - 10
  - 10
inputs:
- id: effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `effectors`.
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
  description: Return value of `_uses_j2_gravity_effector`.
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

# _uses_j2_gravity_effector

## Purpose
Reports whether the run's effector stack includes the oblateness-corrected gravity model, letting callback assembly and solver policy specialise on the presence of J2 perturbations.

## Design & Implementation
Mirrors `_uses_atmospheric_dynamic_effector` in structure: an `@inline` function that walks the `effectors::Tuple` under `@inbounds` and returns `true` at the first element satisfying `effector isa InverseSquaredJ2GravityModel`, otherwise `false`. The single-type test keeps it constant-foldable for a concretely typed tuple.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effectors` | Tuple | n/a | yes | Positional argument `effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_uses_j2_gravity_effector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.runtime_get_density_callback|get_density_callback]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:221-221`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only that one concrete type is recognised; a higher-order geopotential model, or a J2 implementation wrapped in a decorator, returns `false` even though oblateness is present. The predicate says nothing about whether the J2 term is actually enabled inside the model instance — a model configured with a zero J2 coefficient still reports `true`. The linear scan is recomputed at each call rather than cached on the configuration.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/assembly.jl` line 10.
