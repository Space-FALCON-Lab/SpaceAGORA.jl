---
id: dynamics.perturbations__nbody_body_position_from_spice_j2000_m
label: _nbody_body_position_from_spice_j2000_m
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _nbody_body_position_from_spice_j2000_m
  lines:
  - 1029
  - 1029
inputs:
- id: body_name_spice
  type: String
  units: n/a
  required: true
  description: Positional argument `body_name_spice`.
- id: et
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et`.
- id: primary_body_name
  type: String
  units: n/a
  required: true
  description: Positional argument `primary_body_name`.
- id: memo_enabled
  type: Bool
  units: n/a
  required: true
  description: Positional argument `memo_enabled`.
- id: memo
  type: SpiceRhsMemo
  units: n/a
  required: true
  description: Positional argument `memo`.
- id: counter
  type: Base.Threads.Atomic{Int64}
  units: n/a
  required: true
  description: Positional argument `counter`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `_nbody_body_position_from_spice_j2000_m`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# _nbody_body_position_from_spice_j2000_m

## Purpose
Routes a third-body position request to the memoised or the direct SPICE path according to the run's `spice_rhs_memo_enabled` flag, so the effector code does not branch on the flag itself.

## Design & Implementation
A conditional expression returning `_nbody_body_position_from_spice_memoized_j2000_m` when `memo_enabled` and `_nbody_body_position_from_spice_direct_j2000_m` otherwise, forwarding the memo object and counter. Declared `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `body_name_spice` | String | n/a | yes | Positional argument `body_name_spice`. |
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `primary_body_name` | String | n/a | yes | Positional argument `primary_body_name`. |
| in | `memo_enabled` | Bool | n/a | yes | Positional argument `memo_enabled`. |
| in | `memo` | SpiceRhsMemo | n/a | yes | Positional argument `memo`. |
| in | `counter` | Base.Threads.Atomic{Int64} | n/a | yes | Positional argument `counter`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_nbody_body_position_from_spice_j2000_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:937-937`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[dynamics.perturbations__nbody_body_position_from_spice_direct_j2000_m|_nbody_body_position_from_spice_direct_j2000_m]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1039-1039`
- `callees` → [[dynamics.perturbations__nbody_body_position_from_spice_memoized_j2000_m|_nbody_body_position_from_spice_memoized_j2000_m]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1038-1038`
<!-- vulcan:connections:end -->

## Limitations
None beyond the routing choice; both branches take a lock — the memo's or SPICE's — so neither is lock-free.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1029.
