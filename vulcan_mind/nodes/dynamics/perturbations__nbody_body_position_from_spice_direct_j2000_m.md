---
id: dynamics.perturbations__nbody_body_position_from_spice_direct_j2000_m
label: _nbody_body_position_from_spice_direct_j2000_m
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _nbody_body_position_from_spice_direct_j2000_m
  lines:
  - 1042
  - 1042
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
  description: Return value of `_nbody_body_position_from_spice_direct_j2000_m`.
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

# _nbody_body_position_from_spice_direct_j2000_m

## Purpose
Performs a live SPICE position query for one third body relative to the primary and records it in the runtime call counter, the path taken when the ephemeris cache misses and the memo is disabled.

## Design & Implementation
Increments `counter` with `Threads.atomic_add!` and returns `spice_position_j2000_m(body_name_spice, et, primary_body_name)`. Declared `@inline`. Counting before the call means an exception inside SPICE still leaves the attempt recorded.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `body_name_spice` | String | n/a | yes | Positional argument `body_name_spice`. |
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `primary_body_name` | String | n/a | yes | Positional argument `primary_body_name`. |
| in | `counter` | Base.Threads.Atomic{Int64} | n/a | yes | Positional argument `counter`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_nbody_body_position_from_spice_direct_j2000_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__nbody_body_position_from_spice_j2000_m|_nbody_body_position_from_spice_j2000_m]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1039-1039`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[environment.simple_ephemerides_spice_position_j2000_m|spice_position_j2000_m]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1049-1049`
<!-- vulcan:connections:end -->

## Limitations
The callee acquires the process-wide SPICE lock, so concurrent RHS evaluations on different threads serialise on every call through this path.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1042.
