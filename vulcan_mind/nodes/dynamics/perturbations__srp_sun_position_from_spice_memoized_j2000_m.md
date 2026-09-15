---
id: dynamics.perturbations__srp_sun_position_from_spice_memoized_j2000_m
label: _srp_sun_position_from_spice_memoized_j2000_m
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _srp_sun_position_from_spice_memoized_j2000_m
  lines:
  - 1424
  - 1424
inputs:
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
  description: Return value of `_srp_sun_position_from_spice_memoized_j2000_m`.
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

# _srp_sun_position_from_spice_memoized_j2000_m

## Purpose
Serves the Sun position from the single-epoch memo under the key `sun`, querying and storing on a miss.

## Design & Implementation
Under the memo lock, resets if stale, returns the cached entry if present, else counts, queries, stores and returns. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `primary_body_name` | String | n/a | yes | Positional argument `primary_body_name`. |
| in | `memo` | SpiceRhsMemo | n/a | yes | Positional argument `memo`. |
| in | `counter` | Base.Threads.Atomic{Int64} | n/a | yes | Positional argument `counter`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_srp_sun_position_from_spice_memoized_j2000_m`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations__srp_sun_position_from_spice_j2000_m|_srp_sun_position_from_spice_j2000_m]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1411-1411`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- `callees` → [[dynamics.perturbations__spice_rhs_memo_reset_if_stale_bang|_spice_rhs_memo_reset_if_stale!]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1432-1432`
- `callees` → [[environment.simple_ephemerides_spice_position_j2000_m|spice_position_j2000_m]] · `callers` · call · `src/dynamics/coupled/perturbations.jl:1437-1437`
<!-- vulcan:connections:end -->

## Limitations
Shares the memo with third bodies, so an N-body query at a different stage time evicts the Sun entry.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1424.
