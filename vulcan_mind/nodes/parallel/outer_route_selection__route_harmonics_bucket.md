---
id: parallel.outer_route_selection__route_harmonics_bucket
label: _route_harmonics_bucket
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _route_harmonics_bucket
  lines:
  - 49
  - 49
inputs:
- id: L
  type: Int
  units: n/a
  required: true
  description: Positional argument `L`.
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
  type: String
  units: n/a
  description: Return value of `_route_harmonics_bucket`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# _route_harmonics_bucket

## Purpose
Buckets the spherical-harmonics gravity degree `L` into four string classes for the `harm=` signature field, reflecting that harmonic evaluation cost grows roughly quadratically with degree.

## Design & Implementation
`@inline _route_harmonics_bucket(L::Int)::String` returns "0" for `L <= 0` (no harmonics), "1_10" for 1 to 10, "11_20" for 11 to 20, and "21p" beyond. `default_outer_route` and `_feature_heavy_for_process` separately test `f.harmonics_degree >= 20` when deciding process isolation, which straddles the "11_20"/"21p" boundary.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `L` | Int | n/a | yes | Positional argument `L`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_route_harmonics_bucket`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- [[parallel.outer_route_selection__compat_outer_route_signature|_compat_outer_route_signature]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:142-142`
- [[parallel.outer_route_selection__outer_route_signature_hierarchy|_outer_route_signature_hierarchy]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:191-191`
- [[parallel.outer_route_selection_outer_route_signature|outer_route_signature]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:163-163`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Degree 20 is treated as heavy by the routing rules but bucketed with degrees 11-19 in the signature, so feedback for degree-20 runs is pooled with cheaper ones. Only the degree is considered; order is ignored even though it also affects cost.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 49.
