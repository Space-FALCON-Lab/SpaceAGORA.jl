---
id: parallel.outer_route_selection__route_link_bucket
label: _route_link_bucket
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _route_link_bucket
  lines:
  - 16
  - 16
inputs:
- id: n_links
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_links`.
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
  description: Return value of `_route_link_bucket`.
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

# _route_link_bucket

## Purpose
Buckets the total number of inter-spacecraft interaction links (`n_links`) into four string classes for inclusion in routing signatures, so link-count variations that do not change parallel cost do not fragment the feedback history.

## Design & Implementation
`@inline _route_link_bucket(n_links::Int)::String` returns "1" for `<= 1`, "2_4" for 2 to 4, "5_8" for 5 to 8, and "9p" above 8. Consumed by all three signature builders in this file under the `links=` key, and the same thresholds are what `_feature_is_lightweight` compares indirectly through `t.outer_light_link_threshold`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `n_links` | Int | n/a | yes | Positional argument `n_links`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_route_link_bucket`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- [[parallel.outer_route_selection__compat_outer_route_signature|_compat_outer_route_signature]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:138-138`
- [[parallel.outer_route_selection__outer_route_signature_hierarchy|_outer_route_signature_hierarchy]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:186-186`
- [[parallel.outer_route_selection_outer_route_signature|outer_route_signature]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:158-158`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Negative or zero link counts map to "1" rather than a dedicated "0" bucket, so a constellation with no links is indistinguishable from one with a single link. The bucket boundaries are constants and cannot be changed through `OuterRouteTuning`.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 16.
