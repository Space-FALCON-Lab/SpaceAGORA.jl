---
id: parallel.outer_route_selection__route_max_link_bucket
label: _route_max_link_bucket
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _route_max_link_bucket
  lines:
  - 27
  - 27
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
  description: Return value of `_route_max_link_bucket`.
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

# _route_max_link_bucket

## Purpose
Buckets the maximum per-spacecraft link fan-out (`max_links_per_sat`) into five string classes for the `maxlinks=` signature field, capturing how unevenly interaction work is distributed across the constellation.

## Design & Implementation
`@inline _route_max_link_bucket(n_links::Int)::String` returns "1" for `<= 1`, "2" for 2, "3_4" for 3 to 4, "5_8" for 5 to 8, and "9p" otherwise. It differs from `_route_link_bucket` by splitting out the "2" case. It appears in `outer_route_signature` and the mid-level signature of `_outer_route_signature_hierarchy` but not in the legacy `_compat_outer_route_signature`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `n_links` | Int | n/a | yes | Positional argument `n_links`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_route_max_link_bucket`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- [[parallel.outer_route_selection__outer_route_signature_hierarchy|_outer_route_signature_hierarchy]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:187-187`
- [[parallel.outer_route_selection_outer_route_signature|outer_route_signature]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:159-159`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the legacy signature omits this field, histories recorded under the old scheme cannot distinguish fan-out; the hierarchy fallback therefore mixes workloads with different `maxlinks` buckets. Non-positive inputs silently become "1".

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 27.
