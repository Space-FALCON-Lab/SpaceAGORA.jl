---
id: parallel.outer_route_selection__route_count_bucket
label: _route_count_bucket
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _route_count_bucket
  lines:
  - 111
  - 111
inputs:
- id: v
  type: Int
  units: n/a
  required: true
  description: Positional argument `v`.
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
  description: Return value of `_route_count_bucket`.
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

# _route_count_bucket

## Purpose
Buckets a small non-negative count, used for both `control_effector_count` and `dynamic_effector_count`, into five tokens for the `ctrl_eff=` and `eff_cnt=` signature fields.

## Design & Implementation
`@inline _route_count_bucket(v::Int)::String` returns "0" for `v <= 0`, "1" for 1, "2_3" for 2 or 3, "4_6" for 4 to 6, and "7p" above 6. Unlike `_route_link_bucket` it has an explicit zero bucket because the absence of effectors is a meaningful routing feature. Only the full signature includes these fields.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `v` | Int | n/a | yes | Positional argument `v`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_route_count_bucket`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- [[parallel.outer_route_selection_outer_route_signature|outer_route_signature]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:174-174`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Negative counts are treated as zero silently. Seven effectors and seventy share a bucket, so very effector-heavy spacecraft pool routing history with moderately equipped ones; the `eff_cost=` class partly compensates.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 111.
