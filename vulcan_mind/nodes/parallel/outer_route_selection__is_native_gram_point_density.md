---
id: parallel.outer_route_selection__is_native_gram_point_density
label: _is_native_gram_point_density
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _is_native_gram_point_density
  lines:
  - 268
  - 268
inputs:
- id: f
  type: OuterRouteFeatures
  units: n/a
  required: true
  description: Positional argument `f`.
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
  description: Return value of `_is_native_gram_point_density`.
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

# _is_native_gram_point_density

## Purpose
Detects workloads that call the native GRAM atmosphere library point-by-point without a surrogate or static grid, a case that is serialised by an internal lock and therefore only scales through process isolation.

## Design & Implementation
`@inline _is_native_gram_point_density(f::OuterRouteFeatures)::Bool` returns true when `_route_density_bucket(f.density_family) == "gram_pt"` and both `f.gram_surrogate_enabled` and `f.gram_static_grid_enabled` are false. `default_outer_route` short-circuits to `:process` and `outer_route_candidates` restricts candidates to `[:none, :process]` whenever this holds, and `_feature_heavy_for_process` returns true.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | OuterRouteFeatures | n/a | yes | Positional argument `f`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_is_native_gram_point_density`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.default_outer_route|default_outer_route]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:323-323`
- [[parallel.outer_route_selection__feature_heavy_for_process|_feature_heavy_for_process]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:275-275`
- [[parallel.outer_route_selection_outer_route_candidates|outer_route_candidates]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:381-381`

**Downstream**

- `callees` → [[parallel.outer_route_selection__route_density_bucket|_route_density_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:269-269`
<!-- vulcan:connections:end -->

## Limitations
Relies on the string bucket, so an unrecognised GRAM family alias (anything other than "gram_point" or "gram") is not detected and the workload may be routed to threads, where the GRAM lock would serialise it. It does not check whether the GRAM binary is actually available on the machine.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 268.
