---
id: parallel.outer_route_selection__outer_route_signature_hierarchy
label: _outer_route_signature_hierarchy
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: _outer_route_signature_hierarchy
  lines:
  - 181
  - 181
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
  type: Vector{String}
  units: n/a
  description: Return value of `_outer_route_signature_hierarchy`.
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

# _outer_route_signature_hierarchy

## Purpose
Returns an ordered list of progressively coarser signature keys (full, mid, legacy) so `select_outer_route!` can fall back to less specific routing history when the exact workload class has never been observed.

## Design & Implementation
`@inline _outer_route_signature_hierarchy(f::OuterRouteFeatures)::Vector{String}` computes `full = outer_route_signature(f)`, a `mid` signature of 14 fields (the full set minus the four rate buckets, the two GRAM flags, and the two effector counts), and `legacy = _compat_outer_route_signature(f)`, then returns `unique(String[full, mid, legacy])`. The caller iterates in order and takes the first signature with a non-empty stats snapshot.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | OuterRouteFeatures | n/a | yes | Positional argument `f`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{String} | n/a | — | Return value of `_outer_route_signature_hierarchy`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.select_outer_route_select_outer_route_bang|select_outer_route!]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:566-566`
- [[parcore.outer_route_metrics_record_outer_route_feedback_bang|record_outer_route_feedback!]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_metrics.jl:41-41`

**Downstream**

- `callees` → [[parallel.outer_route_selection__compat_outer_route_signature|_compat_outer_route_signature]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:199-199`
- `callees` → [[parallel.outer_route_selection__route_density_bucket|_route_density_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:194-194`
- `callees` → [[parallel.outer_route_selection__route_effector_cost_bucket|_route_effector_cost_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:197-197`
- `callees` → [[parallel.outer_route_selection__route_harmonics_bucket|_route_harmonics_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:191-191`
- `callees` → [[parallel.outer_route_selection__route_link_bucket|_route_link_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:186-186`
- `callees` → [[parallel.outer_route_selection__route_max_link_bucket|_route_max_link_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:187-187`
- `callees` → [[parallel.outer_route_selection__route_mission_bucket|_route_mission_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:188-188`
- `callees` → [[parallel.outer_route_selection__route_sat_bucket|_route_sat_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:185-185`
- `callees` → [[parallel.outer_route_selection__route_solver_bucket|_route_solver_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:195-195`
- `callees` → [[parallel.outer_route_selection_outer_route_signature|outer_route_signature]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:182-182`
<!-- vulcan:connections:end -->

## Limitations
The mid signature is duplicated literally rather than derived from the full one, so the two lists must be kept in sync by hand. `unique` only removes exact duplicates, which never happens because the field counts differ, so the call is effectively a no-op allocation. Feedback recorded is always written under the full key, so coarser keys only ever hold pre-schema history.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 181.
