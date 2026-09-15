---
id: parallel.outer_route_selection_outer_route_signature
label: outer_route_signature
kind: function
source:
  file: src/parallel/routing/outer_route_selection.jl
  symbol: outer_route_signature
  lines:
  - 154
  - 154
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
  type: String
  units: n/a
  description: Return value of `outer_route_signature`.
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

# outer_route_signature

## Purpose
Public function producing the canonical 22-field pipe-delimited routing signature that keys the adaptive feedback history (`OuterRouteState.history`) for a workload described by `OuterRouteFeatures`.

## Design & Implementation
`@inline outer_route_signature(f::OuterRouteFeatures)::String` joins, in fixed order: `cat`, `sat`, `links`, `maxlinks`, `mission`, `nbody`, `srp`, `harm`, `ctrl`, `orient`, `dens`, `solver`, `dt`, `ctrl_rate`, `guid_rate`, `nav_rate`, `gram_srg`, `gram_grid`, `ctrl_eff`, `thermal`, `eff_cnt`, and `eff_cost`. Each numeric field is passed through the corresponding `_route_*_bucket` function and each Bool becomes "1" or "0". The result is the first (most specific) entry of `_outer_route_signature_hierarchy` and the key `select_outer_route!` reports in its trace line.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | OuterRouteFeatures | n/a | yes | Positional argument `f`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `outer_route_signature`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/routing/outer_route_selection.jl`
- [[parallel.outer_route_selection__compat_outer_route_signature|_compat_outer_route_signature]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:149-149`
- [[parallel.outer_route_selection__outer_route_signature_hierarchy|_outer_route_signature_hierarchy]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_selection.jl:182-182`

**Downstream**

- `callees` → [[parallel.outer_route_selection__route_count_bucket|_route_count_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:174-174`
- `callees` → [[parallel.outer_route_selection__route_density_bucket|_route_density_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:166-166`
- `callees` → [[parallel.outer_route_selection__route_effector_cost_bucket|_route_effector_cost_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:177-177`
- `callees` → [[parallel.outer_route_selection__route_harmonics_bucket|_route_harmonics_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:163-163`
- `callees` → [[parallel.outer_route_selection__route_interval_bucket|_route_interval_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:168-168`
- `callees` → [[parallel.outer_route_selection__route_link_bucket|_route_link_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:158-158`
- `callees` → [[parallel.outer_route_selection__route_max_link_bucket|_route_max_link_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:159-159`
- `callees` → [[parallel.outer_route_selection__route_mission_bucket|_route_mission_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:160-160`
- `callees` → [[parallel.outer_route_selection__route_sat_bucket|_route_sat_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:157-157`
- `callees` → [[parallel.outer_route_selection__route_solver_bucket|_route_solver_bucket]] · `callers` · call · `src/parallel/routing/outer_route_selection.jl:167-167`
<!-- vulcan:connections:end -->

## Limitations
Field order is part of the key format; reordering or adding a field invalidates all persisted history. `f.category` is not sanitised for `|` or `=`. Every call allocates 22 interpolated strings plus the join, which is fine per simulation launch but not per time step.

## Provenance
Mapped from `src/parallel/routing/outer_route_selection.jl` line 154.
