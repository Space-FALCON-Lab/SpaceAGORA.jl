---
id: gncy.pso_helpers_rpo_pso_bounds
label: rpo_pso_bounds
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_helpers.jl
  symbol: rpo_pso_bounds
  lines:
  - 2
  - 11
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: endpoints
  type: Tuple
  units: m
  required: true
  description: Start and goal RTN positions together with the PSO configuration supplying
    search margin and spread scale.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: search_box
  type: Tuple
  units: m
  description: Lower and upper corner vectors of the axis-aligned box within which
    PSO waypoints are sampled and clamped.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncy
origin: agent
---

# rpo_pso_bounds

## Purpose
`rpo_pso_bounds` builds the axis-aligned search box in which PSO interior waypoints live. Every particle position is initialised inside this box and clamped to it after each velocity update, so the box is what bounds the reachable set of the swarm.

## Theory & Math
With endpoints $p_s, p_g$, margin $m$, spread $\sigma$ and span $L = \max(\lVert p_g - p_s\rVert, m)$, the bounds are $lo = \min\left(\min(p_s,p_g) - m,\ \tfrac{1}{2}(p_s+p_g) - \sigma L\right)$ and $hi = \max\left(\max(p_s,p_g) + m,\ \tfrac{1}{2}(p_s+p_g) + \sigma L\right)$, elementwise.

## Model & Assumptions
The box is constructed to satisfy two requirements at once. It must contain both endpoints with at least `search_margin_m` of slack, which is the minimum needed for a direct path plus small deviations. It must also allow lateral excursions proportional to the transfer size, because a detour around a station module has to bulge away from the straight line by an amount that scales with the distance travelled, not with a fixed margin. The span used for that second requirement is the larger of the endpoint separation and the search margin, so the box never collapses for coincident endpoints.

## Design & Implementation
Both endpoints are converted to static three-vectors. The first box is the elementwise minimum and maximum of the endpoints, expanded by the search margin in every direction. The second box is centred on the endpoint midpoint and extends `spread_scale * span` in each direction. The final bounds take the elementwise minimum of the two lower corners and the elementwise maximum of the two upper corners, so the result is the union bounding box rather than an intersection. The companion `rpo_pso_warmstart_bounds` instead brackets an RRT warm-start path with `rrt_warmstart_box_margin_m` and validates that the path has three rows and at least two columns.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `endpoints` | Tuple | m | yes | Start and goal RTN positions together with the PSO configuration supplying search margin and spread scale. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `search_box` | Tuple | m | — | Lower and upper corner vectors of the axis-aligned box within which PSO waypoints are sampled and clamped. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_reset_swarm_bang|reset_swarm!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:327-327`
- [[gnc.pso_refinement_rpo_refinement_clamp_path|rpo_refinement_clamp_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:19-19`
- [[gnc.rrt_connect_rpo_rrt_star_plan_path|rpo_rrt_star_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:501-501`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:327-327`
- [[gncy.rrt_connect_rpo_rrt_connect_plan_path|rpo_rrt_connect_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:321-321`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
An axis-aligned box in the RTN frame is a coarse fit to a transfer that runs diagonally, so a large fraction of the sampled volume can be irrelevant and the swarm wastes evaluations there. The box ignores geometry entirely, so it can span straight through a station module. Because bounds are symmetric about the midpoint, an asymmetric obstacle field is not reflected in where the swarm is allowed to search.

## Provenance
Mapped from pso_helpers.jl lines 1-35; include site observed at guidance_hooks.jl line 70.
