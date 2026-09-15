---
id: gncy.pso_adaptive_policy_rpo_adaptive_pso_config
label: rpo_adaptive_pso_config
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl
  symbol: rpo_adaptive_pso_config
  lines:
  - 25
  - 76
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: base_config
  type: RPOPSOConfig
  units: n/a
  required: true
  description: Baseline PSO configuration whose swarm size, iteration budget, and
    weights are rescaled by the geometry probe.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: adapted_config
  type: Tuple
  units: n/a
  description: Validated adapted configuration paired with a diagnostics tuple reporting
    distance, complexity, explore fraction, and whether adaptation ran.
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

# rpo_adaptive_pso_config

## Purpose
`rpo_adaptive_pso_config` sizes the PSO search to the difficulty of the specific transfer before any swarm is created. It runs a cheap straight-line probe between start and goal, converts what it finds into a scalar explore fraction, and rescales waypoint count, swarm size, iteration budget, and the PSO behavioural weights accordingly.

## Theory & Math
The explore fraction is $e = \mathrm{clamp}\left(\frac{w_c c + w_d \hat{d}}{w_c + w_d}, 0, 1\right)$ with complexity $c$ and normalised distance $\hat{d} = \mathrm{clamp}(d / d_{ref}, 0, 1)$. Effort scales as $f = f_{\min} + (f_{\max}-f_{\min})e$.

## Model & Assumptions
Difficulty is measured on the straight-line corridor only. `rpo_estimate_geometry_complexity` samples that corridor, computes clearance statistics, and blends a binary buffer-violation indicator at weight 0.7 with a reciprocal clearance term at weight 0.3, clamped to the unit interval. Distance is normalised against `cost_ref_distance_m`. The two are combined using `adaptive_complexity_weight` and `adaptive_distance_weight`, normalised by their sum with a 1e-9 floor, so the explore fraction stays in the unit interval regardless of how the weights are set.

## Design & Implementation
When `adaptive_enable` is false the base configuration is validated and returned with zeroed diagnostics. Otherwise every adapted quantity is clamped into an explicit band. The `adaptive_allow_downscale` flag decides whether the bands may go below the base values or are floored at them, which lets an operator guarantee that adaptation never weakens a hand-tuned configuration. Waypoint count grows linearly with complexity through `adaptive_waypoint_gain`. Particle count and iteration count are both multiplied by an effort scale interpolated between `adaptive_effort_min_fraction` and `adaptive_effort_max_fraction`. Length weight falls as complexity rises through the factor `1.25 - 0.5 * complexity`, inertia rises as `0.45 + 0.3 * explore`, the cognitive coefficient falls as `1.2 + 0.4 * (1 - explore)`, the social coefficient rises as `1.2 + 0.8 * explore`, and spread scale is multiplied by `0.75 + explore`. The result passes through `validate_rpo_pso_config`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `base_config` | RPOPSOConfig | n/a | yes | Baseline PSO configuration whose swarm size, iteration budget, and weights are rescaled by the geometry probe. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `adapted_config` | Tuple | n/a | — | Validated adapted configuration paired with a diagnostics tuple reporting distance, complexity, explore fraction, and whether adaptation ran. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:190-190`

**Downstream**

- `callees` → [[gnc.pso_adaptive_policy_rpo_estimate_geometry_complexity|rpo_estimate_geometry_complexity]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:31-31`
- `callees` → [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:63-63`
- `callees` → [[gnc.pso_parameters_validate_rpo_pso_config|validate_rpo_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_adaptive_policy.jl:28-28`
<!-- vulcan:connections:end -->

## Limitations
Probing only the straight-line corridor means an obstacle that blocks the required detour but not the direct line is invisible to the difficulty estimate, so a hard case can be sized as easy. `rpo_probe_geometry_metrics` returns a hard-coded `detour_ratio` of one, so that field carries no information. The magic constants in the weight formulas are not exposed as configuration, unlike the clamping bands around them.

## Provenance
Mapped from pso_adaptive_policy.jl lines 1-76; include site observed at guidance_hooks.jl line 69.
