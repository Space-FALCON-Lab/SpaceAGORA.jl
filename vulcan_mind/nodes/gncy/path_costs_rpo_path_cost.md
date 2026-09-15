---
id: gncy.path_costs_rpo_path_cost
label: rpo_path_cost
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_costs.jl
  symbol: rpo_path_cost
  lines:
  - 160
  - 162
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: candidate_path
  type: Matrix
  units: m
  required: true
  description: Three-row matrix of RTN control points or waypoints describing the
    candidate path to be scored.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: total_cost
  type: Float64
  units: n/a
  description: Scalar weighted objective combining normalised squared path length,
    obstacle penalty, and normalised squared fuel proxy.
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

# rpo_path_cost

## Purpose
`rpo_path_cost` is the scalar objective every RPO planner minimises. It reduces the full named-tuple cost breakdown produced by `rpo_normalized_path_cost_components` to a single number so PSO, CHOMP, STOMP, and the RRT refinement stage can all be compared on one scale.

## Theory & Math
The objective is $J = w_{len}\left(\frac{L}{L_{ref}}\right)^2 + w_{obs} J_{obs} + w_{fuel}\left(\frac{F}{F_{ref}}\right)^2$, with $L$ the sampled arc length, $F$ the fuel proxy, and $J_{obs}$ the sigmoid clearance penalty.

## Model & Assumptions
The objective normalises length and fuel against reference values from `rpo_path_cost_normalization_refs` so the three terms are dimensionless and the weights `w_len`, `w_obs`, and `w_fuel` are directly comparable. Length and fuel enter squared, which penalises large excursions superlinearly, while the obstacle score enters linearly because it is already a saturating sigmoid of clearance controlled by `obstacle_sigmoid_k` and `obstacle_sigmoid_tol_m`. The default obstacle weight of one million makes any keep-out intrusion dominate the objective, which is how the planner enforces safety through a soft penalty.

## Design & Implementation
The wrapper forwards `safe_distance_m` and `cost_cutoff` straight through and takes the `total` field. The components function first samples the path with `rpo_sample_path` at the density returned by `rpo_hypr_sampling_density_m`, then computes clearance statistics. Two early exits keep the hot path cheap. If the clearance pass reports `cutoff_exceeded` the function returns infinity immediately without measuring length or fuel. Otherwise a partial cost of the obstacle and length terms is checked against the cutoff, and only if it passes is the fuel proxy computed. Both early exits still report minimum clearance and violation count so a rejected candidate remains diagnosable.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `candidate_path` | Matrix | m | yes | Three-row matrix of RTN control points or waypoints describing the candidate path to be scored. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `total_cost` | Float64 | n/a | — | Scalar weighted objective combining normalised squared path length, obstacle penalty, and normalised squared fuel proxy. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.pso_refinement_rpo_post_refine_path|rpo_post_refine_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:277-277`

**Downstream**

- `callees` → [[gnc.path_costs_rpo_normalized_path_cost_components|rpo_normalized_path_cost_components]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:161-161`
<!-- vulcan:connections:end -->

## Limitations
Because obstacle avoidance is a weighted penalty rather than a hard constraint, a sufficiently long detour can score worse than a marginal keep-out intrusion when weights are mis-tuned. The returned scalar discards the component breakdown, so callers that need to know why a path scored badly must call the components function directly. Cutoff-driven early exit means costs above the cutoff are not comparable to each other, which makes the scalar unsuitable for ranking rejected candidates.

## Provenance
Mapped from path_costs.jl lines 86-162; include site observed at guidance_hooks.jl line 68.
