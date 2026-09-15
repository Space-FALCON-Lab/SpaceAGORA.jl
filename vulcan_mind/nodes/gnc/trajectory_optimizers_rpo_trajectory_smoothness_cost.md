---
id: gnc.trajectory_optimizers_rpo_trajectory_smoothness_cost
label: rpo_trajectory_smoothness_cost
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: rpo_trajectory_smoothness_cost
  lines:
  - 145
  - 145
inputs:
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
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
  type: Any
  units: n/a
  description: Return value of `rpo_trajectory_smoothness_cost`. Returns `acc`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# rpo_trajectory_smoothness_cost

## Purpose
Computes the raw smoothness penalty of a full waypoint trajectory as the sum of squared discrete second differences, the un-normalised acceleration proxy that `rpo_trajectory_soft_objective` later divides by the squared reference length.

## Theory & Math
$$J_{\text{smooth}} = \sum_{i=2}^{N-1} \left\| p_{i-1} - 2p_i + p_{i+1} \right\|^2$$ where $p_i \in \mathbb{R}^3$ is the $i$-th waypoint (RTN metres) and $N$ is the total number of waypoints including endpoints.

## Design & Implementation
Converts `points` to `Matrix{Float64}`, returns `0.0` when fewer than 3 columns exist, and otherwise loops `i` from 2 to `end-1` accumulating `dot(d2, d2)` with `d2 = pts[:, i-1] - 2 pts[:, i] + pts[:, i+1]`. Units are metres squared. Operates on the full trajectory including `start` and `goal`, so the first and last internal waypoints are penalised against the fixed endpoints.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_trajectory_smoothness_cost`. Returns `acc`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.trajectory_optimizers_rpo_trajectory_soft_objective|rpo_trajectory_soft_objective]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:200-200`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Each loop iteration allocates two column slices and a temporary vector; a view-based implementation would avoid that. The penalty is not scaled by waypoint spacing, so trajectories with more waypoints for the same geometry report smaller per-segment differences and a different total, making cross-resolution comparisons non-trivial.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl` line 145.
