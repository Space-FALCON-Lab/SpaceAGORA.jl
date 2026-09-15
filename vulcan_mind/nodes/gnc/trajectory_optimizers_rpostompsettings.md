---
id: gnc.trajectory_optimizers_rpostompsettings
label: RPOSTOMPSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: RPOSTOMPSettings
  lines:
  - 12
  - 12
inputs:
- id: n_iters
  type: Int
  units: n/a
  required: false
  description: Field `n_iters` (default `100`).
- id: n_rollouts
  type: Int
  units: n/a
  required: false
  description: Field `n_rollouts` (default `20`).
- id: noise_std
  type: Float64
  units: n/a
  required: false
  description: Field `noise_std` (default `0.25`).
- id: lambda
  type: Float64
  units: n/a
  required: false
  description: Field `lambda` (default `10.0`).
- id: update_step
  type: Float64
  units: n/a
  required: false
  description: Field `update_step` (default `1.0`).
- id: w_smooth
  type: Float64
  units: n/a
  required: false
  description: Field `w_smooth` (default `1.0`).
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
  type: RPOSTOMPSettings
  units: n/a
  description: Constructed `RPOSTOMPSettings` (keyword constructor via @kwdef).
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

# RPOSTOMPSettings

## Purpose
Immutable `Base.@kwdef` settings struct for the STOMP-like stochastic trajectory optimizer in the RPO comparison suite. It parameterises the number of noisy rollouts per iteration, the exploration noise, the soft-min temperature, and the update step consumed by `rpo_stomp_plan_path`.

## Design & Implementation
Fields: `n_iters::Int = 100` outer iterations; `n_rollouts::Int = 20` (K) perturbed trajectories sampled each iteration; `noise_std::Float64 = 0.25` metres, squared and multiplied into the inverse second-difference metric `R_inv` to build the rollout covariance whose Cholesky factor colours `randn` noise; `lambda::Float64 = 10.0` is the temperature in `exp(-(cost - min_cost)/lambda)` for per-waypoint rollout weighting; `update_step::Float64 = 1.0` scales the smoothed weighted update before backtracking over scales (1, 0.5, 0.25, 0.125); `w_smooth::Float64 = 1.0` weights the second-difference term in both the global objective and `rpo_stomp_waypoint_state_cost`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `n_iters` | Int | n/a | no | Field `n_iters` (default `100`). |
| in | `n_rollouts` | Int | n/a | no | Field `n_rollouts` (default `20`). |
| in | `noise_std` | Float64 | n/a | no | Field `noise_std` (default `0.25`). |
| in | `lambda` | Float64 | n/a | no | Field `lambda` (default `10.0`). |
| in | `update_step` | Float64 | n/a | no | Field `update_step` (default `1.0`). |
| in | `w_smooth` | Float64 | n/a | no | Field `w_smooth` (default `1.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOSTOMPSettings | n/a | — | Constructed `RPOSTOMPSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpoplannercomparisonconfig|RPOPlannerComparisonConfig]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:36-36`
- [[gnc.planner_comparison_runtime_limited_iters|runtime_limited_iters]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:401-401`
- [[gnc.trajectory_optimizers_rpo_stomp_plan_path|rpo_stomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:362-362`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The struct performs no validation; the planner clamps `n_rollouts` to at least 1, floors `noise_std` at 1.0e-9 and `lambda` at 1.0e-9, but does not warn. `noise_std` is isotropic across the three RTN axes and constant over the run, so no annealing of exploration is possible through this struct.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl` line 12.
