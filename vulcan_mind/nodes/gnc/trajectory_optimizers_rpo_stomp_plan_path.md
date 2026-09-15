---
id: gnc.trajectory_optimizers_rpo_stomp_plan_path
label: rpo_stomp_plan_path
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: rpo_stomp_plan_path
  lines:
  - 356
  - 356
inputs:
- id: start_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `start_rtn`.
- id: goal_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `goal_rtn`.
- id: geometry
  type: Any
  units: n/a
  required: true
  description: Positional argument `geometry`.
- id: cfg
  type: RPOPSOConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: safe_distance_m
  type: Real
  units: n/a
  required: false
  description: Keyword argument `safe_distance_m` (default `0.0`).
- id: settings
  type: RPOSTOMPSettings
  units: n/a
  required: false
  description: Keyword argument `settings` (default `RPOSTOMPSettings(n_iters=cfg.n_iters)`).
- id: optimizer
  type: RPOTrajectoryOptimizerSettings
  units: n/a
  required: false
  description: Keyword argument `optimizer` (default `RPOTrajectoryOptimizerSettings()`).
- id: max_runtime_s
  type: Real
  units: n/a
  required: false
  description: Keyword argument `max_runtime_s` (default `Inf`).
- id: rng
  type: Any
  units: n/a
  required: false
  description: Keyword argument `rng` (default `Random.default_rng()`).
- id: initial_path
  type: Any
  units: n/a
  required: false
  description: Keyword argument `initial_path` (default `nothing`).
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
  description: Return value of `rpo_stomp_plan_path`. Returns `(`.
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

# rpo_stomp_plan_path

## Purpose
Runs the STOMP-like stochastic trajectory optimizer for RPO comparison: perturbs internal waypoints with second-difference-correlated Gaussian noise, weights rollouts by exponentiated local cost, applies a smoothed update with backtracking, and returns a HYPR-compatible result named tuple.

## Theory & Math
$$\epsilon_k \sim \mathcal{N}(0, \sigma^2 R^{-1}),\qquad w_{ik} = \frac{e^{-(S_{ik} - \min_k S_{ik})/\lambda}}{\sum_k e^{-(S_{ik} - \min_k S_{ik})/\lambda}},\qquad \delta_i = \sum_k w_{ik}\,\epsilon_{ik},\qquad \theta \leftarrow \theta + \eta\, M\delta$$ with $\sigma$ = `noise_std`, $\lambda$ = `lambda`, $\eta$ = `update_step` times the backtracking scale, and $S_{ik}$ the local cost of rollout $k$ at waypoint $i$.

## Design & Implementation
Signature `rpo_stomp_plan_path(start_rtn, goal_rtn, geometry, cfg::RPOPSOConfig; safe_distance_m=0.0, settings::RPOSTOMPSettings, optimizer::RPOTrajectoryOptimizerSettings, max_runtime_s=Inf, rng=Random.default_rng(), initial_path=nothing)`. Setup: `local_cfg = rpo_pso_config(cfg; w_obs=cfg.w_obs * optimizer.w_obs_scale)`, `n_internal = max(1, cfg.n_waypoints)`, seed via `rpo_trajectory_internal_points_from_seed`, bounds widened to contain the seed, `(R_inv, M)` from `rpo_second_difference_metric`, and `noise_chol = cholesky(Symmetric(noise_std^2 * R_inv + 1e-10 I)).L`. Each iteration draws `K = n_rollouts` perturbations `eps[d,:] = noise_chol * randn(rng, n)`, clamps, records the actual `eps`, evaluates the global objective and per-waypoint `rpo_stomp_waypoint_state_cost`; per waypoint it forms softmax weights `exp(-(cost - min)/lambda)` (uniform fallback when the denominator is `<= 1e-12` or non-finite), builds `delta`, smooths with `M`, then tries `update_step * scale` for scales (1, 0.5, 0.25, 0.125) accepting the first objective decrease; the best raw rollout replaces the candidate if better still. Accepted steps update `theta`, `best_points`, `best_components` (hard cost). Histories, stall counting, runtime checks, and `force_full_iters` mirror CHOMP. Finishes with `rpo_post_refine_path` and returns `(path, raw_path, cost, raw_cost, components, raw_components, config, adaptive=(enabled=false,), refinement_improved, cost_history, history, iterations, objective)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `start_rtn` | Any | n/a | yes | Positional argument `start_rtn`. |
| in | `goal_rtn` | Any | n/a | yes | Positional argument `goal_rtn`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `safe_distance_m` | Real | n/a | no | Keyword argument `safe_distance_m` (default `0.0`). |
| in | `settings` | RPOSTOMPSettings | n/a | no | Keyword argument `settings` (default `RPOSTOMPSettings(n_iters=cfg.n_iters)`). |
| in | `optimizer` | RPOTrajectoryOptimizerSettings | n/a | no | Keyword argument `optimizer` (default `RPOTrajectoryOptimizerSettings()`). |
| in | `max_runtime_s` | Real | n/a | no | Keyword argument `max_runtime_s` (default `Inf`). |
| in | `rng` | Any | n/a | no | Keyword argument `rng` (default `Random.default_rng()`). |
| in | `initial_path` | Any | n/a | no | Keyword argument `initial_path` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_stomp_plan_path`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_runtime_limited_iters|runtime_limited_iters]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:419-419`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:379-379`
- `callees` → [[gnc.path_costs_rpo_normalized_path_cost_components|rpo_normalized_path_cost_components]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:392-392`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:488-488`
- `callees` → [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:368-368`
- `callees` → [[gnc.trajectory_optimizers_rpo_clamp_internal_waypoints|rpo_clamp_internal_waypoints]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:418-418`
- `callees` → [[gnc.trajectory_optimizers_rpo_second_difference_metric|rpo_second_difference_metric]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:376-376`
- `callees` → [[gnc.trajectory_optimizers_rpo_stomp_waypoint_state_cost|rpo_stomp_waypoint_state_cost]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:425-425`
- `callees` → [[gnc.trajectory_optimizers_rpo_trajectory_internal_points_from_seed|rpo_trajectory_internal_points_from_seed]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:372-372`
- `callees` → [[gnc.trajectory_optimizers_rpo_trajectory_points_from_internal|rpo_trajectory_points_from_internal]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:390-390`
- `callees` → [[gnc.trajectory_optimizers_rpo_trajectory_search_bounds|rpo_trajectory_search_bounds]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:373-373`
- `callees` → [[gnc.trajectory_optimizers_rpo_trajectory_soft_objective|rpo_trajectory_soft_objective]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:381-381`
- `callees` → [[gnc.trajectory_optimizers_rpostompsettings|RPOSTOMPSettings]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:362-362`
- `callees` → [[gnc.trajectory_optimizers_rpotrajectoryoptimizersettings|RPOTrajectoryOptimizerSettings]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:363-363`
- `callees` → [[gncy.pso_refinement_rpo_post_refine_path|rpo_post_refine_path]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:502-502`
<!-- vulcan:connections:end -->

## Limitations
Per iteration it allocates `K` rollout matrices plus a `3 x n x K` array and performs `K + 4` global objective evaluations and `K * n` clearance queries. The Cholesky factorisation can throw `PosDefException` if `R_inv` is ill-conditioned despite the `1e-10` jitter. `cost_history` records the hard normalised cost while acceptance uses the soft objective, so the plotted curve is not guaranteed monotone. Results depend on `rng`; the default global RNG makes runs non-reproducible unless seeded by the caller.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl` line 356.
