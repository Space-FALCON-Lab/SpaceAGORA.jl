---
id: gnc.path_costs_rpo_normalized_path_cost_components
label: rpo_normalized_path_cost_components
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_costs.jl
  symbol: rpo_normalized_path_cost_components
  lines:
  - 86
  - 86
inputs:
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
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
- id: cost_cutoff
  type: Real
  units: n/a
  required: false
  description: Keyword argument `cost_cutoff` (default `Inf`).
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
  description: Return value of `rpo_normalized_path_cost_components`. Returns `(`.
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

# rpo_normalized_path_cost_components

## Purpose
Evaluates the full HYPR/RPO objective for one candidate path, returning every cost component (length, obstacle, fuel) in raw and normalized form together with clearance statistics, so that PSO, RRT-Connect, CHOMP and the post-refinement passes all score paths identically.

## Theory & Math
Total cost: $J = w_{len}\left(\frac{L}{L_{ref}}\right)^2 + w_{obs} J_{obs} + w_{fuel}\left(\frac{F}{F_{ref}}\right)^2$ where $L$ is the sampled path length (m), $L_{ref}$ the reference length (m), $J_{obs}=\sum_j \sigma\big(-k(d_j - (d_{safe}-\delta))\big)$ the summed sigmoid clearance penalty over samples with clearance $d_j$ (m), gain $k$ and tolerance $\delta$ (m), and $F$ the fuel proxy (kg) with $F_{ref} = m\,L_{ref}/(t_f\, I_{sp}\, g_0)$ for mass $m$ (kg), horizon $t_f$ (s), specific impulse $I_{sp}$ (s) and $g_0$ (m/s^2).

## Design & Implementation
Takes `points` (3xN control points), the station `geometry`, and an `RPOPSOConfig` `cfg`, with keyword `safe_distance_m` (default 0.0) and `cost_cutoff` (default `Inf`). It first resamples the path with `rpo_sample_path` at the density from `rpo_hypr_sampling_density_m`, then calls `rpo_clearance_stats_from_samples` which accumulates the sigmoid obstacle score `J_obs` using `cfg.obstacle_sigmoid_k` and `cfg.obstacle_sigmoid_tol_m`. Evaluation is staged for early exit: if the clearance pass reports `cutoff_exceeded`, a NamedTuple with `total=Inf` and zeroed length/fuel terms is returned; otherwise `rpo_path_cost_normalization_refs` supplies `len_ref` and `fuel_ref`, `J_len` is the sampled path length and the partial cost `w_obs*J_obs + w_len*J_len_norm^2` is checked against `cost_cutoff` before the comparatively expensive `rpo_fuel_proxy_from_samples` runs. The final `total` is `w_len*J_len_norm^2 + w_obs*J_obs + w_fuel*J_fuel_norm^2`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `safe_distance_m` | Real | n/a | no | Keyword argument `safe_distance_m` (default `0.0`). |
| in | `cost_cutoff` | Real | n/a | no | Keyword argument `cost_cutoff` (default `Inf`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_normalized_path_cost_components`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:614-614`
- [[gnc.pso_path_planning_evaluate_position|evaluate_position]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:370-370`
- [[gnc.pso_refinement_rpo_try_accept_refinement|rpo_try_accept_refinement]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:166-166`
- [[gnc.rrt_connect_rpo_rrt_connect_bezier_plan_path|rpo_rrt_connect_bezier_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:454-454`
- [[gnc.rrt_connect_rpo_rrt_star_plan_path|rpo_rrt_star_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:514-514`
- [[gnc.trajectory_optimizers_rpo_stomp_plan_path|rpo_stomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:392-392`
- [[gncy.path_costs_rpo_path_cost|rpo_path_cost]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:161-161`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:210-210`
- [[gncy.pso_refinement_rpo_post_refine_path|rpo_post_refine_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:275-275`
- [[gncy.rrt_connect_rpo_rrt_connect_plan_path|rpo_rrt_connect_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:334-334`
- [[gncy.trajectory_optimizers_rpo_chomp_plan_path|rpo_chomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:258-258`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:129-129`
- `callees` → [[gnc.path_costs_rpo_clearance_stats_from_samples|rpo_clearance_stats_from_samples]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:101-101`
- `callees` → [[gnc.path_costs_rpo_fuel_proxy_from_samples|rpo_fuel_proxy_from_samples]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:143-143`
- `callees` → [[gnc.path_costs_rpo_path_cost_normalization_refs|rpo_path_cost_normalization_refs]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:125-125`
- `callees` → [[gnc.path_sampling_rpo_path_length|rpo_path_length]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:126-126`
- `callees` → [[gnc.path_sampling_rpo_sample_path|rpo_sample_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:93-93`
- `callees` → [[gnc.pso_parameters_rpo_hypr_sampling_density_m|rpo_hypr_sampling_density_m]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_costs.jl:98-98`
<!-- vulcan:connections:end -->

## Limitations
Length and fuel terms are squared after normalization while the obstacle term is linear, so the relative weighting depends on `len_ref` and `fuel_ref` being sensible; when `cfg.cost_ref_distance_m <= 0` the reference falls back to the straight-line endpoint distance, clamped to at least `sample_ds_m` and 1e-6 m. Early-exit returns report `len_ref=0.0`/`fuel_ref=0.0` or `J_fuel=0.0`, so callers must not interpret those fields as physical values when `total == Inf`. The fuel proxy uses a second-difference acceleration estimate and assumes uniformly spaced samples in time over `cfg.tf_s`.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_costs.jl` line 86.
