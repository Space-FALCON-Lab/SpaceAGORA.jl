---
id: gnc.trajectory_optimizers_rpo_trajectory_soft_objective
label: rpo_trajectory_soft_objective
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: rpo_trajectory_soft_objective
  lines:
  - 182
  - 182
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
  type: Any
  units: n/a
  required: true
  description: Keyword argument `safe_distance_m`.
- id: obstacle_margin_m
  type: Any
  units: n/a
  required: true
  description: Keyword argument `obstacle_margin_m`.
- id: w_smooth
  type: Any
  units: n/a
  required: true
  description: Keyword argument `w_smooth`.
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
  description: Return value of `rpo_trajectory_soft_objective`. Returns `cfg.w_len
    * (J_len / refs.len_ref)^2 +`.
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

# rpo_trajectory_soft_objective

## Purpose
The scalar objective minimised by both CHOMP and STOMP: a normalised path-length term, fuel-proxy term, mean soft obstacle potential, and a length-normalised smoothness term, weighted by the `RPOPSOConfig` weights so results are comparable to HYPR's cost decomposition.

## Theory & Math
$$J = w_{len}\left(\frac{J_{len}}{L_{ref}}\right)^2 + w_{fuel}\left(\frac{J_{fuel}}{F_{ref}}\right)^2 + w_{obs}\,J_{obs}^{soft} + w_{smooth}\,\frac{J_{smooth}}{\max(L_{ref}^2, 10^{-9})}$$ where $L_{ref}$, $F_{ref}$ are the normalisation references from `rpo_path_cost_normalization_refs`, $J_{obs}^{soft}$ is the mean CHOMP potential, and the weights $w_\ast$ come from `cfg`.

## Design & Implementation
Signature `rpo_trajectory_soft_objective(points, geometry, cfg::RPOPSOConfig; safe_distance_m, obstacle_margin_m, w_smooth)`. It densifies `points` with `rpo_sample_path` using `cfg.sample_ds_m` and `cfg.curve_type`, obtains reference scales `refs` from `rpo_path_cost_normalization_refs`, then computes `J_len = rpo_path_length(samples)`, `J_fuel` (only when `cfg.w_fuel > 0`), `J_obs_soft` via `rpo_soft_obstacle_cost_from_samples`, and `J_smooth = rpo_trajectory_smoothness_cost(points) / max(refs.len_ref^2, 1e-9)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `safe_distance_m` | Any | n/a | yes | Keyword argument `safe_distance_m`. |
| in | `obstacle_margin_m` | Any | n/a | yes | Keyword argument `obstacle_margin_m`. |
| in | `w_smooth` | Any | n/a | yes | Keyword argument `w_smooth`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_trajectory_soft_objective`. Returns `cfg.w_len * (J_len / refs.len_ref)^2 +`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.trajectory_optimizers_rpo_stomp_plan_path|rpo_stomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:381-381`
- [[gncy.trajectory_optimizers_rpo_chomp_plan_path|rpo_chomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:247-247`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:204-204`
- `callees` → [[gnc.path_costs_rpo_fuel_proxy_from_samples|rpo_fuel_proxy_from_samples]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:193-193`
- `callees` → [[gnc.path_costs_rpo_path_cost_normalization_refs|rpo_path_cost_normalization_refs]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:191-191`
- `callees` → [[gnc.path_sampling_rpo_path_length|rpo_path_length]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:192-192`
- `callees` → [[gnc.path_sampling_rpo_sample_path|rpo_sample_path]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:183-183`
- `callees` → [[gnc.trajectory_optimizers_rpo_soft_obstacle_cost_from_samples|rpo_soft_obstacle_cost_from_samples]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:194-194`
- `callees` → [[gnc.trajectory_optimizers_rpo_trajectory_smoothness_cost|rpo_trajectory_smoothness_cost]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:200-200`
<!-- vulcan:connections:end -->

## Limitations
Obstacle cost is soft, so a converged optimum may still violate `safe_distance_m`; callers rely on `rpo_post_refine_path` and the hard-cost `rpo_normalized_path_cost_components` for reporting. Every evaluation re-samples and re-queries clearance, with no memoisation. If `refs.fuel_ref` is zero and `w_fuel > 0` the fuel term divides by zero.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl` line 182.
