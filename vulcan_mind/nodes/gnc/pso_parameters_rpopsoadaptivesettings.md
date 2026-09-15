---
id: gnc.pso_parameters_rpopsoadaptivesettings
label: RPOPSOAdaptiveSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: RPOPSOAdaptiveSettings
  lines:
  - 30
  - 30
inputs:
- id: enabled
  type: Bool
  units: n/a
  required: false
  description: Field `enabled` (default `true`).
- id: allow_downscale
  type: Bool
  units: n/a
  required: false
  description: Field `allow_downscale` (default `false`).
- id: complexity_weight
  type: Float64
  units: n/a
  required: false
  description: Field `complexity_weight` (default `0.35`).
- id: distance_weight
  type: Float64
  units: n/a
  required: false
  description: Field `distance_weight` (default `0.65`).
- id: waypoint_gain
  type: Float64
  units: n/a
  required: false
  description: Field `waypoint_gain` (default `3.0`).
- id: effort_min_fraction
  type: Float64
  units: n/a
  required: false
  description: Field `effort_min_fraction` (default `0.75`).
- id: effort_max_fraction
  type: Float64
  units: n/a
  required: false
  description: Field `effort_max_fraction` (default `1.5`).
- id: n_waypoints_min
  type: Int
  units: n/a
  required: false
  description: Field `n_waypoints_min` (default `3`).
- id: n_waypoints_max
  type: Int
  units: n/a
  required: false
  description: Field `n_waypoints_max` (default `8`).
- id: n_particles_min
  type: Int
  units: n/a
  required: false
  description: Field `n_particles_min` (default `140`).
- id: n_particles_max
  type: Int
  units: n/a
  required: false
  description: Field `n_particles_max` (default `320`).
- id: n_iters_min
  type: Int
  units: n/a
  required: false
  description: Field `n_iters_min` (default `14`).
- id: n_iters_max
  type: Int
  units: n/a
  required: false
  description: Field `n_iters_max` (default `55`).
- id: w_len_min
  type: Float64
  units: n/a
  required: false
  description: Field `w_len_min` (default `0.5`).
- id: w_len_max
  type: Float64
  units: n/a
  required: false
  description: Field `w_len_max` (default `1.5`).
- id: w_obs_min
  type: Float64
  units: n/a
  required: false
  description: Field `w_obs_min` (default `1.0e6`).
- id: w_obs_max
  type: Float64
  units: n/a
  required: false
  description: Field `w_obs_max` (default `1.0e6`).
- id: w_inertia_min
  type: Float64
  units: n/a
  required: false
  description: Field `w_inertia_min` (default `0.4`).
- id: w_inertia_max
  type: Float64
  units: n/a
  required: false
  description: Field `w_inertia_max` (default `0.75`).
- id: c1_min
  type: Float64
  units: n/a
  required: false
  description: Field `c1_min` (default `1.2`).
- id: c1_max
  type: Float64
  units: n/a
  required: false
  description: Field `c1_max` (default `1.8`).
- id: c2_min
  type: Float64
  units: n/a
  required: false
  description: Field `c2_min` (default `1.2`).
- id: c2_max
  type: Float64
  units: n/a
  required: false
  description: Field `c2_max` (default `2.2`).
- id: spread_scale_min
  type: Float64
  units: n/a
  required: false
  description: Field `spread_scale_min` (default `0.05`).
- id: spread_scale_max
  type: Float64
  units: n/a
  required: false
  description: Field `spread_scale_max` (default `0.5`).
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
  type: RPOPSOAdaptiveSettings
  units: n/a
  description: Constructed `RPOPSOAdaptiveSettings` (keyword constructor via @kwdef).
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

# RPOPSOAdaptiveSettings

## Purpose
Grouped struct for the adaptive effort-scaling stage of HYPR, which resizes the swarm (waypoints, particles, iterations) and tunes weights and PSO coefficients based on scene complexity and start-goal distance before the main search begins.

## Design & Implementation
Switches: `enabled::Bool = true`, `allow_downscale::Bool = false`. Effort model: `complexity_weight = 0.35`, `distance_weight = 0.65`, `waypoint_gain = 3.0`, `effort_min_fraction = 0.75`, `effort_max_fraction = 1.5`. Clamp ranges for each adapted quantity: `n_waypoints` in [3, 8], `n_particles` in [140, 320], `n_iters` in [14, 55], `w_len` in [0.5, 1.5], `w_obs` fixed at [1e6, 1e6], `w_inertia` in [0.4, 0.75], `c1` in [1.2, 1.8], `c2` in [1.2, 2.2], `spread_scale` in [0.05, 0.5]. Every field is copied to an `adaptive_*` field of `RPOPSOConfig`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `enabled` | Bool | n/a | no | Field `enabled` (default `true`). |
| in | `allow_downscale` | Bool | n/a | no | Field `allow_downscale` (default `false`). |
| in | `complexity_weight` | Float64 | n/a | no | Field `complexity_weight` (default `0.35`). |
| in | `distance_weight` | Float64 | n/a | no | Field `distance_weight` (default `0.65`). |
| in | `waypoint_gain` | Float64 | n/a | no | Field `waypoint_gain` (default `3.0`). |
| in | `effort_min_fraction` | Float64 | n/a | no | Field `effort_min_fraction` (default `0.75`). |
| in | `effort_max_fraction` | Float64 | n/a | no | Field `effort_max_fraction` (default `1.5`). |
| in | `n_waypoints_min` | Int | n/a | no | Field `n_waypoints_min` (default `3`). |
| in | `n_waypoints_max` | Int | n/a | no | Field `n_waypoints_max` (default `8`). |
| in | `n_particles_min` | Int | n/a | no | Field `n_particles_min` (default `140`). |
| in | `n_particles_max` | Int | n/a | no | Field `n_particles_max` (default `320`). |
| in | `n_iters_min` | Int | n/a | no | Field `n_iters_min` (default `14`). |
| in | `n_iters_max` | Int | n/a | no | Field `n_iters_max` (default `55`). |
| in | `w_len_min` | Float64 | n/a | no | Field `w_len_min` (default `0.5`). |
| in | `w_len_max` | Float64 | n/a | no | Field `w_len_max` (default `1.5`). |
| in | `w_obs_min` | Float64 | n/a | no | Field `w_obs_min` (default `1.0e6`). |
| in | `w_obs_max` | Float64 | n/a | no | Field `w_obs_max` (default `1.0e6`). |
| in | `w_inertia_min` | Float64 | n/a | no | Field `w_inertia_min` (default `0.4`). |
| in | `w_inertia_max` | Float64 | n/a | no | Field `w_inertia_max` (default `0.75`). |
| in | `c1_min` | Float64 | n/a | no | Field `c1_min` (default `1.2`). |
| in | `c1_max` | Float64 | n/a | no | Field `c1_max` (default `1.8`). |
| in | `c2_min` | Float64 | n/a | no | Field `c2_min` (default `1.2`). |
| in | `c2_max` | Float64 | n/a | no | Field `c2_max` (default `2.2`). |
| in | `spread_scale_min` | Float64 | n/a | no | Field `spread_scale_min` (default `0.05`). |
| in | `spread_scale_max` | Float64 | n/a | no | Field `spread_scale_max` (default `0.5`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPSOAdaptiveSettings | n/a | — | Constructed `RPOPSOAdaptiveSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_parameters_rpopsoconfigurator|RPOPSOConfigurator]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:172-172`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Validation happens downstream: each `*_min <= *_max` pair, non-negative weights, `effort_min_fraction > 0`, and `effort_max_fraction >= effort_min_fraction` are enforced by `validate_rpo_pso_config`, not here. `complexity_weight + distance_weight` is not required to sum to 1. The identical `w_obs_min`/`w_obs_max` defaults mean adaptation never changes the obstacle weight unless a caller widens the range.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 30.
