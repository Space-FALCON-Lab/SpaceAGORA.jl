---
id: gnc.pso_parameters_rpopsorefinementsettings
label: RPOPSORefinementSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: RPOPSORefinementSettings
  lines:
  - 141
  - 141
inputs:
- id: enabled
  type: Bool
  units: n/a
  required: false
  description: Field `enabled` (default `true`).
- id: start_iter
  type: Int
  units: n/a
  required: false
  description: Field `start_iter` (default `60`).
- id: period
  type: Int
  units: n/a
  required: false
  description: Field `period` (default `10`).
- id: sample_ds_m
  type: Float64
  units: n/a
  required: false
  description: Field `sample_ds_m` (default `0.025`).
- id: merge_distance_m
  type: Float64
  units: n/a
  required: false
  description: Field `merge_distance_m` (default `1.0`).
- id: waypoint_passes
  type: Int
  units: n/a
  required: false
  description: Field `waypoint_passes` (default `24`).
- id: rounds
  type: Int
  units: n/a
  required: false
  description: Field `rounds` (default `16`).
- id: min_abs_cost_improvement
  type: Float64
  units: n/a
  required: false
  description: Field `min_abs_cost_improvement` (default `1.0e-8`).
- id: min_rel_cost_improvement
  type: Float64
  units: n/a
  required: false
  description: Field `min_rel_cost_improvement` (default `1.0e-5`).
- id: insert_straight_waypoints
  type: Bool
  units: n/a
  required: false
  description: Field `insert_straight_waypoints` (default `true`).
- id: straight_max_segment_length_m
  type: Float64
  units: n/a
  required: false
  description: Field `straight_max_segment_length_m` (default `2.0`).
- id: straight_max_inserted
  type: Int
  units: n/a
  required: false
  description: Field `straight_max_inserted` (default `12`).
- id: straight_clearance_margin_m
  type: Float64
  units: n/a
  required: false
  description: Field `straight_clearance_margin_m` (default `0.0`).
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
  type: RPOPSORefinementSettings
  units: n/a
  description: Constructed `RPOPSORefinementSettings` (keyword constructor via @kwdef).
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

# RPOPSORefinementSettings

## Purpose
Grouped struct for the post-PSO refinement stage: periodic in-loop shortcutting, waypoint merging and tightening, straight-segment waypoint insertion, and Bezier refit that polish the best path found by the swarm.

## Design & Implementation
Fields: `enabled = true`; `start_iter = 60` and `period = 10` gate in-loop refinement (with the default `n_iters = 55`, in-loop refinement never fires and only the final pass runs); `sample_ds_m = 0.025` fine collision spacing; `merge_distance_m = 1.0` for collapsing nearby waypoints; `waypoint_passes = 24` and `rounds = 16` iteration counts; `min_abs_cost_improvement = 1e-8` and `min_rel_cost_improvement = 1e-5` convergence thresholds; `insert_straight_waypoints = true`, `straight_max_segment_length_m = 2.0`, `straight_max_inserted = 12`, `straight_clearance_margin_m = 0.0` for straight-segment densification. Mapped to `refinement_*` in `RPOPSOConfig`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `enabled` | Bool | n/a | no | Field `enabled` (default `true`). |
| in | `start_iter` | Int | n/a | no | Field `start_iter` (default `60`). |
| in | `period` | Int | n/a | no | Field `period` (default `10`). |
| in | `sample_ds_m` | Float64 | n/a | no | Field `sample_ds_m` (default `0.025`). |
| in | `merge_distance_m` | Float64 | n/a | no | Field `merge_distance_m` (default `1.0`). |
| in | `waypoint_passes` | Int | n/a | no | Field `waypoint_passes` (default `24`). |
| in | `rounds` | Int | n/a | no | Field `rounds` (default `16`). |
| in | `min_abs_cost_improvement` | Float64 | n/a | no | Field `min_abs_cost_improvement` (default `1.0e-8`). |
| in | `min_rel_cost_improvement` | Float64 | n/a | no | Field `min_rel_cost_improvement` (default `1.0e-5`). |
| in | `insert_straight_waypoints` | Bool | n/a | no | Field `insert_straight_waypoints` (default `true`). |
| in | `straight_max_segment_length_m` | Float64 | n/a | no | Field `straight_max_segment_length_m` (default `2.0`). |
| in | `straight_max_inserted` | Int | n/a | no | Field `straight_max_inserted` (default `12`). |
| in | `straight_clearance_margin_m` | Float64 | n/a | no | Field `straight_clearance_margin_m` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPSORefinementSettings | n/a | — | Constructed `RPOPSORefinementSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_parameters_rpopsoconfigurator|RPOPSOConfigurator]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:181-181`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Validation is external: `period > 0`, `sample_ds_m > 0`, `straight_max_segment_length_m > 0`, and non-negative counts and thresholds. The refinement `sample_ds_m` is replaced by `safe_distance_m` when that is positive (`rpo_hypr_refinement_sampling_density_m`), so the 0.025 m default rarely takes effect in keep-out scenarios.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 141.
