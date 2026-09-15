---
id: gnc.pso_parameters_rpopsoreexploresettings
label: RPOPSOReexploreSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: RPOPSOReexploreSettings
  lines:
  - 118
  - 118
inputs:
- id: enabled
  type: Bool
  units: n/a
  required: false
  description: Field `enabled` (default `true`).
- id: trigger_iter
  type: Int
  units: n/a
  required: false
  description: Field `trigger_iter` (default `10`).
- id: search_margin_scale
  type: Float64
  units: n/a
  required: false
  description: Field `search_margin_scale` (default `2.0`).
- id: waypoint_scale
  type: Float64
  units: n/a
  required: false
  description: Field `waypoint_scale` (default `1.5`).
- id: waypoint_increment
  type: Int
  units: n/a
  required: false
  description: Field `waypoint_increment` (default `2`).
- id: max_waypoints
  type: Int
  units: n/a
  required: false
  description: Field `max_waypoints` (default `max(8 + 4, Int(ceil(1.5 * 8)))`).
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
  type: RPOPSOReexploreSettings
  units: n/a
  description: Constructed `RPOPSOReexploreSettings` (keyword constructor via @kwdef).
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

# RPOPSOReexploreSettings

## Purpose
Grouped struct controlling the re-exploration response when the swarm stagnates with an infeasible best: the search box is widened and more waypoints are added so the swarm can route around obstacles it could not previously clear.

## Design & Implementation
Fields: `enabled::Bool = true`; `trigger_iter::Int = 10` earliest iteration for re-exploration; `search_margin_scale::Float64 = 2.0` multiplier on `search_margin_m`; `waypoint_scale::Float64 = 1.5` multiplicative growth of `n_waypoints`; `waypoint_increment::Int = 2` additive growth; `max_waypoints::Int = max(8 + 4, Int(ceil(1.5 * 8)))`, which evaluates to 12 and caps the grown waypoint count relative to the adaptive maximum of 8. Mapped to `reexplore_*` in `RPOPSOConfig`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `enabled` | Bool | n/a | no | Field `enabled` (default `true`). |
| in | `trigger_iter` | Int | n/a | no | Field `trigger_iter` (default `10`). |
| in | `search_margin_scale` | Float64 | n/a | no | Field `search_margin_scale` (default `2.0`). |
| in | `waypoint_scale` | Float64 | n/a | no | Field `waypoint_scale` (default `1.5`). |
| in | `waypoint_increment` | Int | n/a | no | Field `waypoint_increment` (default `2`). |
| in | `max_waypoints` | Int | n/a | no | Field `max_waypoints` (default `max(8 + 4, Int(ceil(1.5 * 8)))`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPSOReexploreSettings | n/a | — | Constructed `RPOPSOReexploreSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_parameters_rpopsoconfigurator|RPOPSOConfigurator]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:179-179`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The `max_waypoints` default hard-codes the adaptive default `n_waypoints_max = 8`; if a caller raises `RPOPSOAdaptiveSettings.n_waypoints_max` the reexplore cap does not follow automatically. Validation (positive scales, non-negative counts) is performed only in `validate_rpo_pso_config`.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 118.
