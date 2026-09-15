---
id: gnc.pso_parameters_rpopsoswarmsettings
label: RPOPSOSwarmSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: RPOPSOSwarmSettings
  lines:
  - 2
  - 2
inputs:
- id: n_waypoints
  type: Int
  units: n/a
  required: false
  description: Field `n_waypoints` (default `5`).
- id: n_particles
  type: Int
  units: n/a
  required: false
  description: Field `n_particles` (default `200`).
- id: n_iters
  type: Int
  units: n/a
  required: false
  description: Field `n_iters` (default `55`).
- id: spread_scale
  type: Float64
  units: n/a
  required: false
  description: Field `spread_scale` (default `0.2`).
- id: search_margin_m
  type: Float64
  units: n/a
  required: false
  description: Field `search_margin_m` (default `10.0`).
- id: sample_ds_m
  type: Float64
  units: n/a
  required: false
  description: Field `sample_ds_m` (default `0.05`).
- id: curve_type
  type: Symbol
  units: n/a
  required: false
  description: Field `curve_type` (default `:bezier`).
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
  type: RPOPSOSwarmSettings
  units: n/a
  description: Constructed `RPOPSOSwarmSettings` (keyword constructor via @kwdef).
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

# RPOPSOSwarmSettings

## Purpose
Grouped `Base.@kwdef` struct holding the core particle-swarm sizing and geometry-sampling parameters for the RPO HYPR planner: how many waypoints define a candidate path, how many particles and iterations the swarm runs, and how the path is discretised for collision checks.

## Design & Implementation
Fields: `n_waypoints::Int = 5` internal waypoints per particle; `n_particles::Int = 200`; `n_iters::Int = 55`; `spread_scale::Float64 = 0.2` fraction of the search box used to seed initial particle dispersion; `search_margin_m::Float64 = 10.0` metres added around the start-goal bounding box; `sample_ds_m::Float64 = 0.05` metres between collision samples along the curve; `curve_type::Symbol = :bezier` (or `:polyline`). `RPOPSOConfig(configurator)` copies each field into the flattened config with the same names.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `n_waypoints` | Int | n/a | no | Field `n_waypoints` (default `5`). |
| in | `n_particles` | Int | n/a | no | Field `n_particles` (default `200`). |
| in | `n_iters` | Int | n/a | no | Field `n_iters` (default `55`). |
| in | `spread_scale` | Float64 | n/a | no | Field `spread_scale` (default `0.2`). |
| in | `search_margin_m` | Float64 | n/a | no | Field `search_margin_m` (default `10.0`). |
| in | `sample_ds_m` | Float64 | n/a | no | Field `sample_ds_m` (default `0.05`). |
| in | `curve_type` | Symbol | n/a | no | Field `curve_type` (default `:bezier`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPSOSwarmSettings | n/a | — | Constructed `RPOPSOSwarmSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_parameters_rpopsoconfigurator|RPOPSOConfigurator]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:170-170`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No validation at construction; range checks happen only when the values reach `validate_rpo_pso_config` through `RPOPSOConfig(configurator)`. `sample_ds_m` is overridden by `safe_distance_m` whenever the latter is positive (see `_rpo_pso_sync_sample_ds_with_safe_distance`), so the default 0.05 m is frequently not what the planner actually uses.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 2.
