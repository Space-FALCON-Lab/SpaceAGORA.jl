---
id: gnc.pso_parameters_rpoadaptivesamplingsettings
label: RPOAdaptiveSamplingSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: RPOAdaptiveSamplingSettings
  lines:
  - 59
  - 59
inputs:
- id: enabled
  type: Bool
  units: n/a
  required: false
  description: Field `enabled` (default `true`).
- id: max_ds_m
  type: Float64
  units: n/a
  required: false
  description: Field `max_ds_m` (default `0.50`).
- id: far_clearance_m
  type: Float64
  units: n/a
  required: false
  description: Field `far_clearance_m` (default `1.0`).
- id: power
  type: Float64
  units: n/a
  required: false
  description: Field `power` (default `1.0`).
- id: safe_distance_fraction
  type: Float64
  units: n/a
  required: false
  description: Field `safe_distance_fraction` (default `0.5`).
- id: obstacle_guard_fraction
  type: Float64
  units: n/a
  required: false
  description: Field `obstacle_guard_fraction` (default `0.5`).
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
  type: RPOAdaptiveSamplingSettings
  units: n/a
  description: Constructed `RPOAdaptiveSamplingSettings` (keyword constructor via
    @kwdef).
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

# RPOAdaptiveSamplingSettings

## Purpose
Grouped struct controlling adaptive collision-sample spacing along candidate paths, allowing the HYPR evaluator to sample coarsely far from the station and densely near it.

## Design & Implementation
Fields: `enabled::Bool = true`; `max_ds_m::Float64 = 0.50` metres, the coarsest spacing used at large clearance; `far_clearance_m::Float64 = 1.0` metres beyond which spacing saturates at `max_ds_m`; `power::Float64 = 1.0` exponent shaping the clearance-to-spacing curve; `safe_distance_fraction::Float64 = 0.5` and `obstacle_guard_fraction::Float64 = 0.5` scale the keep-out distance to decide when to fall back to the fine base spacing. Flattened into `adaptive_sampling_*` fields of `RPOPSOConfig`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `enabled` | Bool | n/a | no | Field `enabled` (default `true`). |
| in | `max_ds_m` | Float64 | n/a | no | Field `max_ds_m` (default `0.50`). |
| in | `far_clearance_m` | Float64 | n/a | no | Field `far_clearance_m` (default `1.0`). |
| in | `power` | Float64 | n/a | no | Field `power` (default `1.0`). |
| in | `safe_distance_fraction` | Float64 | n/a | no | Field `safe_distance_fraction` (default `0.5`). |
| in | `obstacle_guard_fraction` | Float64 | n/a | no | Field `obstacle_guard_fraction` (default `0.5`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOAdaptiveSamplingSettings | n/a | — | Constructed `RPOAdaptiveSamplingSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_parameters_rpopsoconfigurator|RPOPSOConfigurator]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:173-173`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No validation here; `validate_rpo_pso_config` requires `max_ds_m`, `far_clearance_m`, `power`, and `safe_distance_fraction` to be positive and `obstacle_guard_fraction` in (0, 1]. Coarse sampling near thin station features can miss a collision between samples if `max_ds_m` exceeds the feature size, which nothing in this struct guards against.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 59.
