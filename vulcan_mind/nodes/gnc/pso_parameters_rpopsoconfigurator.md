---
id: gnc.pso_parameters_rpopsoconfigurator
label: RPOPSOConfigurator
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: RPOPSOConfigurator
  lines:
  - 169
  - 169
inputs:
- id: swarm
  type: RPOPSOSwarmSettings
  units: n/a
  required: false
  description: Field `swarm` (default `RPOPSOSwarmSettings()`).
- id: objective
  type: RPOPSOObjectiveSettings
  units: n/a
  required: false
  description: Field `objective` (default `RPOPSOObjectiveSettings()`).
- id: adaptive
  type: RPOPSOAdaptiveSettings
  units: n/a
  required: false
  description: Field `adaptive` (default `RPOPSOAdaptiveSettings()`).
- id: adaptive_sampling
  type: RPOAdaptiveSamplingSettings
  units: n/a
  required: false
  description: Field `adaptive_sampling` (default `RPOAdaptiveSamplingSettings()`).
- id: cull
  type: RPOPSOCullSettings
  units: n/a
  required: false
  description: Field `cull` (default `RPOPSOCullSettings()`).
- id: schedule
  type: RPOPSOScheduleSettings
  units: n/a
  required: false
  description: Field `schedule` (default `RPOPSOScheduleSettings()`).
- id: stagnation
  type: RPOPSOStagnationSettings
  units: n/a
  required: false
  description: Field `stagnation` (default `RPOPSOStagnationSettings()`).
- id: early_stopping
  type: RPOPSOEarlyStoppingSettings
  units: n/a
  required: false
  description: Field `early_stopping` (default `RPOPSOEarlyStoppingSettings()`).
- id: probe
  type: RPOPSOProbeSettings
  units: n/a
  required: false
  description: Field `probe` (default `RPOPSOProbeSettings()`).
- id: reexplore
  type: RPOPSOReexploreSettings
  units: n/a
  required: false
  description: Field `reexplore` (default `RPOPSOReexploreSettings()`).
- id: rrt_warmstart
  type: RPOPSORRTConnectWarmstartSettings
  units: n/a
  required: false
  description: Field `rrt_warmstart` (default `RPOPSORRTConnectWarmstartSettings()`).
- id: refinement
  type: RPOPSORefinementSettings
  units: n/a
  required: false
  description: Field `refinement` (default `RPOPSORefinementSettings()`).
- id: retiming
  type: RPOPSORetimingSettings
  units: n/a
  required: false
  description: Field `retiming` (default `RPOPSORetimingSettings()`).
- id: safe_distance_m
  type: Float64
  units: n/a
  required: false
  description: Field `safe_distance_m` (default `0.0`).
- id: goal_collision_margin_m
  type: Float64
  units: n/a
  required: false
  description: Field `goal_collision_margin_m` (default `0.0`).
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
  type: RPOPSOConfigurator
  units: n/a
  description: Constructed `RPOPSOConfigurator` (keyword constructor via @kwdef).
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

# RPOPSOConfigurator

## Purpose
Hierarchical `Base.@kwdef` container that composes all thirteen grouped RPO PSO setting structs plus the two top-level keep-out distances, giving users a structured way to build a planner configuration before it is flattened into `RPOPSOConfig`.

## Design & Implementation
Fields, each defaulting to the group's own default constructor: `swarm`, `objective`, `adaptive`, `adaptive_sampling`, `cull`, `schedule`, `stagnation`, `early_stopping`, `probe`, `reexplore`, `rrt_warmstart`, `refinement`, `retiming`, followed by `safe_distance_m::Float64 = 0.0` and `goal_collision_margin_m::Float64 = 0.0`. The constructor method `RPOPSOConfig(configurator; kwargs...)` reads every group field explicitly (about 115 keyword assignments) and finishes with `rpo_pso_config(cfg; kwargs...)`, so keyword overrides and validation apply after flattening. `rpo_pso_config(configurator; kwargs...)` is the convenience alias.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `swarm` | RPOPSOSwarmSettings | n/a | no | Field `swarm` (default `RPOPSOSwarmSettings()`). |
| in | `objective` | RPOPSOObjectiveSettings | n/a | no | Field `objective` (default `RPOPSOObjectiveSettings()`). |
| in | `adaptive` | RPOPSOAdaptiveSettings | n/a | no | Field `adaptive` (default `RPOPSOAdaptiveSettings()`). |
| in | `adaptive_sampling` | RPOAdaptiveSamplingSettings | n/a | no | Field `adaptive_sampling` (default `RPOAdaptiveSamplingSettings()`). |
| in | `cull` | RPOPSOCullSettings | n/a | no | Field `cull` (default `RPOPSOCullSettings()`). |
| in | `schedule` | RPOPSOScheduleSettings | n/a | no | Field `schedule` (default `RPOPSOScheduleSettings()`). |
| in | `stagnation` | RPOPSOStagnationSettings | n/a | no | Field `stagnation` (default `RPOPSOStagnationSettings()`). |
| in | `early_stopping` | RPOPSOEarlyStoppingSettings | n/a | no | Field `early_stopping` (default `RPOPSOEarlyStoppingSettings()`). |
| in | `probe` | RPOPSOProbeSettings | n/a | no | Field `probe` (default `RPOPSOProbeSettings()`). |
| in | `reexplore` | RPOPSOReexploreSettings | n/a | no | Field `reexplore` (default `RPOPSOReexploreSettings()`). |
| in | `rrt_warmstart` | RPOPSORRTConnectWarmstartSettings | n/a | no | Field `rrt_warmstart` (default `RPOPSORRTConnectWarmstartSettings()`). |
| in | `refinement` | RPOPSORefinementSettings | n/a | no | Field `refinement` (default `RPOPSORefinementSettings()`). |
| in | `retiming` | RPOPSORetimingSettings | n/a | no | Field `retiming` (default `RPOPSORetimingSettings()`). |
| in | `safe_distance_m` | Float64 | n/a | no | Field `safe_distance_m` (default `0.0`). |
| in | `goal_collision_margin_m` | Float64 | n/a | no | Field `goal_collision_margin_m` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPSOConfigurator | n/a | — | Constructed `RPOPSOConfigurator` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- `callees` → [[gnc.pso_parameters_rpoadaptivesamplingsettings|RPOAdaptiveSamplingSettings]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:173-173`
- `callees` → [[gnc.pso_parameters_rpopsoadaptivesettings|RPOPSOAdaptiveSettings]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:172-172`
- `callees` → [[gnc.pso_parameters_rpopsocullsettings|RPOPSOCullSettings]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:174-174`
- `callees` → [[gnc.pso_parameters_rpopsoearlystoppingsettings|RPOPSOEarlyStoppingSettings]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:177-177`
- `callees` → [[gnc.pso_parameters_rpopsoobjectivesettings|RPOPSOObjectiveSettings]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:171-171`
- `callees` → [[gnc.pso_parameters_rpopsoprobesettings|RPOPSOProbeSettings]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:178-178`
- `callees` → [[gnc.pso_parameters_rpopsoreexploresettings|RPOPSOReexploreSettings]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:179-179`
- `callees` → [[gnc.pso_parameters_rpopsorefinementsettings|RPOPSORefinementSettings]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:181-181`
- `callees` → [[gnc.pso_parameters_rpopsoretimingsettings|RPOPSORetimingSettings]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:182-182`
- `callees` → [[gnc.pso_parameters_rpopsorrtconnectwarmstartsettings|RPOPSORRTConnectWarmstartSettings]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:180-180`
- `callees` → [[gnc.pso_parameters_rpopsoschedulesettings|RPOPSOScheduleSettings]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:175-175`
- `callees` → [[gnc.pso_parameters_rpopsostagnationsettings|RPOPSOStagnationSettings]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:176-176`
- `callees` → [[gnc.pso_parameters_rpopsoswarmsettings|RPOPSOSwarmSettings]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:170-170`
<!-- vulcan:connections:end -->

## Limitations
The flattening is a hand-maintained field list; adding a field to a group struct without updating `RPOPSOConfig(configurator)` silently drops it. `iteration_runtime_limit_s` has no group counterpart and can only be set via keyword. Validation (`ArgumentError` on out-of-range values) occurs only during flattening, not at configurator construction.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 169.
