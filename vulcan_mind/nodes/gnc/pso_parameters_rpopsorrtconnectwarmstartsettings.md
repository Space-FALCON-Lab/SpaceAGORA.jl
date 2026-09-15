---
id: gnc.pso_parameters_rpopsorrtconnectwarmstartsettings
label: RPOPSORRTConnectWarmstartSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: RPOPSORRTConnectWarmstartSettings
  lines:
  - 128
  - 128
inputs:
- id: enabled
  type: Bool
  units: n/a
  required: false
  description: Field `enabled` (default `false`).
- id: n_iters
  type: Int
  units: n/a
  required: false
  description: Field `n_iters` (default `250`).
- id: step_size_m
  type: Float64
  units: n/a
  required: false
  description: Field `step_size_m` (default `0.75`).
- id: goal_sample_rate
  type: Float64
  units: n/a
  required: false
  description: Field `goal_sample_rate` (default `0.05`).
- id: collision_sample_ds_m
  type: Float64
  units: n/a
  required: false
  description: Field `collision_sample_ds_m` (default `0.10`).
- id: connect_max_steps
  type: Int
  units: n/a
  required: false
  description: Field `connect_max_steps` (default `10_000`).
- id: shortcut_iters
  type: Int
  units: n/a
  required: false
  description: Field `shortcut_iters` (default `40`).
- id: runtime_limit_s
  type: Float64
  units: n/a
  required: false
  description: Field `runtime_limit_s` (default `Inf`).
- id: box_margin_m
  type: Float64
  units: n/a
  required: false
  description: Field `box_margin_m` (default `0.75`).
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
  type: RPOPSORRTConnectWarmstartSettings
  units: n/a
  description: Constructed `RPOPSORRTConnectWarmstartSettings` (keyword constructor
    via @kwdef).
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

# RPOPSORRTConnectWarmstartSettings

## Purpose
Grouped struct for an optional RRT-Connect warm start that finds a collision-free seed path before PSO refinement, giving the swarm a feasible starting point in cluttered geometry.

## Design & Implementation
Fields: `enabled::Bool = false`; `n_iters::Int = 250` tree-growth iterations; `step_size_m::Float64 = 0.75` extension length; `goal_sample_rate::Float64 = 0.05` probability of sampling the goal; `collision_sample_ds_m::Float64 = 0.10` metres between edge collision samples; `connect_max_steps::Int = 10_000` cap on greedy connect extensions; `shortcut_iters::Int = 40` post-processing shortcut attempts; `runtime_limit_s::Float64 = Inf`; `box_margin_m::Float64 = 0.75` sampling-box padding. Mapped to `rrt_warmstart_*` in `RPOPSOConfig`, with `n_iters` renamed to `rrt_warmstart_iters`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `enabled` | Bool | n/a | no | Field `enabled` (default `false`). |
| in | `n_iters` | Int | n/a | no | Field `n_iters` (default `250`). |
| in | `step_size_m` | Float64 | n/a | no | Field `step_size_m` (default `0.75`). |
| in | `goal_sample_rate` | Float64 | n/a | no | Field `goal_sample_rate` (default `0.05`). |
| in | `collision_sample_ds_m` | Float64 | n/a | no | Field `collision_sample_ds_m` (default `0.10`). |
| in | `connect_max_steps` | Int | n/a | no | Field `connect_max_steps` (default `10_000`). |
| in | `shortcut_iters` | Int | n/a | no | Field `shortcut_iters` (default `40`). |
| in | `runtime_limit_s` | Float64 | n/a | no | Field `runtime_limit_s` (default `Inf`). |
| in | `box_margin_m` | Float64 | n/a | no | Field `box_margin_m` (default `0.75`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPSORRTConnectWarmstartSettings | n/a | — | Constructed `RPOPSORRTConnectWarmstartSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_parameters_rpopsoconfigurator|RPOPSOConfigurator]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:180-180`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Disabled by default. `validate_rpo_pso_config` enforces `goal_sample_rate` in [0, 1], positive `step_size_m`, `collision_sample_ds_m`, and `connect_max_steps`, and non-negative counts and margins. With `runtime_limit_s = Inf` a pathological scene can consume the full `n_iters * connect_max_steps` budget.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 128.
