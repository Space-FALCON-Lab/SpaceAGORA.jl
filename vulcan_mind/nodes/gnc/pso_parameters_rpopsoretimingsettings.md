---
id: gnc.pso_parameters_rpopsoretimingsettings
label: RPOPSORetimingSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: RPOPSORetimingSettings
  lines:
  - 158
  - 158
inputs:
- id: dt_s
  type: Float64
  units: n/a
  required: false
  description: Field `dt_s` (default `1.0`).
- id: reaction_time_s
  type: Float64
  units: n/a
  required: false
  description: Field `reaction_time_s` (default `0.25`).
- id: a_max_mps2
  type: Float64
  units: n/a
  required: false
  description: Field `a_max_mps2` (default `0.02`).
- id: speed_scale
  type: Float64
  units: n/a
  required: false
  description: Field `speed_scale` (default `0.5`).
- id: min_speed_mps
  type: Float64
  units: n/a
  required: false
  description: Field `min_speed_mps` (default `0.0`).
- id: max_speed_mps
  type: Float64
  units: n/a
  required: false
  description: Field `max_speed_mps` (default `Inf`).
- id: max_steps
  type: Int
  units: n/a
  required: false
  description: Field `max_steps` (default `100_000`).
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
  type: RPOPSORetimingSettings
  units: n/a
  description: Constructed `RPOPSORetimingSettings` (keyword constructor via @kwdef).
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

# RPOPSORetimingSettings

## Purpose
Grouped struct holding the kinematic limits used to convert a geometric RPO path into a time-parameterised trajectory: time step, reaction time, acceleration cap, and speed bounds.

## Design & Implementation
Fields: `dt_s::Float64 = 1.0` retiming step; `reaction_time_s::Float64 = 0.25` lookahead used when limiting speed near obstacles; `a_max_mps2::Float64 = 0.02` acceleration ceiling; `speed_scale::Float64 = 0.5` fraction of the derived speed limit actually commanded; `min_speed_mps::Float64 = 0.0`; `max_speed_mps::Float64 = Inf`; `max_steps::Int = 100_000` cap on retimed samples. Copied into `retime_*` fields of `RPOPSOConfig`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `dt_s` | Float64 | n/a | no | Field `dt_s` (default `1.0`). |
| in | `reaction_time_s` | Float64 | n/a | no | Field `reaction_time_s` (default `0.25`). |
| in | `a_max_mps2` | Float64 | n/a | no | Field `a_max_mps2` (default `0.02`). |
| in | `speed_scale` | Float64 | n/a | no | Field `speed_scale` (default `0.5`). |
| in | `min_speed_mps` | Float64 | n/a | no | Field `min_speed_mps` (default `0.0`). |
| in | `max_speed_mps` | Float64 | n/a | no | Field `max_speed_mps` (default `Inf`). |
| in | `max_steps` | Int | n/a | no | Field `max_steps` (default `100_000`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPSORetimingSettings | n/a | — | Constructed `RPOPSORetimingSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_parameters_rpopsoconfigurator|RPOPSOConfigurator]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:182-182`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`validate_rpo_pso_config` requires `dt_s`, `a_max_mps2`, `speed_scale`, and `max_steps` positive, `min_speed_mps >= 0`, and `max_speed_mps >= min_speed_mps`. `max_steps = 100_000` with `dt_s = 1.0` limits retimed trajectories to roughly 27.8 hours; longer paths are truncated by the retimer rather than rejected here.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 158.
