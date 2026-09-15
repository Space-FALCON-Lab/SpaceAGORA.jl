---
id: gnc.pso_parameters_rpopsoschedulesettings
label: RPOPSOScheduleSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: RPOPSOScheduleSettings
  lines:
  - 78
  - 78
inputs:
- id: enabled
  type: Bool
  units: n/a
  required: false
  description: Field `enabled` (default `true`).
- id: w_end_fraction
  type: Float64
  units: n/a
  required: false
  description: Field `w_end_fraction` (default `0.65`).
- id: c1_end_fraction
  type: Float64
  units: n/a
  required: false
  description: Field `c1_end_fraction` (default `0.75`).
- id: c2_end_fraction
  type: Float64
  units: n/a
  required: false
  description: Field `c2_end_fraction` (default `1.25`).
- id: transition_fraction
  type: Float64
  units: n/a
  required: false
  description: Field `transition_fraction` (default `0.45`).
- id: w_min
  type: Float64
  units: n/a
  required: false
  description: Field `w_min` (default `0.25`).
- id: c_min
  type: Float64
  units: n/a
  required: false
  description: Field `c_min` (default `0.5`).
- id: c_max
  type: Float64
  units: n/a
  required: false
  description: Field `c_max` (default `2.5`).
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
  type: RPOPSOScheduleSettings
  units: n/a
  description: Constructed `RPOPSOScheduleSettings` (keyword constructor via @kwdef).
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

# RPOPSOScheduleSettings

## Purpose
Grouped struct defining the iteration schedule that anneals PSO inertia `w` and acceleration coefficients `c1`, `c2` from their base values toward end-of-run targets, shifting the swarm from exploration to exploitation.

## Design & Implementation
Fields: `enabled::Bool = true`; `w_end_fraction = 0.65`, `c1_end_fraction = 0.75`, `c2_end_fraction = 1.25`, multipliers applied to the base coefficients at the end of the run; `transition_fraction = 0.45`, fraction of `n_iters` over which the interpolation completes; `w_min = 0.25`, `c_min = 0.5`, `c_max = 2.5`, absolute clamps on the scheduled values. Copied into `schedule_*` fields of `RPOPSOConfig`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `enabled` | Bool | n/a | no | Field `enabled` (default `true`). |
| in | `w_end_fraction` | Float64 | n/a | no | Field `w_end_fraction` (default `0.65`). |
| in | `c1_end_fraction` | Float64 | n/a | no | Field `c1_end_fraction` (default `0.75`). |
| in | `c2_end_fraction` | Float64 | n/a | no | Field `c2_end_fraction` (default `1.25`). |
| in | `transition_fraction` | Float64 | n/a | no | Field `transition_fraction` (default `0.45`). |
| in | `w_min` | Float64 | n/a | no | Field `w_min` (default `0.25`). |
| in | `c_min` | Float64 | n/a | no | Field `c_min` (default `0.5`). |
| in | `c_max` | Float64 | n/a | no | Field `c_max` (default `2.5`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPSOScheduleSettings | n/a | — | Constructed `RPOPSOScheduleSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_parameters_rpopsoconfigurator|RPOPSOConfigurator]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:175-175`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`validate_rpo_pso_config` requires `transition_fraction > 0`, `w_min >= 0`, `c_min >= 0`, and `c_min <= c_max`; a `transition_fraction` greater than 1 is accepted and simply never completes the transition. The end fractions are unbounded, so `c2_end_fraction` can push `c2` to `c_max` early in the run.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 78.
