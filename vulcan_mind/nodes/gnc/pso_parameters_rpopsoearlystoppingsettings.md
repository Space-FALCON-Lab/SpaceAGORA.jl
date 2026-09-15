---
id: gnc.pso_parameters_rpopsoearlystoppingsettings
label: RPOPSOEarlyStoppingSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/hypr/pso_parameters.jl
  symbol: RPOPSOEarlyStoppingSettings
  lines:
  - 98
  - 98
inputs:
- id: enabled
  type: Bool
  units: n/a
  required: false
  description: Field `enabled` (default `false`).
- id: patience
  type: Int
  units: n/a
  required: false
  description: Field `patience` (default `10`).
- id: min_iters
  type: Int
  units: n/a
  required: false
  description: Field `min_iters` (default `12`).
- id: min_abs_improvement
  type: Float64
  units: n/a
  required: false
  description: Field `min_abs_improvement` (default `1.0e-8`).
- id: min_rel_improvement
  type: Float64
  units: n/a
  required: false
  description: Field `min_rel_improvement` (default `1.0e-4`).
- id: require_feasible
  type: Bool
  units: n/a
  required: false
  description: Field `require_feasible` (default `true`).
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
  type: RPOPSOEarlyStoppingSettings
  units: n/a
  description: Constructed `RPOPSOEarlyStoppingSettings` (keyword constructor via
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

# RPOPSOEarlyStoppingSettings

## Purpose
Grouped struct governing optional early termination of the PSO loop once the best feasible cost has stopped improving by a meaningful absolute or relative amount.

## Design & Implementation
Fields: `enabled::Bool = false` (off by default so benchmark runs use the full iteration budget); `patience::Int = 10` non-improving iterations tolerated; `min_iters::Int = 12` iterations that must elapse before stopping is considered; `min_abs_improvement::Float64 = 1.0e-8` and `min_rel_improvement::Float64 = 1.0e-4` thresholds below which an improvement does not reset patience; `require_feasible::Bool = true` so stopping only counts iterations whose best path is collision-free. Flattened to `early_stopping_*` in `RPOPSOConfig`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `enabled` | Bool | n/a | no | Field `enabled` (default `false`). |
| in | `patience` | Int | n/a | no | Field `patience` (default `10`). |
| in | `min_iters` | Int | n/a | no | Field `min_iters` (default `12`). |
| in | `min_abs_improvement` | Float64 | n/a | no | Field `min_abs_improvement` (default `1.0e-8`). |
| in | `min_rel_improvement` | Float64 | n/a | no | Field `min_rel_improvement` (default `1.0e-4`). |
| in | `require_feasible` | Bool | n/a | no | Field `require_feasible` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPSOEarlyStoppingSettings | n/a | — | Constructed `RPOPSOEarlyStoppingSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_parameters_rpopsoconfigurator|RPOPSOConfigurator]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl:177-177`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_parameters.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because `enabled` defaults to false, the other fields are inert unless a caller opts in. `validate_rpo_pso_config` only checks non-negativity of `patience`, `min_iters`, and the two improvement thresholds. With `require_feasible = true` a run that never finds a feasible path can never stop early.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_parameters.jl` line 98.
