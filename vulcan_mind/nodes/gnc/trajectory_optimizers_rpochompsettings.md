---
id: gnc.trajectory_optimizers_rpochompsettings
label: RPOCHOMPSettings
kind: struct
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: RPOCHOMPSettings
  lines:
  - 4
  - 4
inputs:
- id: n_iters
  type: Int
  units: n/a
  required: false
  description: Field `n_iters` (default `100`).
- id: learning_rate
  type: Float64
  units: n/a
  required: false
  description: Field `learning_rate` (default `0.06`).
- id: gradient_eps
  type: Float64
  units: n/a
  required: false
  description: Field `gradient_eps` (default `1.0e-3`).
- id: w_smooth
  type: Float64
  units: n/a
  required: false
  description: Field `w_smooth` (default `1.0`).
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
  type: RPOCHOMPSettings
  units: n/a
  description: Constructed `RPOCHOMPSettings` (keyword constructor via @kwdef).
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

# RPOCHOMPSettings

## Purpose
Immutable `Base.@kwdef` settings struct for the CHOMP-like gradient-descent trajectory optimizer used in RPO planner comparisons. It holds the four knobs that `rpo_chomp_plan_path` reads: iteration budget, step size, finite-difference perturbation, and smoothness weight.

## Design & Implementation
Fields with defaults: `n_iters::Int = 100` (outer descent iterations), `learning_rate::Float64 = 0.06` (multiplied by a backtracking scale of 1.0 down to 0.0625 before being applied to the metric-preconditioned gradient), `gradient_eps::Float64 = 1.0e-3` (absolute finite-difference step in metres, floored per component by `rpo_chomp_numeric_gradient`), and `w_smooth::Float64 = 1.0` (weight on the normalised squared second-difference term inside `rpo_trajectory_soft_objective`). `rpo_chomp_plan_path` constructs the default as `RPOCHOMPSettings(n_iters=cfg.n_iters)` so the iteration count follows the PSO config unless overridden.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `n_iters` | Int | n/a | no | Field `n_iters` (default `100`). |
| in | `learning_rate` | Float64 | n/a | no | Field `learning_rate` (default `0.06`). |
| in | `gradient_eps` | Float64 | n/a | no | Field `gradient_eps` (default `1.0e-3`). |
| in | `w_smooth` | Float64 | n/a | no | Field `w_smooth` (default `1.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOCHOMPSettings | n/a | — | Constructed `RPOCHOMPSettings` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpoplannercomparisonconfig|RPOPlannerComparisonConfig]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:35-35`
- [[gnc.planner_comparison_runtime_limited_iters|runtime_limited_iters]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:364-364`
- [[gncy.trajectory_optimizers_rpo_chomp_plan_path|rpo_chomp_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:231-231`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No field validation occurs at construction: a non-positive `n_iters` is silently lifted to 1 by `max(1, Int(settings.n_iters))` in the planner, and a negative `learning_rate` would turn descent into ascent without error. `gradient_eps` is in absolute metres, so very large RTN corridors would benefit from scaling that is not applied here.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl` line 4.
