---
id: gnc.planner_comparison_rpoplannercomparisoncase
label: RPOPlannerComparisonCase
kind: struct
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: RPOPlannerComparisonCase
  lines:
  - 21
  - 21
inputs:
- id: start_rtn
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `start_rtn`.
- id: goal_rtn
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `goal_rtn`.
- id: label
  type: String
  units: n/a
  required: false
  description: Field `label` (default `"case"`).
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
  type: RPOPlannerComparisonCase
  units: n/a
  description: Constructed `RPOPlannerComparisonCase` (keyword constructor via @kwdef).
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

# RPOPlannerComparisonCase

## Purpose
Minimal record describing one start-goal pair in the RTN frame for a planner comparison batch, with a label used in CSV rows, plot legends, and artifact file names.

## Design & Implementation
`Base.@kwdef struct RPOPlannerComparisonCase` with required `start_rtn::SVector{3, Float64}` and `goal_rtn::SVector{3, Float64}` (metres in radial, along-track, cross-track) and optional `label::String = "case"`. `rpo_run_planner_comparison_batch` iterates a collection of these, and the case is attached to each plan tuple as `case` so plotting functions can read `plan.case.start_rtn`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `start_rtn` | SVector{3, Float64} | n/a | yes | Field `start_rtn`. |
| in | `goal_rtn` | SVector{3, Float64} | n/a | yes | Field `goal_rtn`. |
| in | `label` | String | n/a | no | Field `label` (default `"case"`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPlannerComparisonCase | n/a | — | Constructed `RPOPlannerComparisonCase` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The default label "case" is not unique, so batches built without explicit labels produce indistinguishable legend entries and identical `rpo_comparison_artifact_slug` tokens. No check that start and goal differ or lie outside the station keep-out volume; those failures surface as planner errors or failed tracking.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 21.
