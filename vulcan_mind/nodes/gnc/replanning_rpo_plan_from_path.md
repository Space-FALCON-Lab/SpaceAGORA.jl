---
id: gnc.replanning_rpo_plan_from_path
label: rpo_plan_from_path
kind: function
source:
  file: src/gnc/guidance/rpo/hypr_planning/replanning.jl
  symbol: rpo_plan_from_path
  lines:
  - 249
  - 249
inputs:
- id: path_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `path_rtn`.
- id: geometry
  type: Any
  units: n/a
  required: true
  description: Positional argument `geometry`.
- id: cfg
  type: RPOPSOConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: safe_distance_m
  type: Real
  units: n/a
  required: true
  description: Positional argument `safe_distance_m`.
- id: t
  type: Real
  units: n/a
  required: true
  description: Positional argument `t`.
- id: cost
  type: Any
  units: n/a
  required: false
  description: Keyword argument `cost` (default `NaN`).
- id: diagnostics
  type: Any
  units: n/a
  required: false
  description: Keyword argument `diagnostics` (default `NamedTuple()`).
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
  type: RPOPlan
  units: n/a
  description: Return value of `rpo_plan_from_path`. Returns `RPOPlan(`.
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

# rpo_plan_from_path

## Purpose
Wraps a freshly planned or retimed geometric path into a complete `RPOPlan` by generating its time-parameterised position and velocity references, recording the planning time in diagnostics. Both fresh replans and retimes funnel through it so the plan structure is built in one place.

## Design & Implementation
Signature `rpo_plan_from_path(path_rtn, geometry, cfg::RPOPSOConfig, safe_distance_m::Real, t::Real; cost=NaN, diagnostics=NamedTuple())`. It calls `rpo_reference_from_path(path_rtn, geometry, cfg; safe_distance_m)` to obtain `(t_ref, r_ref, v_ref)`, then constructs `RPOPlan(valid=true, t_ref_s=t_ref, r_ref_rtn=r_ref, v_ref_rtn=v_ref, path_rtn=Matrix{Float64}(path_rtn), cost=Float64(cost), diagnostics=merge(diagnostics, (planned_at_s=Float64(t),)))`. The `merge` means a caller-supplied `planned_at_s` is overwritten by the current time.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path_rtn` | Any | n/a | yes | Positional argument `path_rtn`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `safe_distance_m` | Real | n/a | yes | Positional argument `safe_distance_m`. |
| in | `t` | Real | n/a | yes | Positional argument `t`. |
| in | `cost` | Any | n/a | no | Keyword argument `cost` (default `NaN`). |
| in | `diagnostics` | Any | n/a | no | Keyword argument `diagnostics` (default `NamedTuple()`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPlan | n/a | — | Return value of `rpo_plan_from_path`. Returns `RPOPlan(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.replanning_rpo_retime_existing_plan|rpo_retime_existing_plan]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:269-269`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:257-257`
- `callees` → [[gnc.rpo_plan_buffer_rpoplan|RPOPlan]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:251-251`
- `callees` → [[gncz.rpo_reference_trajectory_rpo_reference_from_path|rpo_reference_from_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:250-250`
<!-- vulcan:connections:end -->

## Limitations
`valid` is unconditionally `true`; if `rpo_reference_from_path` returns a degenerate reference (for example a single-point path) the plan is still marked valid. `cost` defaults to `NaN`, which propagates into any downstream cost comparison and will make `<` comparisons false. `path_rtn` is copied into a dense `Matrix{Float64}`, so passing a very long resampled path allocates twice.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr_planning/replanning.jl` line 249.
