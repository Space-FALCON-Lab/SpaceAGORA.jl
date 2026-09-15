---
id: gnc.rpo_plan_buffer_rpoplan
label: RPOPlan
kind: struct
source:
  file: src/gnc/guidance/rpo/rpo_plan_buffer.jl
  symbol: RPOPlan
  lines:
  - 2
  - 2
inputs:
- id: valid
  type: Bool
  units: n/a
  required: false
  description: Field `valid` (default `false`).
- id: t_ref_s
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `t_ref_s` (default `Float64[]`).
- id: r_ref_rtn
  type: Matrix{Float64}
  units: n/a
  required: false
  description: Field `r_ref_rtn` (default `zeros(3, 0)`).
- id: v_ref_rtn
  type: Matrix{Float64}
  units: n/a
  required: false
  description: Field `v_ref_rtn` (default `zeros(3, 0)`).
- id: path_rtn
  type: Matrix{Float64}
  units: n/a
  required: false
  description: Field `path_rtn` (default `zeros(3, 0)`).
- id: cost
  type: Float64
  units: n/a
  required: false
  description: Field `cost` (default `Inf`).
- id: diagnostics
  type: NamedTuple
  units: n/a
  required: false
  description: Field `diagnostics` (default `NamedTuple()`).
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
  description: Constructed `RPOPlan` (keyword constructor via @kwdef).
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

# RPOPlan

## Purpose
Value type holding one rendezvous-and-proximity-operations trajectory plan: the retimed reference states in the target's RTN frame, the geometric path, the planner cost, and free-form diagnostics. It is what an RPO planner produces and what the guidance loop later samples to build tracking commands.

## Design & Implementation
A `Base.@kwdef mutable struct` with seven fields. `valid::Bool` defaults to `false` so a freshly constructed plan is never mistaken for a usable one. `t_ref_s::Vector{Float64}` holds the reference epochs in seconds; `r_ref_rtn` and `v_ref_rtn` are `Matrix{Float64}` sized `3 x N` (metres and metres per second, radial/along-track/cross-track), initialised as `zeros(3, 0)`; `path_rtn` is the same shape for the geometric path. `cost::Float64` defaults to `Inf`, the natural identity for a minimisation, and `diagnostics::NamedTuple` defaults to an empty `NamedTuple()` so planners can attach solver detail without changing the type.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `valid` | Bool | n/a | no | Field `valid` (default `false`). |
| in | `t_ref_s` | Vector{Float64} | n/a | no | Field `t_ref_s` (default `Float64[]`). |
| in | `r_ref_rtn` | Matrix{Float64} | n/a | no | Field `r_ref_rtn` (default `zeros(3, 0)`). |
| in | `v_ref_rtn` | Matrix{Float64} | n/a | no | Field `v_ref_rtn` (default `zeros(3, 0)`). |
| in | `path_rtn` | Matrix{Float64} | n/a | no | Field `path_rtn` (default `zeros(3, 0)`). |
| in | `cost` | Float64 | n/a | no | Field `cost` (default `Inf`). |
| in | `diagnostics` | NamedTuple | n/a | no | Field `diagnostics` (default `NamedTuple()`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPlan | n/a | — | Constructed `RPOPlan` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.replanning_rpo_plan_from_path|rpo_plan_from_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:251-251`
- [[gnc.rpo_guidance_hooks_build_rpo_plan_from_start|build_rpo_plan_from_start]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_guidance_hooks.jl:36-36`
- [[gnc.rpo_plan_buffer_rpoplanbuffer|RPOPlanBuffer]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_plan_buffer.jl:15-15`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/rpo_plan_buffer.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Nothing checks that `t_ref_s`, `r_ref_rtn`, `v_ref_rtn` and `path_rtn` share a column count, so a partially populated plan can be marked valid and read out of bounds downstream. The struct is mutable and unsynchronised, so a plan handed to a running guidance loop while a planner still writes into it is a data race. Because `diagnostics` is an untyped `NamedTuple` field, every access to it is type-unstable.

## Provenance
Mapped from `src/gnc/guidance/rpo/rpo_plan_buffer.jl` line 2.
