---
id: core.runtime_types_controller
label: Controller
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: Controller
  lines:
  - 370
  - 370
inputs:
- id: guidance_t_eval
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `guidance_t_eval` (default `[]`).
- id: count_controller
  type: Int64
  units: n/a
  required: false
  description: Field `count_controller` (default `1`).
- id: count_prev_controller
  type: Int64
  units: n/a
  required: false
  description: Field `count_prev_controller` (default `0`).
- id: stored_state
  type: Int64
  units: n/a
  required: false
  description: Field `stored_state` (default `1`).
- id: prev_time
  type: Float64
  units: n/a
  required: false
  description: Field `prev_time` (default `0.0`).
- id: t
  type: Float64
  units: n/a
  required: false
  description: Field `t` (default `0.0`).
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
  type: Controller
  units: n/a
  description: Constructed `Controller` (keyword constructor via @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# Controller

## Purpose
Small mutable bookkeeping record for the legacy guidance scheduler, tracking when the controller last ran and how many evaluations have occurred.

## Design & Implementation
`@kwdef mutable struct Controller` with `guidance_t_eval::Vector{Float64}` (scheduled evaluation times, s), `count_controller::Int64 = 1`, `count_prev_controller::Int64 = 0`, `stored_state::Int64 = 1`, `prev_time::Float64 = 0.0`, and `t::Float64 = 0.0`. Legacy callbacks increment the counters and compare `t` with `prev_time` to decide whether a new guidance solution is due.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `guidance_t_eval` | Vector{Float64} | n/a | no | Field `guidance_t_eval` (default `[]`). |
| in | `count_controller` | Int64 | n/a | no | Field `count_controller` (default `1`). |
| in | `count_prev_controller` | Int64 | n/a | no | Field `count_prev_controller` (default `0`). |
| in | `stored_state` | Int64 | n/a | no | Field `stored_state` (default `1`). |
| in | `prev_time` | Float64 | n/a | no | Field `prev_time` (default `0.0`). |
| in | `t` | Float64 | n/a | no | Field `t` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Controller | n/a | — | Constructed `Controller` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Counters start at 1 and 0 respectively with no documented meaning for the offset. Nothing prevents `guidance_t_eval` from being unsorted, which the legacy scheduler assumes. The struct is unused by the modern effector-based control path.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 370.
