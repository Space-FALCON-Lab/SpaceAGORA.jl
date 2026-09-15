---
id: parallel.outer_route_state_outerroutestats
label: OuterRouteStats
kind: struct
source:
  file: src/parallel/routing/outer_route_state.jl
  symbol: OuterRouteStats
  lines:
  - 60
  - 60
inputs:
- id: samples
  type: Int
  units: n/a
  required: false
  description: Field `samples` (default `0`).
- id: successes
  type: Int
  units: n/a
  required: false
  description: Field `successes` (default `0`).
- id: failures
  type: Int
  units: n/a
  required: false
  description: Field `failures` (default `0`).
- id: elapsed_sum_s
  type: Float64
  units: n/a
  required: false
  description: Field `elapsed_sum_s` (default `0.0`).
- id: elapsed_sq_sum_s
  type: Float64
  units: n/a
  required: false
  description: Field `elapsed_sq_sum_s` (default `0.0`).
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
  type: OuterRouteStats
  units: n/a
  description: Constructed `OuterRouteStats` (keyword constructor via @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# OuterRouteStats

## Purpose
Running statistics for one (workload signature, route) pair in the adaptive outer-route history. It counts samples, successes and failures and keeps the first and second moments of elapsed wall time in seconds so the router can compute a mean and variance per route without storing individual runs.

## Design & Implementation
`Base.@kwdef mutable struct` with all-zero defaults: `samples::Int`, `successes::Int`, `failures::Int`, `elapsed_sum_s::Float64`, `elapsed_sq_sum_s::Float64`. Instances are stored in `OuterRouteState.history::Dict{String, Dict{Symbol, OuterRouteStats}}` keyed by signature then route symbol (`:none`, `:threads`, `:process`). `load_outer_route_state!` merges persisted stats by adding every field, and `_route_stats_payload` clamps them to non-negative before writing TOML.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `samples` | Int | n/a | no | Field `samples` (default `0`). |
| in | `successes` | Int | n/a | no | Field `successes` (default `0`). |
| in | `failures` | Int | n/a | no | Field `failures` (default `0`). |
| in | `elapsed_sum_s` | Float64 | n/a | no | Field `elapsed_sum_s` (default `0.0`). |
| in | `elapsed_sq_sum_s` | Float64 | n/a | no | Field `elapsed_sq_sum_s` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | OuterRouteStats | n/a | — | Constructed `OuterRouteStats` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.outer_route_state__route_payload_stats|_route_payload_stats]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_state.jl:137-137`
- [[parallel.outer_route_state_load_outer_route_state_bang|load_outer_route_state!]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_state.jl:226-226`
- [[parcore.outer_route_metrics_record_outer_route_feedback_bang|record_outer_route_feedback!]] · `callees` → `callers` · call · `src/parallel/routing/outer_route_metrics.jl:48-48`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The variance recovered from `elapsed_sq_sum_s / samples - mean^2` suffers cancellation for long histories with small spread. No invariant ties `successes + failures` to `samples` inside the type; only the loader enforces `successes <= samples` and `failures <= samples - successes`. Mutation is only safe under the owning `OuterRouteState.lock`.

## Provenance
Mapped from `src/parallel/routing/outer_route_state.jl` line 60.
