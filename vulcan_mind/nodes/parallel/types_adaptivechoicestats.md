---
id: parallel.types_adaptivechoicestats
label: AdaptiveChoiceStats
kind: struct
source:
  file: src/parallel/policy/types.jl
  symbol: AdaptiveChoiceStats
  lines:
  - 68
  - 68
inputs:
- id: samples
  type: Int64
  units: n/a
  required: false
  description: Field `samples` (default `0`).
- id: successes
  type: Int64
  units: n/a
  required: false
  description: Field `successes` (default `0`).
- id: failures
  type: Int64
  units: n/a
  required: false
  description: Field `failures` (default `0`).
- id: elapsed_sum_ns
  type: Float64
  units: n/a
  required: false
  description: Field `elapsed_sum_ns` (default `0.0`).
- id: elapsed_sq_sum_ns
  type: Float64
  units: n/a
  required: false
  description: Field `elapsed_sq_sum_ns` (default `0.0`).
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
  type: AdaptiveChoiceStats
  units: n/a
  description: Constructed `AdaptiveChoiceStats` (keyword constructor via @kwdef).
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

# AdaptiveChoiceStats

## Purpose
Running statistics for one thread-allotment choice under one workload signature in the persistent-hint layer. It records how many times an allotment was tried, how often it beat the alternative (`successes`/`failures`) and the first two moments of elapsed wall time so a mean and variance can be recovered without storing samples.

## Design & Implementation
`Base.@kwdef mutable struct` with all-zero defaults. `samples::Int64` counts observations, `elapsed_sum_ns::Float64` and `elapsed_sq_sum_ns::Float64` accumulate `elapsed` and `elapsed^2` in nanoseconds (`persistent_hints.jl` line 365). Instances are stored as `Dict{String, Dict{Int64, AdaptiveChoiceStats}}` in `_PersistentHintState.history`, keyed by signature then allotment. Serialization to JSON writes each field with `max(0.0, ...)` clamping, and merging on load adds the sums field-wise so histories from multiple runs combine. Mean is `elapsed_sum_ns / samples` and variance `elapsed_sq_sum_ns / samples - mean^2`, clamped at zero.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `samples` | Int64 | n/a | no | Field `samples` (default `0`). |
| in | `successes` | Int64 | n/a | no | Field `successes` (default `0`). |
| in | `failures` | Int64 | n/a | no | Field `failures` (default `0`). |
| in | `elapsed_sum_ns` | Float64 | n/a | no | Field `elapsed_sum_ns` (default `0.0`). |
| in | `elapsed_sq_sum_ns` | Float64 | n/a | no | Field `elapsed_sq_sum_ns` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AdaptiveChoiceStats | n/a | — | Constructed `AdaptiveChoiceStats` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/types.jl`
- [[parallel.persistent_hints__hint_payload_stats|_hint_payload_stats]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:49-49`
- [[parallel.persistent_hints__hint_record_observation_bang|_hint_record_observation!]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:359-359`
- [[parallel.persistent_hints__load_persistent_hint_state_locked_bang|_load_persistent_hint_state_locked!]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:100-100`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Storing sums of squares as `Float64` loses precision when `elapsed_sq_sum_ns` grows large relative to the variance (catastrophic cancellation in `E[x^2] - E[x]^2`), which can report zero variance for long-running histories. There is no decay or windowing, so stale measurements from a different machine weigh as much as fresh ones once merged. `successes + failures` is not enforced to equal `samples`.

## Provenance
Mapped from `src/parallel/policy/types.jl` line 68.
