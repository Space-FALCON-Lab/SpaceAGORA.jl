---
id: parallel.types__hintlayerstatsaccumulator
label: _HintLayerStatsAccumulator
kind: struct
source:
  file: src/parallel/policy/types.jl
  symbol: _HintLayerStatsAccumulator
  lines:
  - 76
  - 76
inputs:
- id: signatures
  type: Set{String}
  units: n/a
  required: false
  description: Field `signatures` (default `Set{String}()`).
- id: choice_count
  type: Int64
  units: n/a
  required: false
  description: Field `choice_count` (default `0`).
- id: samples_total
  type: Int64
  units: n/a
  required: false
  description: Field `samples_total` (default `0`).
- id: successes_total
  type: Int64
  units: n/a
  required: false
  description: Field `successes_total` (default `0`).
- id: failures_total
  type: Int64
  units: n/a
  required: false
  description: Field `failures_total` (default `0`).
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
- id: confidence_sum
  type: Float64
  units: n/a
  required: false
  description: Field `confidence_sum` (default `0.0`).
- id: regret_sum_ns
  type: Float64
  units: n/a
  required: false
  description: Field `regret_sum_ns` (default `0.0`).
- id: signature_metric_count
  type: Int64
  units: n/a
  required: false
  description: Field `signature_metric_count` (default `0`).
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
  type: _HintLayerStatsAccumulator
  units: n/a
  description: Constructed `_HintLayerStatsAccumulator` (keyword constructor via @kwdef).
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

# _HintLayerStatsAccumulator

## Purpose
Aggregation scratch struct used when summarising the persistent-hint history for telemetry. It rolls every `AdaptiveChoiceStats` in one hint layer (a group of signatures) into totals so the reporting code can emit mean elapsed time, variance, average confidence and average regret per layer.

## Design & Implementation
`Base.@kwdef mutable struct` with an empty `Set{String}` of `signatures` and zeroed counters. `persistent_hints.jl` (around lines 440 to 466) pushes each signature into `signatures`, increments `choice_count` per allotment entry, adds `samples_total`, `successes_total`, `failures_total`, `elapsed_sum_ns` and `elapsed_sq_sum_ns`, and, once per signature, adds `confidence_sum` and `regret_sum_ns = max(0, best_mean - observed_best_mean)` while incrementing `signature_metric_count`. Layer variance is `max(0, elapsed_sq_sum_ns / n - mean^2)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `signatures` | Set{String} | n/a | no | Field `signatures` (default `Set{String}()`). |
| in | `choice_count` | Int64 | n/a | no | Field `choice_count` (default `0`). |
| in | `samples_total` | Int64 | n/a | no | Field `samples_total` (default `0`). |
| in | `successes_total` | Int64 | n/a | no | Field `successes_total` (default `0`). |
| in | `failures_total` | Int64 | n/a | no | Field `failures_total` (default `0`). |
| in | `elapsed_sum_ns` | Float64 | n/a | no | Field `elapsed_sum_ns` (default `0.0`). |
| in | `elapsed_sq_sum_ns` | Float64 | n/a | no | Field `elapsed_sq_sum_ns` (default `0.0`). |
| in | `confidence_sum` | Float64 | n/a | no | Field `confidence_sum` (default `0.0`). |
| in | `regret_sum_ns` | Float64 | n/a | no | Field `regret_sum_ns` (default `0.0`). |
| in | `signature_metric_count` | Int64 | n/a | no | Field `signature_metric_count` (default `0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | _HintLayerStatsAccumulator | n/a | — | Constructed `_HintLayerStatsAccumulator` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/types.jl`
- [[parallel.persistent_hints_hint_layer_stats_snapshot|hint_layer_stats_snapshot]] · `callees` → `callers` · call · `src/parallel/policy/persistent_hints.jl:426-426`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The struct is a transient accumulator and is never locked; it must be built entirely inside the caller's `_persistent_hint_lock` region. `Set{String}` allocation means aggregation costs O(number of signatures) memory. Regret is clamped to be non-negative so the layer average can never reveal a case where the observed best beat the predicted best. Division by `signature_metric_count` is the caller's responsibility and is not guarded here.

## Provenance
Mapped from `src/parallel/policy/types.jl` line 76.
