---
id: simulation.model_selection__gram_batch_elapsed_time
label: _gram_batch_elapsed_time
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/model_selection.jl
  symbol: _gram_batch_elapsed_time
  lines:
  - 78
  - 78
inputs:
- id: el_time
  type: Float64
  units: n/a
  required: true
  description: Positional argument `el_time`.
- id: idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `idx`.
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
  type: Float64
  units: n/a
  description: Return value of `_gram_batch_elapsed_time`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _gram_batch_elapsed_time

## Purpose
Normalises the elapsed-time argument of a batched GRAM density evaluation, which may be a single scalar shared by all satellites or a per-satellite vector, into a `Float64` for satellite index `idx`.

## Design & Implementation
Two `@inline` methods with an explicit `::Float64` return annotation. The `(el_time::Float64, idx::Int)` method ignores `idx` and returns `el_time` unchanged. The `(el_time::AbstractVector{<:Real}, idx::Int)` method returns `Float64(el_time[idx])`, converting integer or `Float32` entries. Dispatch happens at compile time on the argument type, so the branch has no runtime cost inside the `threaded_foreach_worker_persistent` loop in `_gram_isolated_pool_batch_eval!`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `el_time` | Float64 | n/a | yes | Positional argument `el_time`. |
| in | `idx` | Int | n/a | yes | Positional argument `idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_gram_batch_elapsed_time`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.model_selection_gram_isolated_pool_batch_eval__gram_isolated_pool_batch_eval_bang|_gram_isolated_pool_batch_eval!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:207-207`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:83-83`
<!-- vulcan:connections:end -->

## Limitations
The vector method performs a bounds-checked index; callers must guarantee `length(el_time) == n`, which `_gram_isolated_pool_batch_eval!` verifies beforehand but other callers may not. No method exists for other time representations such as `DateTime` or `Dates.Second`. Units are assumed to be seconds elapsed since the initial epoch, consistent with GRAM's elapsed-time input.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/model_selection.jl` line 78.
