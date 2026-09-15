---
id: parallel.types_adaptivecontrollerstate
label: AdaptiveControllerState
kind: struct
source:
  file: src/parallel/policy/types.jl
  symbol: AdaptiveControllerState
  lines:
  - 1
  - 1
inputs:
- id: desire
  type: Int64
  units: n/a
  required: false
  description: Field `desire` (default `1`).
- id: window_calls
  type: Int64
  units: n/a
  required: false
  description: Field `window_calls` (default `0`).
- id: window_allotment_sum
  type: Int64
  units: n/a
  required: false
  description: Field `window_allotment_sum` (default `0`).
- id: window_useful_sum
  type: Float64
  units: n/a
  required: false
  description: Field `window_useful_sum` (default `0.0`).
- id: window_deprived_calls
  type: Int64
  units: n/a
  required: false
  description: Field `window_deprived_calls` (default `0`).
- id: last_classification
  type: Symbol
  units: n/a
  required: false
  description: Field `last_classification` (default `:none`).
- id: last_utilization
  type: Float64
  units: n/a
  required: false
  description: Field `last_utilization` (default `1.0`).
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
  type: AdaptiveControllerState
  units: n/a
  description: Constructed `AdaptiveControllerState` (keyword constructor via @kwdef).
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

# AdaptiveControllerState

## Purpose
Per-workload-source state for the adaptive threading controller in `src/parallel/policy`. One instance lives per `source` symbol (for example `:control_callback` or `:density`) inside `PolicyContext.adaptive_state`, and carries the thread-count desire together with the measurement window used to classify recent parallel dispatches as efficient or inefficient.

## Design & Implementation
Declared with `Base.@kwdef mutable struct`, so `AdaptiveControllerState()` gives a serial default (`desire = 1`, `last_utilization = 1.0`, `last_classification = :none`). `desire::Int64` is the controller's requested worker count; `adaptive_decision.jl` clamps it to `[1, desire_cap]` and raises it via bootstrap or control-tail guards, and the dispatched allotment is `max(1, min(desire, budget))`. `window_calls`, `window_allotment_sum` and `window_useful_sum` accumulate over a window of `L` observations in `observation_tracking.jl`, where utilization is `window_useful_sum / max(1, window_allotment_sum)`; `window_deprived_calls` counts calls that received fewer threads than desired. After classification the window counters are reset to zero and `last_classification` becomes `:inefficient`, `:efficient_satisfied`, `:efficient_deprived` or `:measured_reward`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `desire` | Int64 | n/a | no | Field `desire` (default `1`). |
| in | `window_calls` | Int64 | n/a | no | Field `window_calls` (default `0`). |
| in | `window_allotment_sum` | Int64 | n/a | no | Field `window_allotment_sum` (default `0`). |
| in | `window_useful_sum` | Float64 | n/a | no | Field `window_useful_sum` (default `0.0`). |
| in | `window_deprived_calls` | Int64 | n/a | no | Field `window_deprived_calls` (default `0`). |
| in | `last_classification` | Symbol | n/a | no | Field `last_classification` (default `:none`). |
| in | `last_utilization` | Float64 | n/a | no | Field `last_utilization` (default `1.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AdaptiveControllerState | n/a | — | Constructed `AdaptiveControllerState` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/types.jl`
- [[parallel.policy_telemetry__adaptive_state_for|_adaptive_state_for]] · `callees` → `callers` · call · `src/parallel/policy/policy_telemetry.jl:4-4`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The struct has no lock of its own; all mutation must happen under `_policy_telemetry_lock`, which the type cannot enforce. Fields are plain `Int64`/`Float64`, so a window with `window_allotment_sum = 0` relies on the caller's `max(1, ...)` guard to avoid division by zero. `last_utilization` defaults to `1.0` before any observation, which reads as fully efficient at cold start.

## Provenance
Mapped from `src/parallel/policy/types.jl` line 1.
