---
id: parallel.thread_execution_threaded_collect_persistent_bang
label: threaded_collect_persistent!
kind: function
source:
  file: src/parallel/policy/thread_execution.jl
  symbol: threaded_collect_persistent!
  lines:
  - 260
  - 260
inputs:
- id: source
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `source`.
- id: results
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `results`.
- id: num_items
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_items`.
- id: allotment
  type: Int
  units: n/a
  required: true
  description: Positional argument `allotment`.
- id: f
  type: F
  units: n/a
  required: true
  description: Positional argument `f`.
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
  type: Nothing
  units: n/a
  description: 'Return value of `threaded_collect_persistent!`; mutates `source` in
    place. Returns `nothing` or `results`. Type parameters: `{F <: Function}`.'
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

# threaded_collect_persistent!

## Purpose
Persistent-pool counterpart of `threaded_collect!`: evaluates `f(idx)` on long-lived worker tasks keyed by `source` and stores results in `results[idx]`, preserving deterministic caller-side reduction while avoiding per-call task spawn cost.

## Design & Implementation
Wraps `threaded_foreach_worker_persistent(source, num_items, allotment) do _, idx ... end` with the body `@inbounds results[idx] = f(idx)`. The `source::Symbol` selects the pool via `_persistent_worker_pool_for`; if persistent workers are disabled by `SPACEAGORA_CALLBACK_PERSISTENT_WORKERS=0` or `workers <= 1`, the underlying primitive transparently falls back to `threaded_foreach_worker`. Returns `results`. A `do`-block ordering method is provided.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `source` | Symbol | n/a | yes | Positional argument `source`. |
| in | `results` | AbstractVector | n/a | yes | Positional argument `results`. |
| in | `num_items` | Int | n/a | yes | Positional argument `num_items`. |
| in | `allotment` | Int | n/a | yes | Positional argument `allotment`. |
| in | `f` | F | n/a | yes | Positional argument `f`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `threaded_collect_persistent!`; mutates `source` in place. Returns `nothing` or `results`. Type parameters: `{F <: Function}`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:950-950`
- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/thread_execution.jl`
- [[parallel.thread_execution_threaded_foreach_worker|threaded_foreach_worker]] · `callees` → `callers` · feedback · `src/parallel/policy/thread_execution.jl:242-242`

**Downstream**

- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/parallel/policy/thread_execution.jl:262-262`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/parallel/policy/thread_execution.jl:262-262`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/parallel/policy/thread_execution.jl:262-262`
- `callees` → [[parallel.thread_execution_threaded_collect_bang|threaded_collect!]] · `callers` · call · `src/parallel/policy/thread_execution.jl:269-269`
- `callees` → [[parallel.thread_execution_threaded_foreach_worker_persistent|threaded_foreach_worker_persistent]] · `callers` · call · `src/parallel/policy/thread_execution.jl:261-261`
<!-- vulcan:connections:end -->

## Limitations
Same sizing and element-type requirements as `threaded_collect!` apply, with `@inbounds` writes. The pool's `run_lock` serialises batches per `source`, so two concurrent callers with the same `source` in the same policy scope block each other. A closure type that differs per call site defeats specialisation inside the persistent worker loop and incurs dynamic dispatch on each item.

## Provenance
Mapped from `src/parallel/policy/thread_execution.jl` line 260.
