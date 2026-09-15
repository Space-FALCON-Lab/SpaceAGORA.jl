---
id: simulation.model_selection__ensure_gram_isolated_pool_bang
label: _ensure_gram_isolated_pool!
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/model_selection.jl
  symbol: _ensure_gram_isolated_pool!
  lines:
  - 130
  - 130
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: template_model
  type: EnvironmentModels.GRAMAtmosphereModel
  units: n/a
  required: true
  description: Positional argument `template_model`.
- id: workers
  type: Int
  units: n/a
  required: true
  description: Positional argument `workers`.
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
  type: Tuple{Vector{EnvironmentModels.GRAMAtmosphereModel},
  units: n/a
  description: Return value of `_ensure_gram_isolated_pool!`; mutates `p` in place.
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

# _ensure_gram_isolated_pool!

## Purpose
Guarantees that the shared buffers hold exactly `workers` independent deep copies of a GRAM atmosphere model plus one `ReentrantLock` per copy, so each worker thread in an isolated-pool batch evaluation can call the non-reentrant GRAM core without contention.

## Design & Implementation
Reads `p.shared_buffers.gram_isolated_pool_models` and `gram_isolated_pool_locks`. If `workers <= 0` the existing vectors are returned untouched. Otherwise `rebuild` is true when either vector length differs from `workers`; in that case both vectors are emptied, `sizehint!`ed to `workers`, and repopulated with `deepcopy(template_model)` and `ReentrantLock()` in a loop. The return type is `Tuple{Vector{GRAMAtmosphereModel}, Vector{ReentrantLock}}`. The pool is thus resized only on worker-count change, amortising the expensive deep copies across steps.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `template_model` | EnvironmentModels.GRAMAtmosphereModel | n/a | yes | Positional argument `template_model`. |
| in | `workers` | Int | n/a | yes | Positional argument `workers`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{Vector{EnvironmentModels.GRAMAtmosphereModel}, | n/a | — | Return value of `_ensure_gram_isolated_pool!`; mutates `p` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation_a.model_selection_gram_isolated_pool_batch_eval__gram_isolated_pool_batch_eval_bang|_gram_isolated_pool_batch_eval!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:200-200`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/model_selection.jl:148-148`
<!-- vulcan:connections:end -->

## Limitations
`deepcopy` of a GRAM model duplicates its native core state, which can be large; a pool rebuild inside a hot loop is very costly. The rebuild discards existing models even when shrinking by one worker, rather than truncating. If `template_model` changes identity between calls but `workers` stays the same, stale copies are kept. The function is not itself synchronised, so concurrent first-time callers could race on `empty!`/`push!`.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/model_selection.jl` line 130.
