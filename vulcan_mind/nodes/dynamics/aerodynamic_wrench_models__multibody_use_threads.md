---
id: dynamics.aerodynamic_wrench_models__multibody_use_threads
label: _multibody_use_threads
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _multibody_use_threads
  lines:
  - 31
  - 31
inputs:
- id: num_items
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_items`.
- id: heavy_work
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `heavy_work` (default `true`).
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
  type: Bool
  units: n/a
  description: Return value of `_multibody_use_threads`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# _multibody_use_threads

## Purpose
Boolean convenience over `_multibody_thread_decision` answering whether `num_items` links should be processed with threads.

## Design & Implementation
Returns `_multibody_thread_decision(num_items; heavy_work=heavy_work).use_threads`. `heavy_work` defaults to `true`, matching the free-molecular coefficient cost per link.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_items` | Int | n/a | yes | Positional argument `num_items`. |
| in | `heavy_work` | Bool | n/a | no | Keyword argument `heavy_work` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_multibody_use_threads`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- `callees` → [[dynamics.aerodynamic_wrench_models__multibody_thread_decision|_multibody_thread_decision]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:32-32`
<!-- vulcan:connections:end -->

## Limitations
Discards the `allotment` and `mode` fields, so callers that later need the worker count must call the full decision function again and repeat all environment parsing.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 31.
