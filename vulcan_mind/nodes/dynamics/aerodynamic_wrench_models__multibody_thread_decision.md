---
id: dynamics.aerodynamic_wrench_models__multibody_thread_decision
label: _multibody_thread_decision
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _multibody_thread_decision
  lines:
  - 35
  - 35
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
  type: Tuple
  units: n/a
  description: 'Return value of `_multibody_thread_decision`. Returns `(use_threads=use_threads,
    allotment=use_threads ? allotment : 1, mode=mode)`.'
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

# _multibody_thread_decision

## Purpose
Combines environment-driven policy inputs into a single decision NamedTuple `(use_threads, allotment, mode)` governing intra-satellite threading of the per-link aero loop.

## Design & Implementation
Gathers `mode`, `threshold`, `outer_active`, `allow_with_outer` (`SPACEAGORA_MULTIBODY_PARALLEL_ALLOW_WITH_OUTER`, default false), and `heavy_only` (`SPACEAGORA_MULTIBODY_PARALLEL_HEAVY_ONLY`, default true), then calls `ParallelPolicy.thread_policy_decision(num_items; mode, threshold, heavy_work, heavy_only, outer_active, allow_with_outer, source=:multibody)`. The allotment is capped by `_multibody_max_threads()`, and `use_threads` additionally requires `allotment > 1`. When threads are not used the allotment is reported as `1`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `num_items` | Int | n/a | yes | Positional argument `num_items`. |
| in | `heavy_work` | Bool | n/a | no | Keyword argument `heavy_work` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `_multibody_thread_decision`. Returns `(use_threads=use_threads, allotment=use_threads ? allotment : 1, mode=mode)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models__multibody_use_threads|_multibody_use_threads]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:32-32`
- [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:773-773`
- [[dynamics.perturbations_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:925-925`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- `callees` → [[dynamics.aerodynamic_wrench_models__multibody_max_threads|_multibody_max_threads]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:51-51`
- `callees` → [[dynamics.aerodynamic_wrench_models__multibody_outer_parallel_hint|_multibody_outer_parallel_hint]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:38-38`
- `callees` → [[dynamics.aerodynamic_wrench_models__multibody_parallel_mode|_multibody_parallel_mode]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:36-36`
- `callees` → [[dynamics.aerodynamic_wrench_models__multibody_thread_threshold|_multibody_thread_threshold]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:37-37`
- `callees` → [[parcore.adaptive_decision_thread_policy_decision|thread_policy_decision]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:41-41`
- `callees` → [[simulation.config__parse_bool_env|_parse_bool_env]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:39-39`
<!-- vulcan:connections:end -->

## Limitations
Five environment lookups per call inside the RHS hot path. The decision is recomputed for every satellite on every step even though inputs only change when the environment changes. No memoisation of the policy result.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 35.
