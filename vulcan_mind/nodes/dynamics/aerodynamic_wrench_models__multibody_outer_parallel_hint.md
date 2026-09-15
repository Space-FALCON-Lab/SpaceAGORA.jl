---
id: dynamics.aerodynamic_wrench_models__multibody_outer_parallel_hint
label: _multibody_outer_parallel_hint
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: _multibody_outer_parallel_hint
  lines:
  - 26
  - 26
inputs:
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
  description: Return value of `_multibody_outer_parallel_hint`.
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

# _multibody_outer_parallel_hint

## Purpose
Reports whether an outer (per-satellite) parallel region is already active so the per-link aerodynamic loop can avoid nested threading.

## Design & Implementation
Returns `ParallelPolicy.outer_parallel_active()::Bool`. `_multibody_thread_decision` passes it as `outer_active`, and unless `SPACEAGORA_MULTIBODY_PARALLEL_ALLOW_WITH_OUTER` is true the policy disables inner threads when this is set.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_multibody_outer_parallel_hint`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models__multibody_thread_decision|_multibody_thread_decision]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:38-38`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Relies on `ParallelPolicy` correctly tracking outer-region entry and exit; a leaked flag would permanently disable inner threading for the process. It is a hint only and does not itself prevent oversubscription.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 26.
