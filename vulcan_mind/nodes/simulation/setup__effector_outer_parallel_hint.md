---
id: simulation.setup__effector_outer_parallel_hint
label: _effector_outer_parallel_hint
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _effector_outer_parallel_hint
  lines:
  - 405
  - 405
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
  description: Return value of `_effector_outer_parallel_hint`.
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

# _effector_outer_parallel_hint

## Purpose
Reports whether an outer parallel layer (a sweep or Monte Carlo driver running several simulations concurrently) has flagged itself active, so the inner effector loop can avoid oversubscribing the thread pool.

## Design & Implementation
Thin forwarder returning `SimulationModel.ParallelPolicy.outer_parallel_active()`, which parses `SPACEAGORA_OUTER_PARALLEL_ACTIVE` as a boolean with default `false`. Used as the fallback in `_dynamic_effector_thread_decision` when no `PolicyDecisionEnvConfig` snapshot is supplied (`penv === nothing`).

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_effector_outer_parallel_hint`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__dynamic_effector_thread_decision|_dynamic_effector_thread_decision]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:667-667`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The hint is cooperative: an outer driver that forgets to set the variable leaves inner threading enabled and nested contention follows. Because the value is re-read from `ENV` on each fallback call, the decision can change mid-run if the driver toggles the variable.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 405.
