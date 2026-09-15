---
id: simulation.config__callback_outer_parallel_hint
label: _callback_outer_parallel_hint
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _callback_outer_parallel_hint
  lines:
  - 131
  - 131
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
  description: Return value of `_callback_outer_parallel_hint`.
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

# _callback_outer_parallel_hint

## Purpose
Reports whether an outer parallel construct is currently active, so callback-level threading decisions can avoid nesting threads inside an already-parallel driver.

## Design & Implementation
One-line `@inline` forwarder to `ParallelPolicy.outer_parallel_active()`, returning a `Bool`. It is used as the fallback source of `outer_active` in all three thread-decision functions whenever the ODE parameter object carries no `PolicyDecisionEnvConfig` snapshot.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_callback_outer_parallel_hint`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/config.jl`
- [[simulation.config__control_callback_thread_decision|_control_callback_thread_decision]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:305-305`
- [[simulation.config__density_callback_thread_decision|_density_callback_thread_decision]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:260-260`
- [[simulation.config__thermal_callback_thread_decision|_thermal_callback_thread_decision]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:336-336`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The name is accurate: it is a hint, derived from whatever state the parallel policy module tracks, and it cannot report how many workers the outer construct holds or whether that construct is currently idle. A driver that parallelises without registering itself with `ParallelPolicy` will not be detected, and nesting will proceed.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 131.
