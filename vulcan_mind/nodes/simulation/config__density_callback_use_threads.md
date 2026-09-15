---
id: simulation.config__density_callback_use_threads
label: _density_callback_use_threads
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _density_callback_use_threads
  lines:
  - 288
  - 288
inputs:
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
- id: num_sats
  type: Int
  units: n/a
  required: true
  description: Positional argument `num_sats`.
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
  description: Return value of `_density_callback_use_threads`.
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

# _density_callback_use_threads

## Purpose
Boolean convenience accessor answering only whether the density callback should run threaded for a given satellite count.

## Design & Implementation
Calls `_density_callback_thread_decision(args, num_sats)` and returns the `use_threads` field of the resulting named tuple, discarding `allotment`, `mode` and `policy_applied`. Provided for call sites and tests that need the yes/no answer without the worker budget.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_density_callback_use_threads`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/config.jl`

**Downstream**

- `callees` → [[simulation.config__density_callback_thread_decision|_density_callback_thread_decision]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:289-289`
<!-- vulcan:connections:end -->

## Limitations
Discarding `allotment` means a caller that threads on the strength of this answer has no information about how many workers it is entitled to and may oversubscribe. Each call re-runs the whole decision, including a possible `ENV` re-parse when the parameter object carries no snapshot.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 288.
