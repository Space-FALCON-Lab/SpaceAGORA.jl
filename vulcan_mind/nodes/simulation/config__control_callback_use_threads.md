---
id: simulation.config__control_callback_use_threads
label: _control_callback_use_threads
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/config.jl
  symbol: _control_callback_use_threads
  lines:
  - 324
  - 324
inputs:
- id: control_model
  type: Any
  units: n/a
  required: true
  description: Positional argument `control_model`.
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
  description: Return value of `_control_callback_use_threads`.
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

# _control_callback_use_threads

## Purpose
Boolean accessor answering only whether the control callback should run threaded for a given controller and satellite count.

## Design & Implementation
Calls `_control_callback_thread_decision(control_model, num_sats)` and returns the `use_threads` field, discarding the worker allotment, the requested mode, and the `policy_applied` flag that distinguishes a thread-safety veto from a policy-driven serial decision.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `control_model` | Any | n/a | yes | Positional argument `control_model`. |
| in | `num_sats` | Int | n/a | yes | Positional argument `num_sats`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_control_callback_use_threads`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/config.jl`

**Downstream**

- `callees` → [[simulation.config__control_callback_thread_decision|_control_callback_thread_decision]] · `callers` · call · `src/simulation/callbacks/density_callbacks/config.jl:325-325`
<!-- vulcan:connections:end -->

## Limitations
Because `policy_applied` is dropped, a caller cannot tell whether a `false` answer came from the model being unsafe or from the satellite count sitting below the threshold, which makes misconfiguration hard to diagnose. The full decision is recomputed on each call.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/config.jl` line 324.
