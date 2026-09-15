---
id: simulation.setup__dynamic_effectors_parallel_supported
label: _dynamic_effectors_parallel_supported
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _dynamic_effectors_parallel_supported
  lines:
  - 487
  - 487
inputs:
- id: dynamic_effectors
  type: Tuple
  units: n/a
  required: true
  description: Positional argument `dynamic_effectors`.
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
  description: Return value of `_dynamic_effectors_parallel_supported`.
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

# _dynamic_effectors_parallel_supported

## Purpose
Decides whether an entire effector tuple can be evaluated with inner threading, requiring every element to be thread-safe and at most one `AerodynamicCoefficientfM` because that model shares density workspace buffers.

## Design & Implementation
Iterates `dynamic_effectors::Tuple` with `@inbounds`, counting `aero_fm_count` for elements that `isa SimulationModel.AerodynamicCoefficientfM` and returning `false` immediately if `_dynamic_effector_threadsafe(effector)` is false. After the loop returns `aero_fm_count <= 1`. Pure, allocation-free, and consulted early in `_dynamic_effector_thread_decision` before any cost estimation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_dynamic_effectors_parallel_supported`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__dynamic_effector_thread_decision|_dynamic_effector_thread_decision]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:652-652`
- [[simulation.setup__rhs_flat_supported|_rhs_flat_supported]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:926-926`

**Downstream**

- `callees` → [[simulation.setup__dynamic_effector_threadsafe|_dynamic_effector_threadsafe]] · `callers` · call · `src/simulation/engine/setup.jl:493-493`
<!-- vulcan:connections:end -->

## Limitations
The one-aero rule is hard-coded rather than derived from whether the two aero models actually share buffers. An empty tuple returns `true` even though there is nothing to parallelise; callers must check `n_effectors` separately (and do). Because the loop short-circuits on the first unsafe effector, the aero count may be incomplete, which is harmless.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 487.
