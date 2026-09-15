---
id: simulation.setup__rhs_single_invsq_flat_supported
label: _rhs_single_invsq_flat_supported
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_single_invsq_flat_supported
  lines:
  - 953
  - 953
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
  description: Return value of `_rhs_single_invsq_flat_supported`.
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

# _rhs_single_invsq_flat_supported

## Purpose
Tests whether the effector set is exactly one inverse-square or J2 gravity model without gravity gradient, which can be evaluated in the flat queue without scratch buffers.

## Design & Implementation
Requires one effector of either inverse-square type, `gravity_gradient` false, and thread safety. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `dynamic_effectors` | Tuple | n/a | yes | Positional argument `dynamic_effectors`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_rhs_single_invsq_flat_supported`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__rhs_execution_plan_uncached|_rhs_execution_plan_uncached]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1176-1176`
- [[simulation.setup__rhs_flat_supported|_rhs_flat_supported]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:924-924`

**Downstream**

- `callees` → [[simulation.setup__dynamic_effector_threadsafe|_dynamic_effector_threadsafe]] · `callers` · call · `src/simulation/engine/setup.jl:959-959`
<!-- vulcan:connections:end -->

## Limitations
Gravity-gradient torque requires attitude state that the flat path does not stage, hence the exclusion; the comment in source notes there is no coefficient table to share, unlike harmonics.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 953.
