---
id: simulation.setup__rhs_flat_batch_privileged_effector
label: _rhs_flat_batch_privileged_effector
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _rhs_flat_batch_privileged_effector
  lines:
  - 732
  - 732
inputs:
- id: effector
  type: Any
  units: n/a
  required: true
  description: Positional argument `effector`.
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
  description: Return value of `_rhs_flat_batch_privileged_effector`.
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

# _rhs_flat_batch_privileged_effector

## Purpose
Identifies effectors that gain an algorithmic benefit from the flat-constellation route independent of effector count, so a two-effector harmonics-plus-drag constellation is not locked out of the flat queue by the minimum-effector gate.

## Design & Implementation
Returns `true` for `GravitationalHarmonicsModel` (receives the SIMD batch Pines kernel across the whole satellite batch), and for `InverseSquaredGravityModel` or `InverseSquaredJ2GravityModel` only when `!effector.gravity_gradient` (those take the batched pre-pass `_accumulate_invsq_flat_batch!` / `_accumulate_invsq_j2_flat_batch!`). Everything else, including N-body and SRP whose ephemeris lookups are already hoisted by `_prefill_shared_body_samples!`, returns `false`. `@inline` and pure.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `effector` | Any | n/a | yes | Positional argument `effector`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_rhs_flat_batch_privileged_effector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__rhs_flat_has_batch_privileged_effector|_rhs_flat_has_batch_privileged_effector]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:741-741`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The set of privileged types is hard-coded in this function and must be updated whenever `dynamics_rhs.jl` adds a new batched kernel. Gravity-gradient variants are excluded because their flat kernel does not exist, not because batching would be incorrect; the exclusion is a code-coverage limitation.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 732.
