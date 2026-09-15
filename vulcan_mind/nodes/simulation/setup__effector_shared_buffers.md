---
id: simulation.setup__effector_shared_buffers
label: _effector_shared_buffers
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _effector_shared_buffers
  lines:
  - 506
  - 506
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
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
  type: Nothing
  units: n/a
  description: Return value of `_effector_shared_buffers`. Returns `nothing` or `getproperty(p,
    :shared_buffers)`.
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

# _effector_shared_buffers

## Purpose
Safely extracts the `shared_buffers` object from an ODE parameter bundle, returning `nothing` when the bundle is absent or hand-built without that field so cost-model code can degrade gracefully in tests.

## Design & Implementation
`_effector_shared_buffers(p)` returns `nothing` if `p === nothing` or `!hasproperty(p, :shared_buffers)`, otherwise `getproperty(p, :shared_buffers)`. Using `hasproperty` rather than a type check allows named tuples and mock structs in unit tests. `@inline`, no allocation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_effector_shared_buffers`. Returns `nothing` or `getproperty(p, :shared_buffers)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__dynamic_effector_thread_decision|_dynamic_effector_thread_decision]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:678-678`
- [[simulation.setup__policy_env_config|_policy_env_config]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:901-901`
- [[simulation.setup__rhs_env_config_from_buffers|_rhs_env_config_from_buffers]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:898-898`
- [[simulation.setup__rhs_execution_plan_uncached|_rhs_execution_plan_uncached]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1231-1231`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A `p` whose `shared_buffers` field is itself `nothing` is passed through as `nothing`, which downstream code handles, but a field of an unexpected type is passed through unchecked and fails later with `hasproperty` returning false on individual cost fields.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 506.
