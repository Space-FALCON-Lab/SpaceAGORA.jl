---
id: simulation.setup__validate_orientation_inertia_bang
label: _validate_orientation_inertia!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _validate_orientation_inertia!
  lines:
  - 35
  - 35
inputs:
- id: args
  type: Any
  units: n/a
  required: true
  description: Positional argument `args`.
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
  description: Return value of `_validate_orientation_inertia!`; mutates `args` in
    place. Returns `nothing`.
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

# _validate_orientation_inertia!

## Purpose
Guarantees that every spacecraft's inertia tensor is usable by the attitude integrator before a run with `orientation_sim=true` begins, failing fast instead of producing NaN quaternions mid-simulation.

## Design & Implementation
Returns `nothing` immediately when `args.mission_configuration.orientation_sim` is false. Otherwise it iterates `enumerate(args.dynamics_model.spacecraft)`, converts each `sc.inertia_tensor` to a dense `Matrix`, and applies three checks in order: `all(isfinite, ...)`, `issymmetric(...)`, and `isposdef(Symmetric(...))`. Each failure throws an `ArgumentError` naming the spacecraft index and the violated property. Despite the `!` suffix nothing is mutated; the name signals a validation pass with side effects limited to throwing.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_validate_orientation_inertia!`; mutates `args` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:164-164`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`issymmetric` uses exact equality, so tensors assembled with floating-point round-off (for example a rotated tensor) can be rejected even though they are physically symmetric; no tolerance is applied. `isposdef` performs a Cholesky factorisation, which is fine for 3×3 but is repeated for every spacecraft on every setup call. Zero-length spacecraft lists pass trivially.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 35.
