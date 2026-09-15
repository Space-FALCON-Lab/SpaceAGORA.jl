---
id: simulation.dynamics_rhs__prefill_atmosphere_samples_bang
label: _prefill_atmosphere_samples!
kind: function
source:
  file: src/simulation/engine/dynamics_rhs.jl
  symbol: _prefill_atmosphere_samples!
  lines:
  - 1291
  - 1291
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: sc_state
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_state`.
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
  description: Return value of `_prefill_atmosphere_samples!`; mutates `p` in place.
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

# _prefill_atmosphere_samples!

## Purpose
Convenience wrapper that prefills both planet-frame and atmosphere samples for all satellites before a flat-queue RHS call.

## Design & Implementation
Calls `_prefill_environment_samples!(p, t, sc_state; atmosphere=true)` and returns `nothing`. Exists so call sites read as intent rather than as a keyword flag.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `sc_state` | Any | n/a | yes | Positional argument `sc_state`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_prefill_atmosphere_samples!`; mutates `p` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/dynamics_rhs.jl`

**Downstream**

- `callees` → [[simulation.dynamics_rhs__prefill_environment_samples_bang|_prefill_environment_samples!]] · `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1292-1292`
<!-- vulcan:connections:end -->

## Limitations
None beyond those of the underlying prefill.

## Provenance
Mapped from `src/simulation/engine/dynamics_rhs.jl` line 1291.
