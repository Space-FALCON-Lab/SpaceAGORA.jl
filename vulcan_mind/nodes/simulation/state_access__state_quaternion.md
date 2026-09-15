---
id: simulation.state_access__state_quaternion
label: _state_quaternion
kind: function
source:
  file: src/simulation/engine/state_access.jl
  symbol: _state_quaternion
  lines:
  - 83
  - 83
inputs:
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
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
  type: Union{Nothing, SVector}
  units: n/a
  description: Return value of `_state_quaternion`. Returns `nothing` or `SVector{4,
    Float64}(u.sc[sat_idx].q)`.
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

# _state_quaternion

## Purpose
Retrieves the attitude quaternion of spacecraft `sat_idx` as an `SVector{4, Float64}`, or `nothing` when the active state layout does not propagate attitude.

## Design & Implementation
Returns `nothing` immediately if `_is_gravity_backbone_state(u)` holds, since the second-order translational backbone carries no attitude block, or if `u.sc[sat_idx]` has no `:q` property. Otherwise it wraps `u.sc[sat_idx].q` into a statically sized four-element vector, converting to `Float64`. Returning `nothing` rather than an identity quaternion lets callers distinguish absent attitude from a genuinely unrotated body.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, SVector} | n/a | — | Return value of `_state_quaternion`. Returns `nothing` or `SVector{4, Float64}(u.sc[sat_idx].q)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/state_access.jl`
- [[simulation.state_access__state_has_quaternion|_state_has_quaternion]] · `callees` → `callers` · call · `src/simulation/engine/state_access.jl:91-91`

**Downstream**

- `callees` → [[simulation.state_access__is_gravity_backbone_state|_is_gravity_backbone_state]] · `callers` · call · `src/simulation/engine/state_access.jl:84-84`
<!-- vulcan:connections:end -->

## Limitations
No normalisation or sign-convention check is performed, so a drifting integrated quaternion is returned unnormalised and the scalar-first versus scalar-last ordering is an unstated contract with the attitude dynamics. The `hasproperty` probe touches `u.sc[sat_idx]` before the backbone guard can protect it in the non-backbone branch, so a malformed flat state throws instead of returning `nothing`.

## Provenance
Mapped from `src/simulation/engine/state_access.jl` line 83.
