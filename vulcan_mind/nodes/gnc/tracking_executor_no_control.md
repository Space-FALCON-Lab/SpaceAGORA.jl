---
id: gnc.tracking_executor_no_control
label: no_control
kind: function
source:
  file: src/gnc/control/aerobraking/tracking_executor.jl
  symbol: no_control
  lines:
  - 25
  - 25
inputs:
- id: ip
  type: Any
  units: n/a
  required: true
  description: Positional argument `ip`.
- id: m
  type: Any
  units: n/a
  required: true
  description: Positional argument `m`.
- id: args
  type: Any
  units: n/a
  required: false
  description: Positional argument `args` (default `0`).
- id: index_ratio
  type: Any
  units: n/a
  required: false
  description: Positional argument `index_ratio` (default `0`).
- id: state
  type: Any
  units: n/a
  required: false
  description: Positional argument `state` (default `0`).
- id: t
  type: Any
  units: n/a
  required: false
  description: Positional argument `t` (default `0`).
- id: position
  type: Any
  units: n/a
  required: false
  description: Positional argument `position` (default `0`).
- id: current_position
  type: Any
  units: n/a
  required: false
  description: Positional argument `current_position` (default `0`).
- id: heat_rate_control
  type: Any
  units: n/a
  required: false
  description: Positional argument `heat_rate_control` (default `true`).
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
  type: Any
  units: n/a
  description: Return value of `no_control`. Returns `α`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# no_control

## Purpose
Null controller that returns the configured fixed solar-panel angle without any feedback, used when a scenario disables active panel control but the control hook still has to produce an angle each step.

## Design & Implementation
Signature `no_control(ip, m, args=0, index_ratio=0, state=0, t=0, position=0, current_position=0, heat_rate_control=true)`; every argument after `m` has a default so it can be called with the same positional convention as the active controllers. It reads `m.aerodynamics.α` (rad) and returns it unchanged. Nothing is mutated and no lock is taken.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `ip` | Any | n/a | yes | Positional argument `ip`. |
| in | `m` | Any | n/a | yes | Positional argument `m`. |
| in | `args` | Any | n/a | no | Positional argument `args` (default `0`). |
| in | `index_ratio` | Any | n/a | no | Positional argument `index_ratio` (default `0`). |
| in | `state` | Any | n/a | no | Positional argument `state` (default `0`). |
| in | `t` | Any | n/a | no | Positional argument `t` (default `0`). |
| in | `position` | Any | n/a | no | Positional argument `position` (default `0`). |
| in | `current_position` | Any | n/a | no | Positional argument `current_position` (default `0`). |
| in | `heat_rate_control` | Any | n/a | no | Positional argument `heat_rate_control` (default `true`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `no_control`. Returns `α`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/aerobraking/tracking_executor.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The trailing arguments are accepted purely for signature compatibility and are ignored, so misuse is not detected. The returned angle is whatever the mission configuration holds; there is no clamp to `[0, π/2]`. Because it bypasses the control bridge state, `cnf_state.α_past` is never updated, which can affect a later switch to a Newton-seeded controller.

## Provenance
Mapped from `src/gnc/control/aerobraking/tracking_executor.jl` line 25.
