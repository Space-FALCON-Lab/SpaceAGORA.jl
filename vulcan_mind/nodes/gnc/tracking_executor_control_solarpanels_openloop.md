---
id: gnc.tracking_executor_control_solarpanels_openloop
label: control_solarpanels_openloop
kind: function
source:
  file: src/gnc/control/aerobraking/tracking_executor.jl
  symbol: control_solarpanels_openloop
  lines:
  - 260
  - 260
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
  required: true
  description: Positional argument `args`.
- id: index_ratio
  type: Any
  units: n/a
  required: true
  description: Positional argument `index_ratio`.
- id: state
  type: Any
  units: n/a
  required: true
  description: Positional argument `state`.
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
- id: gram_atmosphere
  type: Any
  units: n/a
  required: false
  description: Positional argument `gram_atmosphere` (default `nothing`).
- id: cnf
  type: Any
  units: n/a
  required: false
  description: Keyword argument `cnf` (default `nothing`).
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
  description: Return value of `control_solarpanels_openloop`. Returns `_control_solarpanels_openloop_impl(ip,
    m, args, index_ratio, state, t, position,`.
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

# control_solarpanels_openloop

## Purpose
Public, lock-protected entry point for the open-loop combined controller that mixes the heat-load switch schedule with the heat-rate limiting root solve, returning the panel angle for the current integration step.

## Design & Implementation
Signature `control_solarpanels_openloop(ip, m, args, index_ratio, state, t=0, position=0, current_position=0, heat_rate_control=true, gram_atmosphere=nothing; cnf=nothing)`. Note the positional order differs from `control_solarpanels_heatload`: here `heat_rate_control` precedes `gram_atmosphere`. It acquires `CONTROL_BRIDGE_STATE_LOCK`, calls `_control_solarpanels_openloop_impl` inside `try`/`finally`, and unlocks. `state` is required (no default) because the heat-rate solve needs `T_p`, `ρ` and `S` from it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `ip` | Any | n/a | yes | Positional argument `ip`. |
| in | `m` | Any | n/a | yes | Positional argument `m`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `index_ratio` | Any | n/a | yes | Positional argument `index_ratio`. |
| in | `state` | Any | n/a | yes | Positional argument `state`. |
| in | `t` | Any | n/a | no | Positional argument `t` (default `0`). |
| in | `position` | Any | n/a | no | Positional argument `position` (default `0`). |
| in | `current_position` | Any | n/a | no | Positional argument `current_position` (default `0`). |
| in | `heat_rate_control` | Any | n/a | no | Positional argument `heat_rate_control` (default `true`). |
| in | `gram_atmosphere` | Any | n/a | no | Positional argument `gram_atmosphere` (default `nothing`). |
| in | `cnf` | Any | n/a | no | Keyword argument `cnf` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `control_solarpanels_openloop`. Returns `_control_solarpanels_openloop_impl(ip, m, args, index_ratio, state, t, position,`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/aerobraking/tracking_executor.jl`

**Downstream**

- `callees` → [[gnc.tracking_executor__control_solarpanels_openloop_impl|_control_solarpanels_openloop_impl]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:263-263`
<!-- vulcan:connections:end -->

## Limitations
The differing positional order between this function and `control_solarpanels_heatload` is an easy source of argument mix-ups since both accept untyped arguments. The implementation re-acquires the same lock through `control_solarpanels_heatload`, so the lock must be re-entrant. Nothing validates `state` length before `control_solarpanels_heatrate` indexes `state[1:3]`.

## Provenance
Mapped from `src/gnc/control/aerobraking/tracking_executor.jl` line 260.
