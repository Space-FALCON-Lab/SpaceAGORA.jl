---
id: gnc.tracking_executor_control_solarpanels_heatload
label: control_solarpanels_heatload
kind: function
source:
  file: src/gnc/control/aerobraking/tracking_executor.jl
  symbol: control_solarpanels_heatload
  lines:
  - 192
  - 192
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
- id: gram_atmosphere
  type: Any
  units: n/a
  required: false
  description: Positional argument `gram_atmosphere` (default `nothing`).
- id: heat_rate_control
  type: Any
  units: n/a
  required: false
  description: Positional argument `heat_rate_control` (default `false`).
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
  description: Return value of `control_solarpanels_heatload`. Returns `_control_solarpanels_heatload_impl(ip,
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

# control_solarpanels_heatload

## Purpose
Public, lock-protected entry point for the heat-load bang-bang solar-panel controller used during aerobraking passes; it serialises access to the shared control bridge state before delegating to `_control_solarpanels_heatload_impl`.

## Design & Implementation
Signature `control_solarpanels_heatload(ip, m, args, index_ratio, state=0, t=0, position=0, current_position=0, gram_atmosphere=nothing, heat_rate_control=false; cnf=nothing)`. It calls `lock(CONTROL_BRIDGE_STATE_LOCK)`, then in a `try`/`finally` invokes the implementation with all arguments forwarded and guarantees `unlock` on both normal return and exception. The returned value is the commanded panel angle `α` in radians (0 or `m.aerodynamics.α`). The many positional defaults exist to keep the legacy call signature shared with `no_control`, `control_solarpanels_openloop` and the other controllers registered by name.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `ip` | Any | n/a | yes | Positional argument `ip`. |
| in | `m` | Any | n/a | yes | Positional argument `m`. |
| in | `args` | Any | n/a | yes | Positional argument `args`. |
| in | `index_ratio` | Any | n/a | yes | Positional argument `index_ratio`. |
| in | `state` | Any | n/a | no | Positional argument `state` (default `0`). |
| in | `t` | Any | n/a | no | Positional argument `t` (default `0`). |
| in | `position` | Any | n/a | no | Positional argument `position` (default `0`). |
| in | `current_position` | Any | n/a | no | Positional argument `current_position` (default `0`). |
| in | `gram_atmosphere` | Any | n/a | no | Positional argument `gram_atmosphere` (default `nothing`). |
| in | `heat_rate_control` | Any | n/a | no | Positional argument `heat_rate_control` (default `false`). |
| in | `cnf` | Any | n/a | no | Keyword argument `cnf` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `control_solarpanels_heatload`. Returns `_control_solarpanels_heatload_impl(ip, m, args, index_ratio, state, t, position,`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.tracking_executor__control_solarpanels_openloop_impl|_control_solarpanels_openloop_impl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:272-272`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/aerobraking/tracking_executor.jl`

**Downstream**

- `callees` → [[gnc.tracking_executor__control_solarpanels_heatload_impl|_control_solarpanels_heatload_impl]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:195-195`
<!-- vulcan:connections:end -->

## Limitations
The lock is held for the whole guidance dispatch and multibody rotation, so concurrent spacecraft serialise on it. Re-entrant use from `_control_solarpanels_openloop_impl` requires `CONTROL_BRIDGE_STATE_LOCK` to be re-entrant. Positional arguments `state`, `position` and `current_position` default to the integer 0 rather than typed placeholders, so type errors surface only inside the policy dispatch.

## Provenance
Mapped from `src/gnc/control/aerobraking/tracking_executor.jl` line 192.
