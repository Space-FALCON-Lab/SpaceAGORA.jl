---
id: gnc.tracking_executor__control_solarpanels_openloop_impl
label: _control_solarpanels_openloop_impl
kind: function
source:
  file: src/gnc/control/aerobraking/tracking_executor.jl
  symbol: _control_solarpanels_openloop_impl
  lines:
  - 269
  - 269
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
  description: Return value of `_control_solarpanels_openloop_impl`. Returns `α`.
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

# _control_solarpanels_openloop_impl

## Purpose
Unlocked implementation of the open-loop combined controller: it first runs the heat-load scheduler to refresh the switch times, then inside or outside the switch window either commands a fully feathered panel (angle 0) or delegates to the heat-rate root-finding controller.

## Design & Implementation
Invoked by `control_solarpanels_openloop` under `CONTROL_BRIDGE_STATE_LOCK`. It reads `cnf_state` and `heat_load_sol`, then calls `control_solarpanels_heatload(ip, m, args, index_ratio, 0, t, position, current_position, gram_atmosphere, heat_rate_control; cnf=cnf_state)`, which recomputes `cnf_state.time_switch_1/2` and rotates the links. With `heat_load_sol` 0 or 3, `t` within `[time_switch_1, time_switch_2]` gives `α = 0`, otherwise `α = control_solarpanels_heatrate(...)`; with 1 or 2 the branches are swapped. The heat-rate controller receives `state` (`T_p`, `ρ`, `S`) and uses `cnf_state.α_past` as its Newton seed. The angle in radians is returned; the link rotation done by the heat-load call is not repeated here (that block is commented out).

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
| out | `result` | Any | n/a | — | Return value of `_control_solarpanels_openloop_impl`. Returns `α`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.tracking_executor_control_solarpanels_openloop|control_solarpanels_openloop]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:263-263`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/aerobraking/tracking_executor.jl`

**Downstream**

- `callees` → [[gnc.bridge_helpers__bridge_get_cnf|_bridge_get_cnf]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:270-270`
- `callees` → [[gnc.bridge_helpers__bridge_required_field|_bridge_required_field]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:271-271`
- `callees` → [[gnc.tracking_executor_control_solarpanels_heatload|control_solarpanels_heatload]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:272-272`
- `callees` → [[gncx.tracking_executor_control_solarpanels_heatrate|control_solarpanels_heatrate]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:278-278`
<!-- vulcan:connections:end -->

## Limitations
Because `control_solarpanels_heatload` acquires `CONTROL_BRIDGE_STATE_LOCK` again, correctness depends on that lock being a `ReentrantLock`; with a plain lock this path would deadlock. The window test uses inclusive `>=`/`<=` whereas the heat-load implementation uses strict inequalities, so at exactly the switch instants the two controllers disagree. `heat_load_sol` outside 0-3 leaves `α` undefined and throws. The heat-rate branch ignores the panel rotation applied by the heat-load call, so the model's link angles reflect the bang-bang command rather than the returned `α`.

## Provenance
Mapped from `src/gnc/control/aerobraking/tracking_executor.jl` line 269.
