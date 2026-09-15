---
id: gnc.tracking_executor__control_solarpanels_heatload_impl
label: _control_solarpanels_heatload_impl
kind: function
source:
  file: src/gnc/control/aerobraking/tracking_executor.jl
  symbol: _control_solarpanels_heatload_impl
  lines:
  - 201
  - 201
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
  description: Return value of `_control_solarpanels_heatload_impl`. Returns `α`.
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

# _control_solarpanels_heatload_impl

## Purpose
Unlocked implementation of the heat-load solar-panel scheduler: it asks the aerobraking guidance policy for the two switch times of the current pass, converts them into a bang-bang panel angle command, shifts the heat-load history buffers, and physically rotates the panel links in the multibody model.

## Design & Implementation
Called only from `control_solarpanels_heatload` while `CONTROL_BRIDGE_STATE_LOCK` is held. It fetches the shared `cnf_state` via `_bridge_get_cnf(args; cnf)` and the integer `heat_load_sol` via `_bridge_required_field(args, :heat_load_sol)`, builds an `AerobrakingGuidanceInput` (normalising `index_ratio` to `Vector{Int}` and `t` to `Float64`) and calls `dispatch_aerobraking_guidance(DefaultAerobrakingPolicySelector(), AerobrakingPolicyConfig(), input)`. The returned `time_switch_1`, `time_switch_2` and `security_mode` are written into `cnf_state`. For `heat_load_sol` 0 or 3 the angle is 0 inside the switch window and `m.aerodynamics.α` outside; for 1 or 2 the polarity is reversed. It then lazily allocates `cnf_state.heat_load_ppast` to `length(m.body.links)` zeros and copies `heat_load_past` into it. Finally it traverses `m.body` from `roots[1]` with `config.traverse_bodies` and calls `config.rotate_link(body, abs.(body.r), -α + root.α)` for every non-root link.

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
| out | `result` | Any | n/a | — | Return value of `_control_solarpanels_heatload_impl`. Returns `α`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.tracking_executor_control_solarpanels_heatload|control_solarpanels_heatload]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:195-195`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/aerobraking/tracking_executor.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:212-212`
- `callees` → [[gnc.bridge_helpers__bridge_get_cnf|_bridge_get_cnf]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:202-202`
- `callees` → [[gnc.bridge_helpers__bridge_required_field|_bridge_required_field]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:203-203`
- `callees` → [[gnc.interfaces_aerobrakingguidanceinput|AerobrakingGuidanceInput]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:206-206`
- `callees` → [[gncx.dispatcher_dispatch_aerobraking_guidance|dispatch_aerobraking_guidance]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:219-219`
- `callees` → [[misc.policy_types_aerobrakingpolicyconfig|AerobrakingPolicyConfig]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:205-205`
- `callees` → [[mission.selector_stub_defaultaerobrakingpolicyselector|DefaultAerobrakingPolicySelector]] · `callers` · call · `src/gnc/control/aerobraking/tracking_executor.jl:204-204`
<!-- vulcan:connections:end -->

## Limitations
If `heat_load_sol` is any value other than 0-3, `α` is never assigned and the function throws an `UndefVarError` at the rotation loop. A fresh `DefaultAerobrakingPolicySelector` and `AerobrakingPolicyConfig` are constructed on every call. The panel geometry assumes the standard two-panel-one-bus layout with rotation axis `abs.(body.r)`, and the rotation is applied cumulatively to the live model (mutation of `m.body`) every call. `α` is an `Int` 0 in one branch and a `Float64` in the other.

## Provenance
Mapped from `src/gnc/control/aerobraking/tracking_executor.jl` line 201.
