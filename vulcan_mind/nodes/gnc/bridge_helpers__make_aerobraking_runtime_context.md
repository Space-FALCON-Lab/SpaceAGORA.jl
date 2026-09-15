---
id: gnc.bridge_helpers__make_aerobraking_runtime_context
label: _make_aerobraking_runtime_context
kind: function
source:
  file: src/gnc/internal/bridge_helpers.jl
  symbol: _make_aerobraking_runtime_context
  lines:
  - 166
  - 166
inputs:
- id: mission
  type: Any
  units: n/a
  required: true
  description: Keyword argument `mission`.
- id: index_phase_aerobraking
  type: Any
  units: n/a
  required: true
  description: Keyword argument `index_phase_aerobraking`.
- id: ip
  type: Any
  units: n/a
  required: true
  description: Keyword argument `ip`.
- id: aerobraking_phase
  type: Any
  units: n/a
  required: true
  description: Keyword argument `aerobraking_phase`.
- id: t_prev
  type: Any
  units: n/a
  required: true
  description: Keyword argument `t_prev`.
- id: date_initial
  type: Any
  units: n/a
  required: true
  description: Keyword argument `date_initial`.
- id: time_0
  type: Any
  units: n/a
  required: true
  description: Keyword argument `time_0`.
- id: args
  type: Any
  units: n/a
  required: true
  description: Keyword argument `args`.
- id: initial_state
  type: Any
  units: n/a
  required: true
  description: Keyword argument `initial_state`.
- id: gram_atmosphere
  type: Any
  units: n/a
  required: true
  description: Keyword argument `gram_atmosphere`.
- id: gram
  type: Any
  units: n/a
  required: true
  description: Keyword argument `gram`.
- id: cnf
  type: Any
  units: n/a
  required: false
  description: Keyword argument `cnf` (default `nothing`).
- id: solution
  type: Any
  units: n/a
  required: false
  description: Keyword argument `solution` (default `nothing`).
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
  description: Return value of `_make_aerobraking_runtime_context`. Returns `(`.
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

# _make_aerobraking_runtime_context

## Purpose
Assembles the immutable NamedTuple context passed through the legacy aerobraking control bridge, bundling mission data, phase indices, timing, initial state, GRAM atmosphere handles, and a nested resolved `settings` tuple.

## Design & Implementation
A keyword-only constructor requiring `mission`, `index_phase_aerobraking`, `ip`, `aerobraking_phase`, `t_prev`, `date_initial`, `time_0`, `args`, `initial_state`, `gram_atmosphere`, and `gram`, with optional `cnf=nothing` and `solution=nothing`. All inputs are stored verbatim; the only computed field is `settings=_make_aerobraking_runtime_settings(args, mission)`, which eagerly resolves every `_bridge_aerobraking_*` accessor (topography, entry/exit interface in metres, body shape, heat load, max heat rate, SRP, control mode, struct control, dry mass, thrust phi, control-in-loop, integrator name, drag passage). Callers extend the tuple later with `_with_control_gain` and `_with_time_switch`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `mission` | Any | n/a | yes | Keyword argument `mission`. |
| in | `index_phase_aerobraking` | Any | n/a | yes | Keyword argument `index_phase_aerobraking`. |
| in | `ip` | Any | n/a | yes | Keyword argument `ip`. |
| in | `aerobraking_phase` | Any | n/a | yes | Keyword argument `aerobraking_phase`. |
| in | `t_prev` | Any | n/a | yes | Keyword argument `t_prev`. |
| in | `date_initial` | Any | n/a | yes | Keyword argument `date_initial`. |
| in | `time_0` | Any | n/a | yes | Keyword argument `time_0`. |
| in | `args` | Any | n/a | yes | Keyword argument `args`. |
| in | `initial_state` | Any | n/a | yes | Keyword argument `initial_state`. |
| in | `gram_atmosphere` | Any | n/a | yes | Keyword argument `gram_atmosphere`. |
| in | `gram` | Any | n/a | yes | Keyword argument `gram`. |
| in | `cnf` | Any | n/a | no | Keyword argument `cnf` (default `nothing`). |
| in | `solution` | Any | n/a | no | Keyword argument `solution` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_make_aerobraking_runtime_context`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.constraint_tracking_time_switch_func_affect_bang|time_switch_func_affect!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:353-353`
- [[gnc.control_commands_time_switch_func_affect_bang|time_switch_func_affect!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:365-365`
- [[gnc.eom_predictor_shooting_residual_bang|shooting_residual!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1014-1014`
- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:353-353`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:365-365`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1014-1014`

**Downstream**

- `callees` → [[gncz.bridge_helpers__make_aerobraking_runtime_settings|_make_aerobraking_runtime_settings]] · `callers` · call · `src/gnc/internal/bridge_helpers.jl:193-193`
<!-- vulcan:connections:end -->

## Limitations
Because `settings` is built eagerly, constructing the context throws `ArgumentError` if `EI` or `dry_mass` is missing even when the caller never uses them. Field types are unconstrained, so a wrong `initial_state` shape is not detected until dynamics run. The NamedTuple is copied on every `_with_*` extension, which allocates per call.

## Provenance
Mapped from `src/gnc/internal/bridge_helpers.jl` line 166.
