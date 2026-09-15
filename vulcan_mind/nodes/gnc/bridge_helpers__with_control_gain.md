---
id: gnc.bridge_helpers__with_control_gain
label: _with_control_gain
kind: function
source:
  file: src/gnc/internal/bridge_helpers.jl
  symbol: _with_control_gain
  lines:
  - 201
  - 201
inputs:
- id: context
  type: NamedTuple
  units: n/a
  required: true
  description: Positional argument `context`.
- id: control_gain
  type: Any
  units: n/a
  required: true
  description: Positional argument `control_gain`.
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
  type: Tuple
  units: n/a
  description: Return value of `_with_control_gain`. Returns `(; context..., control_gain)`.
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

# _with_control_gain

## Purpose
Returns a copy of an aerobraking runtime `context::NamedTuple` extended with (or overriding) a `control_gain` field, so guidance code can carry the active gain without mutating shared state.

## Design & Implementation
Implemented as the one-line splat `(; context..., control_gain)`. Julia NamedTuple construction with a later duplicate key overrides the earlier one, so calling this on a context that already has `control_gain` replaces the value. The function is `@inline` and the result is a new immutable NamedTuple whose type depends on the field set and the type of `control_gain`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `context` | NamedTuple | n/a | yes | Positional argument `context`. |
| in | `control_gain` | Any | n/a | yes | Positional argument `control_gain`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `_with_control_gain`. Returns `(; context..., control_gain)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.constraint_tracking_time_switch_func_affect_bang|time_switch_func_affect!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:352-352`
- [[gnc.control_commands_time_switch_func_affect_bang|time_switch_func_affect!]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:364-364`
- [[gnc.eom_predictor_shooting_residual_bang|shooting_residual!]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1013-1013`
- [[gncx.constraint_tracking_asim_ctrl_plot|asim_ctrl_plot]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/constraint_tracking.jl:352-352`
- [[gncx.control_commands_asim_ctrl|asim_ctrl]] · `callees` → `callers` · call · `src/gnc/control/aerobraking/control_commands.jl:364-364`
- [[gncy.eom_predictor_asim_ctrl_targeting_plot|asim_ctrl_targeting_plot]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/t_edg/eom_predictor.jl:1013-1013`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Each call allocates a new NamedTuple containing every field of `context`, and a changing type of `control_gain` across calls produces distinct NamedTuple types, which can trigger dynamic dispatch in callers. No validation of `control_gain` (sign, magnitude, units) is performed.

## Provenance
Mapped from `src/gnc/internal/bridge_helpers.jl` line 201.
