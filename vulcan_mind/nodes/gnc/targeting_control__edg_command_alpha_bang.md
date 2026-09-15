---
id: gnc.targeting_control__edg_command_alpha_bang
label: _edg_command_alpha!
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_command_alpha!
  lines:
  - 173
  - 173
inputs:
- id: model
  type: AerobrakingEnergyDepletionControlModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: env
  type: Any
  units: n/a
  required: true
  description: Positional argument `env`.
- id: spacecraft
  type: Any
  units: n/a
  required: true
  description: Positional argument `spacecraft`.
- id: base_alpha
  type: Float64
  units: n/a
  required: true
  description: Positional argument `base_alpha`.
- id: heat_load_j_cm2
  type: Float64
  units: n/a
  required: true
  description: Positional argument `heat_load_j_cm2`.
- id: heat_load_low_drag_active
  type: Bool
  units: n/a
  required: true
  description: Positional argument `heat_load_low_drag_active`.
- id: i
  type: Int
  units: n/a
  required: true
  description: Positional argument `i`.
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
  description: Return value of `_edg_command_alpha!`; mutates `model` in place. Returns
    `alpha`.
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

# _edg_command_alpha!

## Purpose
Applies the heat-rate, structural-load and heat-load constraints to the base angle and records the resulting command and diagnostics in the state.

## Design & Implementation
Clamps the base angle, then, unless the heat-load low-drag window is active, computes `alpha_hr` through the Maxwellian heat-rate root solve warm-started from the last command, and `alpha_struct` through the structural-load solve over the controlled links. The command is the minimum of whichever constraints are active. If the heat-load sub-mode is on and the accumulated load has reached its limit, the command drops to the minimum angle. It then writes the final, heat-rate and structural angles, the resulting heat rate, the heat load and the dynamic pressure into the state vectors at `i`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | AerobrakingEnergyDepletionControlModel | n/a | yes | Positional argument `model`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `env` | Any | n/a | yes | Positional argument `env`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `base_alpha` | Float64 | n/a | yes | Positional argument `base_alpha`. |
| in | `heat_load_j_cm2` | Float64 | n/a | yes | Positional argument `heat_load_j_cm2`. |
| in | `heat_load_low_drag_active` | Bool | n/a | yes | Positional argument `heat_load_low_drag_active`. |
| in | `i` | Int | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_command_alpha!`; mutates `model` in place. Returns `alpha`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:267-267`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[gnc.heat_rate_control__edg_maxwellian_heat_rate|_edg_maxwellian_heat_rate]] · `callers` · call · `src/gnc/control/targeting_control.jl:223-223`
- `callees` → [[gncx.heat_rate_control__edg_heat_rate_alpha|_edg_heat_rate_alpha]] · `callers` · call · `src/gnc/control/targeting_control.jl:194-194`
- `callees` → [[gncx.struct_load_control__edg_structural_alpha|_edg_structural_alpha]] · `callers` · call · `src/gnc/control/targeting_control.jl:199-199`
<!-- vulcan:connections:end -->

## Limitations
The two constraint solves each cost several root iterations per tick; and since the minimum of the two is taken, a conservative structural model can mask the thermal constraint entirely, which the diagnostics reveal but nothing flags.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 173.
