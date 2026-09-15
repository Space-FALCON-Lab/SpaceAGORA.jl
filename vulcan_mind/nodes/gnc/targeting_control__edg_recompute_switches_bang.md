---
id: gnc.targeting_control__edg_recompute_switches_bang
label: _edg_recompute_switches!
kind: function
source:
  file: src/gnc/control/targeting_control.jl
  symbol: _edg_recompute_switches!
  lines:
  - 96
  - 96
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
- id: pos
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos`.
- id: vel
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel`.
- id: mass
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mass`.
- id: heat_load_j_cm2
  type: Float64
  units: n/a
  required: true
  description: Positional argument `heat_load_j_cm2`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  type: Nothing
  units: n/a
  description: Return value of `_edg_recompute_switches!`; mutates `model` in place.
    Returns `nothing`.
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

# _edg_recompute_switches!

## Purpose
Solves the mode-specific switch times once per drag passage — the targeting switch, or the heat-load low-drag window — and records when they were solved.

## Design & Implementation
In targeting mode with targeting active, it returns if a finite switch is already stored or the vehicle is not yet in the passage, otherwise solves `_edg_solve_targeting_switch` with heat-rate and structural sub-modes taken from the config and stores the result. In max-energy-depletion mode with the heat-load sub-mode, it clears the solved flag on leaving the passage, sets the passage-active flag on entry, and solves `_edg_solve_heat_load_switches` once per passage. Mutates several `state` vectors at index `i`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | AerobrakingEnergyDepletionControlModel | n/a | yes | Positional argument `model`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `env` | Any | n/a | yes | Positional argument `env`. |
| in | `spacecraft` | Any | n/a | yes | Positional argument `spacecraft`. |
| in | `pos` | SVector{3, Float64} | n/a | yes | Positional argument `pos`. |
| in | `vel` | SVector{3, Float64} | n/a | yes | Positional argument `vel`. |
| in | `mass` | Float64 | n/a | yes | Positional argument `mass`. |
| in | `heat_load_j_cm2` | Float64 | n/a | yes | Positional argument `heat_load_j_cm2`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `i` | Int | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_edg_recompute_switches!`; mutates `model` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.targeting_control_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:264-264`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/targeting_control.jl`

**Downstream**

- `callees` → [[gnc.targeting_control__edg_in_drag_passage|_edg_in_drag_passage]] · `callers` · call · `src/gnc/control/targeting_control.jl:112-112`
- `callees` → [[gnc.targeting_control__edg_solve_targeting_switch|_edg_solve_targeting_switch]] · `callers` · call · `src/gnc/control/targeting_control.jl:113-113`
- `callees` → [[gncx.heat_load_control__edg_solve_heat_load_switches|_edg_solve_heat_load_switches]] · `callers` · call · `src/gnc/control/targeting_control.jl:138-138`
<!-- vulcan:connections:end -->

## Limitations
The targeting switch is solved once and never revised within a passage, so a density mis-prediction early in the pass is not corrected; the heat-load flags reset only on exit, so a passage that ends the simulation leaves them set.

## Provenance
Mapped from `src/gnc/control/targeting_control.jl` line 96.
