---
id: gnc.momentum_manager_calccontroleffect_bang
label: calcControlEffect!
kind: function
source:
  file: src/gnc/control/momentum_manager.jl
  symbol: calcControlEffect!
  lines:
  - 63
  - 63
inputs:
- id: model
  type: MagneticMomentumManagerModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  type: Nothing
  units: n/a
  description: Return value of `calcControlEffect!`; mutates `model` in place. Returns
    `nothing`.
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

# calcControlEffect!

## Purpose
Advances the magnetic momentum manager's internal wheel-momentum estimate each control tick, so desaturation can be commanded against a current momentum state.

## Design & Implementation
Returns immediately unless the satellite index matches the model's own, and again if either the commanded-torque callback or the cached body-frame magnetic field is absent. It extracts position, velocity, attitude quaternion and body rate from the satellite's state, increments `ticks`, and evaluates the commanded torque callback. On the first call it seeds `h_wheels` from `h_wheels_0`; afterwards it integrates wheel momentum backwards by the commanded torque over the elapsed `dt`, guarded so that only finite positive steps are applied.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | MagneticMomentumManagerModel | n/a | yes | Positional argument `model`. |
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `calcControlEffect!`; mutates `model` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.propulsive_maneuvers_calccontrolforcetorque|calcControlForceTorque]] · `callees` → `callers` · call · `src/gnc/control/propulsive_maneuvers.jl:425-425`
- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:37-37`
- [[simulation.control_callbacks__schedule_thruster_control_bang|_schedule_thruster_control!]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:37-37`
- [[simulation_a.control_callbacks_get_control_callbacks|get_control_callbacks]] · `callees` → `callers` · call · `src/simulation/callbacks/control_callbacks.jl:100-100`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:377-377`

**Downstream**

- `callees` → [[core.quaternion_utils_rot|rot]] · `callers` · call · `src/gnc/control/momentum_manager.jl:86-86`
<!-- vulcan:connections:end -->

## Limitations
Momentum is integrated with a first-order rectangular step keyed off wall-clock differences between calls, so an irregular or repeated tick sequence biases the estimate; the guard rejects non-positive steps rather than correcting them.

## Provenance
Mapped from `src/gnc/control/momentum_manager.jl` line 63.
