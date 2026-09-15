---
id: gnc.robot_arm_control__robot_arm_control_spacecraft_state
label: _robot_arm_control_spacecraft_state
kind: function
source:
  file: src/gnc/control/robot_arm_control.jl
  symbol: _robot_arm_control_spacecraft_state
  lines:
  - 178
  - 178
inputs:
- id: u
  type: Any
  units: n/a
  required: true
  description: Positional argument `u`.
- id: idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `idx`.
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
  description: 'Return value of `_robot_arm_control_spacecraft_state`. Returns `hasproperty(u,
    :sc) && length(u.sc) >= idx ? u.sc[idx] : nothing`.'
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

# _robot_arm_control_spacecraft_state

## Purpose
Safely selects the per-spacecraft sub-view `u.sc[idx]` from the full simulation state so the arm controller can read link attitudes and rates, returning `nothing` when the state does not contain that spacecraft.

## Design & Implementation
Signature `_robot_arm_control_spacecraft_state(u, idx::Int)`. It returns `u.sc[idx]` when `hasproperty(u, :sc)` and `length(u.sc) >= idx`, otherwise `nothing`. With a `ComponentVector` state, `u.sc[idx]` is a view into the underlying storage, so no copy is made. `calcControlEffect!` passes the result to `robot_arm_measured_joint_state`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `u` | Any | n/a | yes | Positional argument `u`. |
| in | `idx` | Int | n/a | yes | Positional argument `idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_control_spacecraft_state`. Returns `hasproperty(u, :sc) && length(u.sc) >= idx ? u.sc[idx] : nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.robot_arm_control_calccontroleffect_bang|calcControlEffect!]] · `callees` → `callers` · call · `src/gnc/control/robot_arm_control.jl:190-190`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/robot_arm_control.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`idx <= 0` passes the length check and then throws a `BoundsError` from indexing rather than returning `nothing`. Plain `Vector{Float64}` states (as used by some solver modes that flatten the state) have no `sc` property, so the controller silently falls back to open-loop reference tracking in those modes. The function assumes `u.sc` supports `length` and integer indexing.

## Provenance
Mapped from `src/gnc/control/robot_arm_control.jl` line 178.
