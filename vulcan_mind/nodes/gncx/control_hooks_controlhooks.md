---
id: gncx.control_hooks_controlhooks
label: ControlHooks
kind: struct
source:
  file: src/gnc/control/control_hooks.jl
  symbol: ControlHooks
  lines:
  - 1
  - 53
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Parent GNC namespace from which ControlHooks imports structure, configuration,
    environment, ephemeris, and guidance types.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: control_api
  type: Module
  units: n/a
  description: 'Exported control-effector interface: force/torque, mass-flow, reaction-wheel,
    and MPC constructors for every control model.'
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncx
origin: agent
---
# ControlHooks

## Purpose
`ControlHooks` is the module that assembles the entire GNC control layer into one namespace and defines the interface the dynamics right-hand side calls. It is the single boundary between the simulation core and every control effector: aerobraking energy depletion, propulsive maneuvers, rendezvous-and-proximity-operations model predictive control, robot-arm joint control, and magnetic momentum management.

## Model & Assumptions
The module fixes a four-function control contract, exporting `calcControlForceTorque`, `calcControlEffect!`, `calcControlMassFlowRate`, and `calcReactionWheelTorque`. Each concrete effector supplies methods for that contract, and a fallback defined in the propulsive-maneuver file covers effectors that do not consume propellant or drive wheels. Dependencies are imported explicitly by name rather than by blanket `using`, which makes the coupling to configuration types, abstract types, thruster models, guidance models, gravity and aerodynamic effectors, environment and ephemeris models, reference systems, and kinematics visible at the top of the file.

## Design & Implementation
Implementation is organised as an ordered include list. Shared bridge helpers and quaternion utilities load first because later files use them; propulsive maneuvers load next to establish the fallback methods; then the aerobraking chain of heat-rate, heat-load, structural, and targeting control; then the RPO subdirectory in dependency order with types before allocators before the model; then robot-arm and momentum-manager effectors; and finally the aerobraking command, constraint-tracking, and tracking-executor files. `SparseArrays` and `OSQP` are pulled in at module scope because the LQ-MPC controllers set up a sparse quadratic program at construction time. `const config = Structure` gives the aerobraking code a short alias for the spacecraft structure queries.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Parent GNC namespace from which ControlHooks imports structure, configuration, environment, ephemeris, and guidance types. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `control_api` | Module | n/a | — | Exported control-effector interface: force/torque, mass-flow, reaction-wheel, and MPC constructors for every control model. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/control/control_hooks.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Include order is load-bearing and not enforced by any mechanism other than the ordering in this file, so moving an include can break method resolution silently. Every effector shares one flat namespace, so exported symbol names must stay globally unique across all control models. The module takes a hard dependency on OSQP even for simulations that never instantiate an MPC controller.

## Provenance
Mapped from `src/gnc/control/control_hooks.jl:1-53`.
