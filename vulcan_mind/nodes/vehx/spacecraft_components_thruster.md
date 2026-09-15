---
id: vehx.spacecraft_components_thruster
label: Thruster
kind: struct
source:
  file: src/vehicle/spacecraft/components.jl
  symbol: Thruster
  lines:
  - 19
  - 31
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Vehicle package namespace that loads this file and exports the symbol.
- id: geometry
  type: MVector{3,Float64}
  units: m
  required: true
  description: Nozzle location and thrust direction expressed in the link frame.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: thruster
  type: Thruster
  units: n/a
  description: Mutable nozzle record holding hardware limits, trigger thresholds and
    live firing state.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
- components
- thruster
charts:
- vehx
origin: agent
---

# Thruster

## Purpose
`Thruster` is the hardware record for a single reaction control nozzle attached to a link. It carries both the invariant description of the device and the mutable state that the pulse scheduler updates during a run, so the allocation hook, the impulse integrator and the dynamics effectors all read and write the same object rather than passing parallel arrays around.

## Theory & Math
The Schmitt trigger converts a normalised demand $u = \kappa f / f_{max}$ into a valve state that switches on at $u > \lambda_{on}$ and off at $u < \lambda_{off}$, so the hysteresis band $\lambda_{on} - \lambda_{off}$ sets the limit-cycle amplitude. Delivered thrust follows a first-order lag $\dot f = \omega_c (f_{cmd} - f)$ with cutoff $\omega_c$.

## Model & Assumptions
Geometry is given by a location relative to the link centre of mass and a direction unit vector, both in link axes. Performance is bounded by a maximum thrust and characterised by a specific impulse and a first-order cutoff frequency that shapes the ramp between commanded and delivered thrust. Pulse modulation is described by two dimensionless Schmitt trigger thresholds together with a minimum firing time, which prevents commands shorter than the valve can physically execute. A thrust factor completes the actuation model. The current thrust magnitude and the scheduled stop time are the two live state fields.

## Design & Implementation
The type is declared `@kwdef mutable struct` so every field has a documented default and a scenario need only override what differs from the nominal nozzle. Defaults are deliberately conservative: unit maximum thrust, zero location and direction, thresholds of 0.75 and 0.25, and zero current thrust. Location and direction use `MVector{3,Float64}` because `update_thrusters!` normalises the direction in place, which requires mutability while retaining stack allocation. Every field carries a trailing comment giving its physical units, which is the unit contract the rest of the actuator chain relies on.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Vehicle package namespace that loads this file and exports the symbol. |
| in | `geometry` | MVector{3,Float64} | m | yes | Nozzle location and thrust direction expressed in the link frame. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `thruster` | Thruster | n/a | — | Mutable nozzle record holding hardware limits, trigger thresholds and live firing state. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/components.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Nothing enforces that the direction is a unit vector at construction, that `level_off` is below `level_on`, or that the thrust never exceeds `max_thrust` once the allocator has written it. Propellant consumption is not tracked on the nozzle, plume impingement and thermal soak-back are absent, and mutability means a thruster shared between two links would alias state.

## Provenance
Mapped from `src/vehicle/spacecraft/components.jl:19-31`, alongside the `Facet`, `Magnet` and `ReactionWheelAssembly` records in the same module.
