---
id: vehicle.thruster_hooks_thrusterhooks
label: ThrusterHooks
kind: module
source:
  file: src/vehicle/actuators/thruster/thruster_hooks.jl
  symbol: ThrusterHooks
  lines:
  - 1
  - 1
inputs:
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
  description: Value produced by this symbol.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# ThrusterHooks

## Purpose

`ThrusterHooks` is the module that turns a commanded body torque into per-thruster on-times and average thrusts for a `Link`. It exports `update_thrusters!`, `thrust_calculation_schmitt_trigger!`, `schmitt_trigger` and `integrate_impulse!`, and depends on `LinearAlgebra` for the pseudoinverse and cross products, plus `CSV` and `DataFrames` for optional debug logging.

## Design & Implementation

The module pulls `Link`, `Thruster` and `rotate_to_body` from `SpacecraftModels`, `Components` and `Kinematics`. The entry point `update_thrusters!` rebuilds the 3xN torque Jacobian `link.J_thruster` from each thruster's body-frame location and unit direction, solves `pinv(J) * torque` for a thrust vector, shifts and clamps it to non-negative values, then hands each element to the Schmitt-trigger pulse model. Debug output is gated behind `thruster_debug_enabled()`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

Thrust allocation uses a minimum-norm pseudoinverse with no thrust-magnitude bounds, so a commanded torque larger than the cluster can produce is silently clipped by the per-thruster saturation in `thrust_calculation_schmitt_trigger!`. The non-negativity shift changes the realized torque whenever the raw solution has a negative entry. All thrusters in a `Link` share a single `attitude_control_rate` control period.

## Provenance
Mapped from `src/vehicle/actuators/thruster/thruster_hooks.jl` line 1.
