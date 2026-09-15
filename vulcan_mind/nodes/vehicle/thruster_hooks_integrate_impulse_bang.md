---
id: vehicle.thruster_hooks_integrate_impulse_bang
label: integrate_impulse!
kind: function
source:
  file: src/vehicle/actuators/thruster/thruster_hooks.jl
  symbol: integrate_impulse!
  lines:
  - 89
  - 89
inputs:
- id: link
  type: Link
  units: n/a
  required: true
  description: Positional argument `link`.
- id: thruster
  type: Thruster
  units: n/a
  required: true
  description: Positional argument `thruster`.
- id: on_time_request
  type: Float64
  units: n/a
  required: true
  description: Positional argument `on_time_request`.
- id: time
  type: Float64
  units: n/a
  required: true
  description: Positional argument `time`.
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
  description: Return value of `integrate_impulse!`; mutates `link` in place. Returns
    `total_integrated_thrust`.
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

# integrate_impulse!

## Purpose

Integrates the delivered impulse of one `Thruster` over a single attitude control period, accounting for the first-order thrust build-up and decay of the valve, and mutates `thruster.κ` to carry the end-of-period thrust fraction into the next call. It returns the total integrated thrust in newton-seconds accumulated over `link.attitude_control_rate` seconds.

## Design & Implementation

The requested on-time is clamped with `clamp(on_time_request, 0.0, link.attitude_control_rate)`. The cutoff frequency is taken as `abs(thruster.cutoff_frequency)` and floored at `1e-9` rad/s when non-finite or smaller, preventing division by zero. `expm1` is used for the exponential terms so small products stay accurate. If the on-time is shorter than the control period, the remaining `ramp_down_dt` contributes a tail-off impulse and the thrust fraction is further decayed. The `time` argument is accepted but unused here.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `link` | Link | n/a | yes | Positional argument `link`. |
| in | `thruster` | Thruster | n/a | yes | Positional argument `thruster`. |
| in | `on_time_request` | Float64 | n/a | yes | Positional argument `on_time_request`. |
| in | `time` | Float64 | n/a | yes | Positional argument `time`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `integrate_impulse!`; mutates `link` in place. Returns `total_integrated_thrust`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl`
- [[vehicle.thruster_hooks_thrust_calculation_schmitt_trigger_bang|thrust_calculation_schmitt_trigger!]] · `callees` → `callers` · call · `src/vehicle/actuators/thruster/thruster_hooks.jl:71-71`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Theory & Math

The valve thrust fraction obeys a first-order lag toward unity while commanded on, $\dot{\kappa} = \omega(1-\kappa)$, and toward zero while off, $\dot{\kappa} = -\omega\kappa$, with $\omega$ the cutoff frequency in rad/s and $\kappa$ the dimensionless thrust fraction. Integrating over the on-time $t_i$ with initial fraction $\kappa_0$ gives the impulse

$$ I_{\text{on}} = F_{\max}\left(t_i + \frac{\kappa_0 - 1}{\omega}\left(1 - e^{-\omega t_i}\right)\right) $$

and end-of-pulse fraction $\kappa_1 = 1 + (\kappa_0-1)e^{-\omega t_i}$. Over the remaining ramp-down interval $t_d = \Delta t_c - t_i$ the tail impulse is

$$ I_{\text{off}} = \frac{F_{\max}\kappa_1}{\omega}\left(1 - e^{-\omega t_d}\right), \qquad \kappa_2 = \kappa_1 e^{-\omega t_d} $$

where $F_{\max}$ is `thruster.max_thrust` in newtons and $\Delta t_c$ is `link.attitude_control_rate` in seconds.

## Limitations

The model assumes a single first-order pole with the same time constant for rise and fall, and one on-pulse placed at the start of each control period; multiple pulses per period cannot be represented. Clamping the on-time to the control period silently discards any longer request. The `1e-9` rad/s floor makes a nominally instantaneous valve behave as an extremely slow one rather than raising an error, and the thrust fraction is never bounded to the unit interval if a caller seeds it outside that range.

## Provenance
Mapped from `src/vehicle/actuators/thruster/thruster_hooks.jl` line 89.
