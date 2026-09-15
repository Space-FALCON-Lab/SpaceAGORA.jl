---
id: gnc.target_energy_bracketing__edg_target_energy_from_reachable_bracket
label: _edg_target_energy_from_reachable_bracket
kind: function
source:
  file: src/gnc/guidance/target_energy_bracketing.jl
  symbol: _edg_target_energy_from_reachable_bracket
  lines:
  - 226
  - 226
inputs:
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: target_apoapsis_radius_m
  type: Float64
  units: n/a
  required: true
  description: Positional argument `target_apoapsis_radius_m`.
- id: energy_min
  type: Float64
  units: n/a
  required: true
  description: Positional argument `energy_min`.
- id: energy_max
  type: Float64
  units: n/a
  required: true
  description: Positional argument `energy_max`.
- id: periapsis_at_min
  type: Float64
  units: n/a
  required: true
  description: Positional argument `periapsis_at_min`.
- id: periapsis_at_max
  type: Float64
  units: n/a
  required: true
  description: Positional argument `periapsis_at_max`.
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
  type: Roots.find_zero
  units: n/a
  description: 'Return value of `_edg_target_energy_from_reachable_bracket`. Returns
    `exit_energy - desired_energy` or `Roots.find_zero(residual, (energy_min, energy_max),
    Roots.Brent(); rtol=1e-10)` or `energy_min - residual_min * (energy_max - energy_min)
    / (residual_max - residual` or `abs(residual_min) <= abs(residual_max) ? energy_min
    : energy_max`.'
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

# _edg_target_energy_from_reachable_bracket

## Purpose
`_edg_target_energy_from_reachable_bracket` solves for the exit specific energy at which the spacecraft would reach the configured target apoapsis radius, given the two reachable bracket endpoints (low-drag and maximum-depletion outcomes) and their associated periapsis radii. The result becomes `state.target_energy_jkg` and determines whether targeting is feasible on this passage.

## Theory & Math
Find $\varepsilon^*$ such that $\varepsilon^* = \varepsilon_{des}\big(r_a^{target}, r_p(\varepsilon^*)\big)$ where $r_p(\varepsilon)$ is the linearly interpolated periapsis and, from the vis-viva relation, $\varepsilon_{des} = -\mu / (r_a + r_p)$ with $\mu$ the planet gravitational parameter (m^3/s^2), $r_a$ the target apoapsis radius (m) and $r_p$ the periapsis radius (m). The residual $R(\varepsilon) = \varepsilon - \varepsilon_{des}(\varepsilon)$ is driven to zero by Brent's method when $R(\varepsilon_{min}) R(\varepsilon_{max}) \le 0$.

## Design & Implementation
Signature `(planet, target_apoapsis_radius_m::Float64, energy_min::Float64, energy_max::Float64, periapsis_at_min::Float64, periapsis_at_max::Float64)`. It defines the closure `residual(exit_energy)` that interpolates the periapsis with `_edg_interpolate_bracket_value`, asks the control module for `_edg_target_energy_from_apoapsis(planet, target_apoapsis_radius_m, periapsis)` (the vis-viva energy of the orbit with that apoapsis and periapsis), and returns `exit_energy - desired_energy`. It evaluates the residual at both endpoints; if both are finite and their product is `<= 0` it calls `Roots.find_zero(residual, (energy_min, energy_max), Roots.Brent(); rtol=1e-10)`. If both are finite but same-signed and differ by more than `eps(Float64)`, it returns the secant extrapolation `energy_min - residual_min * (energy_max - energy_min) / (residual_max - residual_min)`. Otherwise it returns whichever endpoint has the smaller absolute residual.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `target_apoapsis_radius_m` | Float64 | n/a | yes | Positional argument `target_apoapsis_radius_m`. |
| in | `energy_min` | Float64 | n/a | yes | Positional argument `energy_min`. |
| in | `energy_max` | Float64 | n/a | yes | Positional argument `energy_max`. |
| in | `periapsis_at_min` | Float64 | n/a | yes | Positional argument `periapsis_at_min`. |
| in | `periapsis_at_max` | Float64 | n/a | yes | Positional argument `periapsis_at_max`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Roots.find_zero | n/a | — | Return value of `_edg_target_energy_from_reachable_bracket`. Returns `exit_energy - desired_energy` or `Roots.find_zero(residual, (energy_min, energy_max), Roots.Brent(); rtol=1e-10)` or `energy_min - residual_min * (energy_max - energy_min) / (residual_max - residual` or `abs(residual_min) <= abs(residual_max) ? energy_min : energy_max`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.target_energy_bracketing__edg_run_target_energy_bracketing_bang|_edg_run_target_energy_bracketing!]] · `callees` → `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:309-309`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The secant extrapolation can return an energy far outside the bracket; the caller's reachability test then rejects it, but the intermediate value is still stored as `target_energy_jkg`. If `Roots.find_zero` fails to converge it throws and the exception is not caught here. A `NaN` `target_apoapsis_radius_m` (the config default) yields `NaN` residuals and the function silently returns `energy_min`. The 1e-10 relative tolerance on Brent applies to the energy root, not to the resulting apoapsis error.

## Provenance
Mapped from `src/gnc/guidance/target_energy_bracketing.jl` line 226.
