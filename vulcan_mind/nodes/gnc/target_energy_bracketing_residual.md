---
id: gnc.target_energy_bracketing_residual
label: residual
kind: function
source:
  file: src/gnc/guidance/target_energy_bracketing.jl
  symbol: residual
  lines:
  - 234
  - 234
inputs:
- id: exit_energy
  type: Any
  units: n/a
  required: true
  description: Positional argument `exit_energy`.
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
  description: Return value of `residual`. Returns `exit_energy - desired_energy`.
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

# residual

## Purpose
`residual` is the closure inside `_edg_target_energy_from_reachable_bracket` whose zero is the exit specific energy consistent with the target apoapsis. It is evaluated at the bracket endpoints to decide the solution strategy and, when the signs differ, handed to `Roots.find_zero` with Brent's method.

## Design & Implementation
Signature `residual(exit_energy)`, capturing `energy_min`, `energy_max`, `periapsis_at_min`, `periapsis_at_max`, `planet` and `target_apoapsis_radius_m` from the enclosing scope. It computes `periapsis = _edg_interpolate_bracket_value(exit_energy, energy_min, energy_max, periapsis_at_min, periapsis_at_max)` (m), then `desired_energy = _control_module()._edg_target_energy_from_apoapsis(planet, target_apoapsis_radius_m, periapsis)` (J/kg), and returns `exit_energy - desired_energy`. A positive value means the candidate exit energy is higher than the orbit through the target apoapsis and interpolated periapsis would have. The closure allocates nothing beyond the module lookup.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `exit_energy` | Any | n/a | yes | Positional argument `exit_energy`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `residual`. Returns `exit_energy - desired_energy`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.heat_load_control__edg_solve_heat_load_switches|_edg_solve_heat_load_switches]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:648-648`

**Downstream**

- `callees` → [[gnc.guidance_hooks__control_module|_control_module]] · `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:236-236`
- `callees` → [[gnc.target_energy_bracketing__edg_interpolate_bracket_value|_edg_interpolate_bracket_value]] · `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:235-235`
<!-- vulcan:connections:end -->

## Limitations
Each evaluation performs a `_control_module()` lookup, a dynamic call that is not inlined and adds overhead inside Brent iterations. The residual is only as smooth as the linear periapsis interpolation, so it is piecewise-linear in `exit_energy` and Brent may terminate on a kink. `NaN` inputs (unset target apoapsis) yield `NaN` residuals, which the enclosing function handles by falling back to an endpoint rather than raising.

## Provenance
Mapped from `src/gnc/guidance/target_energy_bracketing.jl` line 234.
