---
id: gnc.target_energy_bracketing__edg_interpolate_bracket_value
label: _edg_interpolate_bracket_value
kind: function
source:
  file: src/gnc/guidance/target_energy_bracketing.jl
  symbol: _edg_interpolate_bracket_value
  lines:
  - 219
  - 219
inputs:
- id: exit_energy
  type: Float64
  units: n/a
  required: true
  description: Positional argument `exit_energy`.
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
- id: value_at_min
  type: Float64
  units: n/a
  required: true
  description: Positional argument `value_at_min`.
- id: value_at_max
  type: Float64
  units: n/a
  required: true
  description: Positional argument `value_at_max`.
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
  description: Return value of `_edg_interpolate_bracket_value`. Returns `value_at_min
    + fraction * (value_at_max - value_at_min)`.
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

# _edg_interpolate_bracket_value

## Purpose
`_edg_interpolate_bracket_value` linearly interpolates a quantity (in practice the post-passage periapsis radius) between the two bracket endpoints as a function of exit specific energy. `_edg_target_energy_from_reachable_bracket` uses it to estimate the periapsis that would accompany a candidate exit energy.

## Theory & Math
For exit energy $\varepsilon$ between bracket energies $\varepsilon_{min}$ and $\varepsilon_{max}$ with endpoint values $y_{min}$, $y_{max}$:

$$y(\varepsilon) = y_{min} + \frac{\varepsilon - \varepsilon_{min}}{\varepsilon_{max} - \varepsilon_{min}}\,(y_{max} - y_{min})$$

falling back to $\tfrac{1}{2}(y_{min}+y_{max})$ when $|\varepsilon_{max}-\varepsilon_{min}| < \epsilon_{mach}$.

## Design & Implementation
Signature `(exit_energy::Float64, energy_min::Float64, energy_max::Float64, value_at_min::Float64, value_at_max::Float64)`. It computes `width = energy_max - energy_min`; when `abs(width) < eps(Float64)` it returns the midpoint `0.5 * (value_at_min + value_at_max)` to avoid division by zero. Otherwise `fraction = (exit_energy - energy_min) / width` and the result is `value_at_min + fraction * (value_at_max - value_at_min)`. The function is pure and does not clamp `fraction`, so it extrapolates linearly outside the bracket.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `exit_energy` | Float64 | n/a | yes | Positional argument `exit_energy`. |
| in | `energy_min` | Float64 | n/a | yes | Positional argument `energy_min`. |
| in | `energy_max` | Float64 | n/a | yes | Positional argument `energy_max`. |
| in | `value_at_min` | Float64 | n/a | yes | Positional argument `value_at_min`. |
| in | `value_at_max` | Float64 | n/a | yes | Positional argument `value_at_max`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_edg_interpolate_bracket_value`. Returns `value_at_min + fraction * (value_at_max - value_at_min)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.target_energy_bracketing_residual|residual]] · `callees` → `callers` · call · `src/gnc/guidance/target_energy_bracketing.jl:235-235`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/target_energy_bracketing.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The degenerate-width threshold is the absolute machine epsilon (about 2.2e-16 J/kg), far below realistic energy resolution, so nearly-equal brackets still divide by a tiny width and amplify noise in the endpoint values. Linear dependence of periapsis on exit energy is an approximation; the real relation through the drag passage is nonlinear. No finiteness checks exist, so `NaN` or `Inf` inputs propagate silently.

## Provenance
Mapped from `src/gnc/guidance/target_energy_bracketing.jl` line 219.
