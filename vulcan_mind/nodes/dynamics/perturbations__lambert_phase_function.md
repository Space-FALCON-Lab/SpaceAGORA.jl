---
id: dynamics.perturbations__lambert_phase_function
label: _lambert_phase_function
kind: function
source:
  file: src/dynamics/coupled/perturbations.jl
  symbol: _lambert_phase_function
  lines:
  - 1175
  - 1175
inputs:
- id: alpha_rad
  type: Float64
  units: n/a
  required: true
  description: Positional argument `alpha_rad`.
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
  type: Float64
  units: n/a
  description: Return value of `_lambert_phase_function`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# _lambert_phase_function

## Purpose
The Lambertian sphere phase function giving the fraction of reflected sunlight seen at a phase angle, used to scale planetary albedo pressure.

## Theory & Math
$$
\Phi(\alpha) = \frac{\sin\alpha + (\pi - \alpha)\cos\alpha}{\pi}
$$

## Design & Implementation
Returns zero for non-finite input, clamps the angle into `[0, π]`, and evaluates `(sin α + (π - α) cos α) / π` floored at zero. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `alpha_rad` | Float64 | n/a | yes | Positional argument `alpha_rad`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_lambert_phase_function`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.perturbations_planetary_albedo_accel|planetary_albedo_accel]] · `callees` → `callers` · call · `src/dynamics/coupled/perturbations.jl:1211-1211`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/perturbations.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Treats the planet as a uniform Lambertian sphere; real albedo varies with surface and cloud cover.

## Provenance
Mapped from `src/dynamics/coupled/perturbations.jl` line 1175.
