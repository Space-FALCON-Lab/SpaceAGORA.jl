---
id: vehicle.thermal_models_heatrate_convective
label: heatrate_convective
kind: function
source:
  file: src/vehicle/thermal/thermal_models.jl
  symbol: heatrate_convective
  lines:
  - 12
  - 12
inputs:
- id: S
  type: Any
  units: n/a
  required: true
  description: Positional argument `S`.
- id: T
  type: Any
  units: n/a
  required: true
  description: Positional argument `T`.
- id: m
  type: Any
  units: n/a
  required: true
  description: Positional argument `m`.
- id: rho
  type: Any
  units: n/a
  required: true
  description: Positional argument `ρ`.
- id: v
  type: Any
  units: n/a
  required: true
  description: Positional argument `v`.
- id: alpha
  type: Any
  units: n/a
  required: true
  description: Positional argument `α`.
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
  description: Return value of `heatrate_convective`. Returns `q_conv`.
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

# heatrate_convective

## Purpose

`heatrate_convective(S, T, m, ρ, v, α)` returns the stagnation-point convective heat rate on a blunt body using a Sutton-Graves style correlation. It is the continuum-regime counterpart to `MaxwellianHeat`, driven by free-stream density and relative velocity rather than by a speed ratio, and it is the only term that `heatrate_convective_radiative` currently propagates.

## Design & Implementation

The function reads the effective nose radius as `m.body.nose_radius + 0.25` (metres; the fixed 0.25 m offset is hard-coded, not derived) and the Sutton-Graves coefficient `k` from `m.planet.k`. It then evaluates `q_conv = k * sqrt(ρ / rn) * v^3 * 1e-4` and returns that scalar. Of the six arguments only `m`, `ρ` and `v` participate: surface area `S`, temperature `T` and angle of attack `α` are accepted for signature uniformity with `heatrate_radiative` and `heatrate_convective_radiative` but are unused. The `1e-4` factor converts W/m² to W/cm².

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `S` | Any | n/a | yes | Positional argument `S`. |
| in | `T` | Any | n/a | yes | Positional argument `T`. |
| in | `m` | Any | n/a | yes | Positional argument `m`. |
| in | `rho` | Any | n/a | yes | Positional argument `ρ`. |
| in | `v` | Any | n/a | yes | Positional argument `v`. |
| in | `alpha` | Any | n/a | yes | Positional argument `α`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `heatrate_convective`. Returns `q_conv`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/thermal/thermal_models.jl`
- [[vehicle.thermal_models_heatrate_convective_radiative|heatrate_convective_radiative]] · `callees` → `callers` · call · `src/vehicle/thermal/thermal_models.jl:89-89`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

The correlation is a stagnation-point fit and gives no distribution over the body, so surface area `S` cannot influence the result. The 0.25 m nose-radius offset biases every vehicle, and for a sharp body it dominates `nose_radius` entirely. Wall temperature is ignored, so the result is a cold-wall flux with no hot-wall enthalpy correction. The units returned are W/cm² despite the in-body comment block claiming W/m². The correlation assumes continuum, laminar, equilibrium air-like flow and is not valid in the rarefied regime.

## Provenance
Mapped from `src/vehicle/thermal/thermal_models.jl` line 12.
