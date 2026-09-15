---
id: vehicle.thermal_models_heatrate_radiative
label: heatrate_radiative
kind: function
source:
  file: src/vehicle/thermal/thermal_models.jl
  symbol: heatrate_radiative
  lines:
  - 37
  - 37
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
  description: Return value of `heatrate_radiative`. Returns `q_rad`.
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

# heatrate_radiative

## Purpose

`heatrate_radiative(S, T, m, ρ, v, α)` computes the stagnation-point radiative heat rate on a blunt body using a Tauber-Sutton style correlation, in which a tabulated velocity function is combined with power laws in nose radius and free-stream density. It complements the convective correlation in the same file for high-speed entry heating estimates.

## Design & Implementation

Two hard-coded 20-element vectors are built each call: `vf`, velocities from 16000 m/s down to 0, and `fv`, the matching velocity function values from 2040 down to 0. `linear_interpolation(sort(vf), sort(fv))` builds the interpolant and `f = fn(v)` samples it. The density exponent is fixed at `b = 1.22` and the coefficient at `C = 4.736e4`, while the nose-radius exponent is itself velocity- and density-dependent: `a = 1.072e6 * v^(-1.88) * ρ^(-0.325)`. The returned value is `C * rn^a * ρ^b * f` with `rn = m.body.nose_radius` taken without the offset that `heatrate_convective` applies. Arguments `S`, `T` and `α` are unused.

## Theory & Math

The implemented correlation is

$$q_{rad} = C\, r_n^{\,a}\, \rho^{\,b}\, f(v), \qquad a = 1.072\times10^{6}\, v^{-1.88}\, \rho^{-0.325}$$

with $C = 4.736\times10^{4}$, $b = 1.22$, $r_n$ the nose radius (m), $\rho$ the free-stream density (kg/m$^3$), $v$ the relative velocity (m/s) and $f(v)$ the piecewise-linear velocity function interpolated from the tabulated pairs. The resulting $q_{rad}$ carries the units implied by the fit constants, W/cm$^2$.

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
| out | `result` | Any | n/a | — | Return value of `heatrate_radiative`. Returns `q_rad`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/thermal/thermal_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

Sorting `vf` and `fv` independently is only safe because both tables happen to be monotonically decreasing together; a table edit that breaks that co-ordering would silently mispair velocity and function value. Evaluating outside 0-16000 m/s relies on `Interpolations` extrapolation behaviour rather than an explicit guard. The tables and exponents are fitted for air, so use with a non-Earth `planet` is unvalidated. The reallocation of both vectors and reconstruction of the interpolant on every call makes this unsuitable for a tight integration loop.

## Provenance
Mapped from `src/vehicle/thermal/thermal_models.jl` line 37.
