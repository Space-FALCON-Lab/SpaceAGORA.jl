---
id: environment.density_models__tab_flight_interp
label: _tab_flight_interp
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _tab_flight_interp
  lines:
  - 690
  - 690
inputs:
- id: alts
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `alts`.
- id: logs
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `logs`.
- id: sigs
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `sigs`.
- id: h
  type: Float64
  units: n/a
  required: true
  description: Positional argument `h`.
- id: sigma_scale
  type: Float64
  units: n/a
  required: true
  description: Positional argument `sigma_scale`.
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
  type: Tuple{Float64,
  units: n/a
  description: Return value of `_tab_flight_interp`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- environment
origin: agent
---

# _tab_flight_interp

## Purpose
Interpolates one flight profile at an altitude, returning shifted log-density and the local scale height, with exponential extrapolation beyond the profile's coverage.

## Theory & Math
$$
H = \frac{h_{j+1} - h_j}{\ln\rho_j - \ln\rho_{j+1}},\qquad \ln\rho(h) = \ln\rho_j + t\,(\ln\rho_{j+1} - \ln\rho_j) + s\,\sigma(h)
$$

## Design & Implementation
Returns `NaN` pairs for a `NaN` altitude. Below the first or above the last sample it extrapolates with the end segment's scale height clamped to 2 km through 12 km, defaulting to 8 km for a single-point profile. Inside, it finds the bracketing samples with `searchsortedlast`, interpolates log-density and sigma linearly, and computes the scale height from the segment slope guarded by 1e-9. The sigma shift is added last.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `alts` | Vector{Float64} | n/a | yes | Positional argument `alts`. |
| in | `logs` | Vector{Float64} | n/a | yes | Positional argument `logs`. |
| in | `sigs` | Vector{Float64} | n/a | yes | Positional argument `sigs`. |
| in | `h` | Float64 | n/a | yes | Positional argument `h`. |
| in | `sigma_scale` | Float64 | n/a | yes | Positional argument `sigma_scale`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{Float64, | n/a | — | Return value of `_tab_flight_interp`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.get_density|getDensity]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:727-727`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The scale-height guard against a non-decreasing density produces a huge `H` that the clamp then hides, so an inverted profile segment reports the 12 km maximum instead of an error.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 690.
