---
id: environment.density_models__exponential_density
label: _exponential_density
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: _exponential_density
  lines:
  - 474
  - 474
inputs:
- id: rho_ref
  type: Float64
  units: n/a
  required: true
  description: Positional argument `ρ_ref`.
- id: h_ref
  type: Float64
  units: n/a
  required: true
  description: Positional argument `h_ref`.
- id: H
  type: Float64
  units: n/a
  required: true
  description: Positional argument `H`.
- id: h
  type: Float64
  units: n/a
  required: true
  description: Positional argument `h`.
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
  description: Return value of `_exponential_density`.
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

# _exponential_density

## Purpose
The scalar exponential density law shared by the single-layer and multi-layer analytic models and by their batch evaluation methods.

## Theory & Math
$$
\rho(h) = \rho_{\text{ref}}\, e^{(h_{\text{ref}} - h)/H}
$$

## Design & Implementation
Returns `ρ_ref * exp((h_ref - h) / H)` for reference density, reference altitude, scale height and query altitude all in SI units. Declared `@inline` so the batch loops in `getDensityBatch!` fuse it into a single pass with no function-call overhead per satellite.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `rho_ref` | Float64 | n/a | yes | Positional argument `ρ_ref`. |
| in | `h_ref` | Float64 | n/a | yes | Positional argument `h_ref`. |
| in | `H` | Float64 | n/a | yes | Positional argument `H`. |
| in | `h` | Float64 | n/a | yes | Positional argument `h`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_exponential_density`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.density_models_getdensitybatch_bang|getDensityBatch!]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:939-939`
- [[environment.density_models_timetabulatedatmospheremodel|TimeTabulatedAtmosphereModel]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:813-813`
- [[environment.get_density|getDensity]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:813-813`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
There is no guard against a zero or negative scale height, which would yield infinite or growing density; the model constructors validate `H > 0` so this is safe only through them.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 474.
