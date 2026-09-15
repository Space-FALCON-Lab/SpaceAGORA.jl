---
id: flow.environment
label: Environment sampling
kind: group
inputs:
- id: kernels
  type: SPICE kernels
  units: n/a
  description: Furnished kernel set.
- id: gram_assets
  type: GRAM datasets
  units: n/a
  description: Atmosphere data.
- id: coefficient_file
  type: CSV
  units: n/a
  description: Harmonics coefficients.
outputs:
- id: atmosphere
  type: (rho, T, wind)
  units: n/a
  description: Density, temperature and wind at a point and time.
- id: frames_and_bodies
  type: l_pi, body positions
  units: n/a
  description: Planet orientation and third-body positions.
tags:
- master-flow
charts:
- master
origin: agent
opens: environment
---

# Environment sampling

## Purpose
Answers the questions the dynamics ask about the world: what is the atmosphere here, where is the Sun and every other body, and how is the planet oriented — through interchangeable models from analytic exponentials to native GRAM, and from a simple rotating frame to SPICE.

## Design & Implementation
`density_models.jl` provides `NoAtmosphere`, exponential, piecewise, polynomial-fit, NRLMSISE-00, tabulated flight and time series, and the GRAM wrappers, all behind `getDensity`/`getDensityBatch!`; `ephemerides/` provides planet constants, shapes, the simple analytic ephemeris and the SPICE-backed one; `gravity/` provides the field constants. Runtime caching — track caches, vacuum-predicted splines, ephemeris tables and the one-entry memo — sits in the simulation layer above.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `kernels` | SPICE kernels | n/a | — | Furnished kernel set. |
| in | `gram_assets` | GRAM datasets | n/a | — | Atmosphere data. |
| in | `coefficient_file` | CSV | n/a | — | Harmonics coefficients. |
| out | `atmosphere` | (rho, T, wind) | n/a | — | Density, temperature and wind at a point and time. |
| out | `frames_and_bodies` | l_pi, body positions | n/a | — | Planet orientation and third-body positions. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[input.gram_data|GRAM atmosphere data]] · `gram_assets` → `gram_assets` · dataflow · `src/environment/atmosphere/density_models.jl`
- [[input.gravity_coefficients|Gravity harmonics coefficients]] · `coefficient_file` → `coefficient_file` · dataflow · `src/dynamics/coupled/perturbations.jl`
- [[input.spice_kernels|SPICE kernels]] · `kernels` → `kernels` · dataflow · `src/environment/ephemerides/planets.jl`

**Downstream**

- `atmosphere` → [[flow.callbacks|Integration callbacks]] · `callback_hooks` · dataflow · `src/simulation/callbacks/density_callbacks/runtime.jl`
- `frames_and_bodies` → [[flow.forces|Force & torque models]] · `force_requests` · dataflow · `src/simulation/engine/effector_sampling.jl`
<!-- vulcan:connections:end -->

## Limitations
Every native GRAM or CSPICE call serialises on a process-wide lock, so fidelity is bought with thread contention; the analytic models have no latitude, longitude or time dependence.
