---
id: simulation_a.vacuum_predicted_gram_build_vacuum_gram_cache__build_vacuum_gram_cache_bang
label: _build_vacuum_gram_cache!
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl
  symbol: _build_vacuum_gram_cache!
  lines:
  - 181
  - 236
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: seed_state
  type: Tuple{VacuumPredictedGRAMCache, SVector{3,Float64}, SVector{3,Float64}, Float64}
  units: m, m/s, s
  required: true
  description: Cache object to fill plus the current inertial position, inertial velocity
    and time from which the vacuum reference trajectory is propagated.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: cache_valid
  type: VacuumPredictedGRAMCache
  units: n/a
  description: Cache populated with knot times, vacuum altitudes and positions, log-density
    and temperature spline coefficients and wind vectors, with `valid` set true on
    success.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation-a
origin: agent
---
# _build_vacuum_gram_cache!

## Purpose
`_build_vacuum_gram_cache!` precomputes an atmosphere profile along a predicted vacuum trajectory so that repeated GRAM queries inside a solver step can be answered by spline evaluation instead of native calls. It is the expensive half of the vacuum-predicted GRAM cache; `_query_vacuum_gram_cache!` is the cheap half that decides whether the existing cache still applies.

## Theory & Math
The reference trajectory is propagated in vacuum under J2-perturbed gravity with a classical fourth-order Runge-Kutta step of size $h = T/(n-1)$:

$$y_{k+1} = y_k + \tfrac{h}{6}\left(k_1 + 2k_2 + 2k_3 + k_4\right)$$

Density is stored logarithmically, $y_i = \ln \max(\rho_i, 10^{-40})$, because atmospheric density is close to exponential in altitude and the logarithm is far better conditioned for polynomial interpolation. A natural cubic spline is then fitted over the uniform knots by solving the tridiagonal system

$$M_{i-1} + 4M_i + M_{i+1} = \frac{6}{h^2}\left(y_{i+1} - 2y_i + y_{i-1}\right), \qquad M_1 = M_n = 0$$

with the Thomas algorithm, giving $O(n)$ construction. Evaluation uses the standard piecewise form with $b_i = (y_{i+1}-y_i)/h - h(M_{i+1} + 2M_i)/6$, $c_i = M_i/2$ and $d_i = (M_{i+1}-M_i)/(6h)$. Recovering density as $\exp(\cdot)$ guarantees positivity for any interpolated value.

## Model & Assumptions
The prediction assumes the spacecraft follows an unperturbed J2 gravity trajectory over the cache horizon, which holds well at orbital altitudes where drag acceleration is negligible compared with gravity. The cache is only trusted while the true inertial position stays within a configured deviation radius of the interpolated vacuum position; beyond that the profile is rebuilt from the current state. Uniform knot spacing is required by both the spline solver and the constant-time index lookup used at query time.

## Design & Implementation
The routine marks the cache invalid up front and returns immediately for fewer than two points, so a partially built cache is never queryable. All eight cache arrays are resized once, then a single forward pass propagates the state with `_vacuum_rk4_step`, rotates each predicted position into the planet frame with `_planet_lpi_at`, converts to altitude, latitude and longitude via `rtolatlong`, and issues one `getDensity` call per knot. Spline coefficients for log-density and temperature are built by `_natural_cubic_spline_build!`; wind and vacuum altitude are interpolated linearly instead, since directional components do not benefit from higher-order fitting. Time bounds, spacing and the validity flag are written last.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `seed_state` | Tuple{VacuumPredictedGRAMCache, SVector{3,Float64}, SVector{3,Float64}, Float64} | m, m/s, s | yes | Cache object to fill plus the current inertial position, inertial velocity and time from which the vacuum reference trajectory is propagated. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `cache_valid` | VacuumPredictedGRAMCache | n/a | — | Cache populated with knot times, vacuum altitudes and positions, log-density and temperature spline coefficients and wind vectors, with `valid` set true on success. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.vacuum_predicted_gram__query_vacuum_gram_cache__query_vacuum_gram_cache_bang|_query_vacuum_gram_cache!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:261-261`

**Downstream**

- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:213-213`
- `callees` → [[environment.get_density|getDensity]] · `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:218-218`
- `callees` → [[simulation.planet_frame__planet_lpi_at|_planet_lpi_at]] · `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:211-211`
- `callees` → [[simulation.vacuum_predicted_gram__natural_cubic_spline_build_bang|_natural_cubic_spline_build!]] · `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:228-228`
- `callees` → [[simulation.vacuum_predicted_gram__vacuum_rk4_step|_vacuum_rk4_step]] · `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:224-224`
<!-- vulcan:connections:end -->

## Limitations
The build cost is exactly `n_pts` GRAM evaluations, so an aggressive point count can outweigh the savings the cache provides. Accuracy degrades once real drag or thrust bends the trajectory away from the vacuum reference, which the deviation check detects only in position, not in velocity. The fallback path at the end of `_query_vacuum_gram_cache!` issues a direct query if the rebuild fails to produce a valid cache.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:180-236`.
