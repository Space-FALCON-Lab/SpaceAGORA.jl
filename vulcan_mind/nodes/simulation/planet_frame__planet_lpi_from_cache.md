---
id: simulation.planet_frame__planet_lpi_from_cache
label: _planet_lpi_from_cache
kind: function
source:
  file: src/simulation/callbacks/density_callbacks/planet_frame.jl
  symbol: _planet_lpi_from_cache
  lines:
  - 8
  - 8
inputs:
- id: cache
  type: PlanetFrameEphemerisCache
  units: n/a
  required: true
  description: Positional argument `cache`.
- id: et
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et`.
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
  type: Union{Nothing,
  units: n/a
  description: Return value of `_planet_lpi_from_cache`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _planet_lpi_from_cache

## Purpose
`_planet_lpi_from_cache(cache::PlanetFrameEphemerisCache, et)` serves the planet-fixed-from-inertial rotation from a precomputed table of attitude quaternions, avoiding a backend call on the hot per-step path. It returns `nothing` whenever the request cannot be satisfied from the table, which is the caller's signal to fall through to the backend.

## Theory & Math
Given bracketing samples $(t_0, q_0)$ and $(t_1, q_1)$ with $t_0 \le t \le t_1$, let $\alpha = (t - t_0)/(t_1 - t_0) \in [0,1]$ and take $q_1 \leftarrow -q_1$ whenever $\langle q_0, q_1 \rangle < 0$. The interpolant is normalised LERP (NLERP),

$$q(t) = \frac{(1-\alpha) q_0 + \alpha q_1}{\lVert (1-\alpha) q_0 + \alpha q_1 \rVert},$$

and $L_{PI} = R(q(t))$ is the corresponding direction-cosine matrix. NLERP follows the same great-circle path on $S^3$ as SLERP but traverses it non-uniformly; the angular-rate error relative to SLERP grows with the half-angle $\theta$ between $q_0$ and $q_1$, vanishing as $\theta \to 0$.

## Design & Implementation
Four guards precede any work: the table must hold at least two samples (`n_samples >= 2`), `et` must lie within `[ets[1], ets[n_samples]]`, and the bracketing index from `searchsortedlast(ets, et)` must be in range — `idx <= 0` returns `nothing`, `idx >= n_samples` returns `rot(cache.quaternions[n_samples])` exactly at the upper endpoint, and a degenerate interval (`et1 <= et0`) returns the left sample without interpolating. Otherwise it forms the interpolation parameter `α = (et - et0) / (et1 - et0)`, takes `q0` and `q1`, and — critically — negates `q1` when `dot(q0, q1) < 0.0` so both quaternions lie on the same hemisphere of the double cover; without this the blend would sweep the long way round and produce a visible attitude discontinuity. The blended quaternion is the normalised linear interpolation `normalize((1 - α) * q0 + α * q1)`, converted to a matrix by `rot`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cache` | PlanetFrameEphemerisCache | n/a | yes | Positional argument `cache`. |
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_planet_lpi_from_cache`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl`
- [[simulation.planet_frame__planet_lpi_at|_planet_lpi_at]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:47-47`

**Downstream**

- `callees` → [[core.quaternion_utils_rot|rot]] · `callers` · call · `src/simulation/callbacks/density_callbacks/planet_frame.jl:20-20`
<!-- vulcan:connections:end -->

## Limitations
NLERP is not constant-angular-velocity, so the interpolated rotation lags and leads the true one within each interval; the error is only negligible while consecutive samples are closely spaced. The hemisphere fix is applied pairwise, which keeps each interval continuous but does not guarantee a globally consistent sign convention across the table. Requests outside `[ets[1], ets[end]]` are rejected rather than extrapolated, silently pushing the caller onto the expensive backend for any time outside coverage. `searchsortedlast` assumes `ets` is sorted ascending — an unsorted table produces a wrong bracket with no error. No check is made that `length(quaternions) == length(ets)`, so a ragged cache yields a `BoundsError` deep in the interpolation.

## Provenance
Mapped from `src/simulation/callbacks/density_callbacks/planet_frame.jl` line 8.
