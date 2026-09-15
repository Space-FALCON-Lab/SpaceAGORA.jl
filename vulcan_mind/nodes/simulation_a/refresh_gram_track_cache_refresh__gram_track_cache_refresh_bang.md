---
id: simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang
label: _gram_track_cache_refresh!
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/refresh.jl
  symbol: _gram_track_cache_refresh!
  lines:
  - 57
  - 293
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: refresh_request
  type: Tuple{GramTrackCache, Any, SVector{3,Float64}, SVector{3,Float64}, Float64}
  units: m, m/s, s
  required: true
  description: Cache to rebuild, the active density model, and the current inertial
    position, velocity and time, with the regime horizon, point count and tolerances
    that shape the new segment.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: density_state
  type: Tuple{Float64, Float64, SVector{3,Float64}}
  units: kg/m^3, K, m/s
  description: Density, temperature and wind at the requested state, produced from
    the rebuilt cache or from a direct GRAM query when the rebuild fails.
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
# _gram_track_cache_refresh!

## Purpose
`_gram_track_cache_refresh!` rebuilds the GRAM track cache after a miss and returns the atmosphere state at the requested point. It is the orchestrator of the track-cache subsystem: it picks a segment endpoint, decides how many samples that segment needs, populates the cache from a propagated trajectory, and guarantees an answer even when every one of those steps fails.

## Theory & Math
Segment length is chosen from two competing constraints and the larger wins. The time constraint takes the per-sample spacing implied by the configured horizon, $\Delta t = T_{\mathrm{base}}/(n_{\mathrm{base}}-1)$, and requires $n_t = \lceil \Delta t_{\mathrm{seg}}/\Delta t \rceil + 1$. The arc-length constraint estimates the ground-track distance to the predicted endpoint and requires $n_\ell = \lceil L/s \rceil + 1$, where the target spacing $s$ is derived from the altitude and angular tolerances at the local radius. The final count is $n = \min\!\left(\max(n_{\mathrm{base}}, n_t, n_\ell),\, n_{\max}\right)$, so a fast-moving entry arc is sampled densely enough to stay inside the interpolation tolerance while an orbital arc is not oversampled.

## Model & Assumptions
Endpoint selection is regime-dependent. Inside the atmospheric band, defined as altitude below the entry-interface altitude plus the transition band, the routine first tries `_gram_periapsis_target` so the segment terminates at the drag-pass low point; for entry-like or open trajectories it falls back to Allen-Eggers endpoint targeting when that mode is enabled, using the spacecraft mass and reference area to form a ballistic coefficient, and finally to solver-endpoint propagation. For orbit missions it targets a full orbital period through `_gram_orbit_period_target`. For time missions it propagates to the remaining solver span when that is finite, and otherwise to one orbital period. Keplerian targeting with optional J2 secular rates is the base predictor, with a linear extrapolation fallback when the Kepler solve fails.

## Design & Implementation
The whole body runs inside a `try` block so any failure — a targeting solve that does not converge, a GRAM call that throws — is caught, the cache is marked invalid, the failure counter is incremented, a one-shot `@warn` is emitted through the `_gram_track_cache_warning_emitted` latch, and a direct `getDensity` call supplies the answer. Statistics collection is compiled out cheaply when profiling is disabled: `stats_enabled` is checked before every `time_ns` call and before every `_gram_runtime_stats_update!`. Refresh call count, total and maximum sample counts, and elapsed seconds all feed the `GramRuntimeStats` record.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `refresh_request` | Tuple{GramTrackCache, Any, SVector{3,Float64}, SVector{3,Float64}, Float64} | m, m/s, s | yes | Cache to rebuild, the active density model, and the current inertial position, velocity and time, with the regime horizon, point count and tolerances that shape the new segment. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `density_state` | Tuple{Float64, Float64, SVector{3,Float64}} | kg/m^3, K, m/s | — | Density, temperature and wind at the requested state, produced from the rebuilt cache or from a direct GRAM query when the rebuild fails. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.runtime__density_state_from_kinematics_bang|_density_state_from_kinematics!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:174-174`

**Downstream**

- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:236-236`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:237-237`
- `callees` → [[environment.density_models_getdensitybatch_bang|getDensityBatch!]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:255-255`
- `callees` → [[environment.get_density|getDensity]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:291-291`
- `callees` → [[simulation.config__gram_entry_target_cd|_gram_entry_target_cd]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:115-115`
- `callees` → [[simulation.config__gram_entry_target_mode|_gram_entry_target_mode]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:99-99`
- `callees` → [[simulation.config__gram_track_trajectory_supported|_gram_track_trajectory_supported]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:195-195`
- `callees` → [[simulation.refresh__gram_kepler_or_linear_target|_gram_kepler_or_linear_target]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:89-89`
- `callees` → [[simulation.refresh__gram_track_cache_fill_from_trajectory_bang|_gram_track_cache_fill_from_trajectory!]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:226-226`
- `callees` → [[simulation.registry__gram_runtime_stats_enabled|_gram_runtime_stats_enabled]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:77-77`
- `callees` → [[simulation.registry__gram_runtime_stats_update_bang|_gram_runtime_stats_update!]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:177-177`
- `callees` → [[simulation.targeting__angle_delta_rad|_angle_delta_rad]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:200-200`
- `callees` → [[simulation.targeting__gram_entry_mass_kg|_gram_entry_mass_kg]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:106-106`
- `callees` → [[simulation.targeting__gram_entry_reference_area_m2|_gram_entry_reference_area_m2]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:107-107`
- `callees` → [[simulation.targeting__gram_expected_track_length_m|_gram_expected_track_length_m]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:162-162`
- `callees` → [[simulation.targeting__gram_orbit_period_target|_gram_orbit_period_target]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:136-136`
- `callees` → [[simulation.targeting__gram_periapsis_target|_gram_periapsis_target]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:93-93`
- `callees` → [[simulation.targeting__gram_track_cache_max_npos|_gram_track_cache_max_npos]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:174-174`
- `callees` → [[simulation.targeting__gram_track_cache_periapsis_split_enabled|_gram_track_cache_periapsis_split_enabled]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:91-91`
- `callees` → [[simulation.targeting__gram_track_cache_target_spacing_m|_gram_track_cache_target_spacing_m]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:171-171`
- `callees` → [[simulation_a.model_selection_gram_isolated_pool_batch_eval__gram_isolated_pool_batch_eval_bang|_gram_isolated_pool_batch_eval!]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:243-243`
- `callees` → [[simulation_a.targeting_gram_entry_target_allen_eggers|_gram_entry_target_allen_eggers]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:108-108`
<!-- vulcan:connections:end -->

## Limitations
Refresh is the dominant cost of the track cache and the reason the feature ships disabled by default; a recorded entry benchmark shows the cached path an order of magnitude slower than direct point-to-point sampling. The warning latch fires once per process, so repeated distinct failures are not individually reported. Endpoint predictors are unperturbed apart from optional J2, so a thrusting vehicle can invalidate the segment faster than the tolerance test expects.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/refresh.jl:56-293`.
