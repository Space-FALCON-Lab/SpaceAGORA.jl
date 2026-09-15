---
id: simulation_a.interpolation_gram_track_cache_ready
label: _gram_track_cache_ready
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/interpolation.jl
  symbol: _gram_track_cache_ready
  lines:
  - 92
  - 114
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: query
  type: Tuple{GramTrackCache, Float64, Float64, Float64, Float64}
  units: s, m, rad, rad
  required: true
  description: Cache to test plus the query time, altitude, latitude and longitude
    of the current spacecraft state, with the altitude and angular tolerances that
    define a hit.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: segment_hit
  type: Union{Nothing, Tuple{Int, Float64}}
  units: n/a
  description: Segment index and normalised in-segment position when the cached track
    still represents the queried state, or `nothing` to force a direct atmosphere
    query.
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
# _gram_track_cache_ready

## Purpose
`_gram_track_cache_ready` is the gate every cached GRAM query passes through. It answers one question: does the precomputed ground track still describe where the spacecraft actually is, closely enough that interpolating the cached atmosphere is defensible instead of calling GRAM again?

## Theory & Math
Given a bracketing segment $[t_i, t_{i+1}]$ the normalised abscissa is $x = (t - t_i)/(t_{i+1} - t_i)$, and each tracked quantity is interpolated linearly, $\hat{y} = y_i + x(y_{i+1} - y_i)$. Angles use short-arc interpolation: the difference $d = b - a$ is wrapped into $(-\pi, \pi]$ before the blend, and the result is renormalised, avoiding a spurious sweep across the antimeridian. The hit test is the conjunction

$$\lvert h - \hat{h} \rvert \le \epsilon_h \;\wedge\; \delta(\varphi, \hat{\varphi}) \le \epsilon_\theta \;\wedge\; \delta(\lambda, \hat{\lambda}) \le \epsilon_\theta$$

where $\delta(a,b) = \lvert \operatorname{atan2}(\sin(b-a), \cos(b-a)) \rvert$ is the absolute short-arc angular separation, computed by the cheap branch when the raw difference already lies in $(-\pi,\pi]$.

## Model & Assumptions
The test treats the cache as valid only when the state matches in all three of altitude, latitude and longitude within the regime-specific tolerances supplied by `GramTrackCacheConfig`. Angular comparisons use short-arc separation so a track crossing the antimeridian or the prime meridian is handled without a false miss. Linear interpolation between knots is assumed adequate, which is why the refresh logic sizes segments so that consecutive points stay inside the tolerance envelope.

## Design & Implementation
The function first calls `_gram_track_cache_segment`, which locates the bracketing interval. That lookup is optimised for the integrator's access pattern: because time advances monotonically, it starts from the stored `index_hint`, tries the current segment, then the immediately following one, and only falls back to `searchsortedlast` when time moved backwards after a rejected step. The successful index is written back into `cache.index_hint`, making the common case constant-time. When `ignore_time_window` is set the query time is clamped into the cached span rather than rejected. The altitude blend uses `_lerp` while latitude and longitude use `_lerp_angle_rad`, and the companion `_gram_track_cache_eval` reuses the same index and abscissa to produce density, temperature and a wind vector assembled with `muladd`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `query` | Tuple{GramTrackCache, Float64, Float64, Float64, Float64} | s, m, rad, rad | yes | Cache to test plus the query time, altitude, latitude and longitude of the current spacecraft state, with the altitude and angular tolerances that define a hit. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `segment_hit` | Union{Nothing, Tuple{Int, Float64}} | n/a | — | Segment index and normalised in-segment position when the cached track still represents the queried state, or `nothing` to force a direct atmosphere query. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.runtime__density_state_from_kinematics_bang|_density_state_from_kinematics!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:167-167`

**Downstream**

- `callees` → [[simulation.config__gram_track_cache_ignore_time_window|_gram_track_cache_ignore_time_window]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/interpolation.jl:100-100`
- `callees` → [[simulation.interpolation__angdiff_rad|_angdiff_rad]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/interpolation.jl:109-109`
- `callees` → [[simulation.interpolation__gram_track_cache_segment|_gram_track_cache_segment]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/interpolation.jl:102-102`
- `callees` → [[simulation.interpolation__lerp|_lerp]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/interpolation.jl:105-105`
- `callees` → [[simulation.interpolation__lerp_angle_rad|_lerp_angle_rad]] · `callers` · call · `src/simulation/callbacks/gram_track_cache/interpolation.jl:106-106`
<!-- vulcan:connections:end -->

## Limitations
A hit certifies geometric proximity only; a genuinely different atmospheric state at the same coordinates, such as a time-varying dust or solar condition, is not detected. Caches shorter than two points and invalidated caches always miss. Mutating `index_hint` makes the cache object non-thread-safe, so each spacecraft owns its own instance in `p.shared_buffers.gram_density_cache`.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/interpolation.jl:92-114`.
