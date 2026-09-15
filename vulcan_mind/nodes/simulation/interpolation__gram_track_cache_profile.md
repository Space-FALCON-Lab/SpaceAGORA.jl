---
id: simulation.interpolation__gram_track_cache_profile
label: _gram_track_cache_profile
kind: function
source:
  file: src/simulation/callbacks/gram_track_cache/interpolation.jl
  symbol: _gram_track_cache_profile
  lines:
  - 6
  - 6
inputs:
- id: cfg
  type: GramTrackCacheConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: alt
  type: Float64
  units: n/a
  required: true
  description: Positional argument `alt`.
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
  description: Return value of `_gram_track_cache_profile`. Returns `cfg.entry_horizon_s,
    cfg.entry_alt_tol_m, cfg.entry_ang_tol_rad, cfg.entry_point` or `cfg.orbit_horizon_s,
    cfg.orbit_alt_tol_m, cfg.orbit_ang_tol_rad, cfg.orbit_point`.
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

# _gram_track_cache_profile

## Purpose
Selects the entry-phase or orbit-phase set of cache tolerances for the current altitude, so the track cache is tight and short-horizoned during atmospheric entry and loose and long-horizoned on orbit.

## Design & Implementation
Takes `cfg::GramTrackCacheConfig`, the parameter container `p`, and `alt::Float64` in metres. It converts the entry-interface altitude with `EI_m = p.args.environment_model.EI * 1e3`, since `EI` is stored in kilometres while `alt` is in metres. When `alt <= EI_m + cfg.transition_band_m` it returns the four entry-profile fields `(entry_horizon_s, entry_alt_tol_m, entry_ang_tol_rad, entry_points)`; otherwise it returns the matching `orbit_*` quadruple. The additive `transition_band_m` deliberately extends the entry profile above the interface so the switch happens before drag matters.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | GramTrackCacheConfig | n/a | yes | Positional argument `cfg`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `alt` | Float64 | n/a | yes | Positional argument `alt`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_gram_track_cache_profile`. Returns `cfg.entry_horizon_s, cfg.entry_alt_tol_m, cfg.entry_ang_tol_rad, cfg.entry_point` or `cfg.orbit_horizon_s, cfg.orbit_alt_tol_m, cfg.orbit_ang_tol_rad, cfg.orbit_point`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/callbacks/gram_track_cache/interpolation.jl`
- [[simulation.runtime__density_state_from_kinematics_bang|_density_state_from_kinematics!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:124-124`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The switch is a hard threshold with no hysteresis, so a trajectory loitering at exactly `EI_m + transition_band_m` flips profiles between successive calls and repeatedly invalidates the cache. The kilometre-to-metre conversion is hard-coded as `* 1e3` and silently produces wrong thresholds if `environment_model.EI` is ever changed to SI metres. Reaching through `p.args.environment_model` on every call couples this routine to the parameter object's nesting, and no check confirms that field exists.

## Provenance
Mapped from `src/simulation/callbacks/gram_track_cache/interpolation.jl` line 6.
