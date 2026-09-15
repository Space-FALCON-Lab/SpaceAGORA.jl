---
id: core.runtime_types_gramtrackcacheconfig
label: GramTrackCacheConfig
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: GramTrackCacheConfig
  lines:
  - 590
  - 590
inputs:
- id: mode
  type: Symbol
  units: n/a
  required: true
  description: Field `mode`.
- id: entry_horizon_s
  type: Float64
  units: n/a
  required: true
  description: Field `entry_horizon_s`.
- id: entry_alt_tol_m
  type: Float64
  units: n/a
  required: true
  description: Field `entry_alt_tol_m`.
- id: entry_ang_tol_rad
  type: Float64
  units: n/a
  required: true
  description: Field `entry_ang_tol_rad`.
- id: entry_points
  type: Int
  units: n/a
  required: true
  description: Field `entry_points`.
- id: orbit_horizon_s
  type: Float64
  units: n/a
  required: true
  description: Field `orbit_horizon_s`.
- id: orbit_alt_tol_m
  type: Float64
  units: n/a
  required: true
  description: Field `orbit_alt_tol_m`.
- id: orbit_ang_tol_rad
  type: Float64
  units: n/a
  required: true
  description: Field `orbit_ang_tol_rad`.
- id: orbit_points
  type: Int
  units: n/a
  required: true
  description: Field `orbit_points`.
- id: transition_band_m
  type: Float64
  units: n/a
  required: true
  description: Field `transition_band_m`.
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
  type: GramTrackCacheConfig
  units: n/a
  description: Constructed `GramTrackCacheConfig`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# GramTrackCacheConfig

## Purpose
Immutable snapshot of the tolerances governing when the GRAM along-track density cache may be used, for both the atmospheric-entry regime and the orbital regime.

## Design & Implementation
`struct GramTrackCacheConfig` with `mode::Symbol`, then two parallel sets of knobs: `entry_horizon_s`, `entry_alt_tol_m`, `entry_ang_tol_rad`, `entry_points::Int` and `orbit_horizon_s`, `orbit_alt_tol_m`, `orbit_ang_tol_rad`, `orbit_points::Int`, plus `transition_band_m` (m), the altitude band in which the two regimes blend. Built by `SimulationCallbacks._gram_track_cache_config` from `SPACEAGORA_GRAM_*` variables and embedded in `CallbackEnvConfig.gram_track_cache`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `mode` | Symbol | n/a | yes | Field `mode`. |
| in | `entry_horizon_s` | Float64 | n/a | yes | Field `entry_horizon_s`. |
| in | `entry_alt_tol_m` | Float64 | n/a | yes | Field `entry_alt_tol_m`. |
| in | `entry_ang_tol_rad` | Float64 | n/a | yes | Field `entry_ang_tol_rad`. |
| in | `entry_points` | Int | n/a | yes | Field `entry_points`. |
| in | `orbit_horizon_s` | Float64 | n/a | yes | Field `orbit_horizon_s`. |
| in | `orbit_alt_tol_m` | Float64 | n/a | yes | Field `orbit_alt_tol_m`. |
| in | `orbit_ang_tol_rad` | Float64 | n/a | yes | Field `orbit_ang_tol_rad`. |
| in | `orbit_points` | Int | n/a | yes | Field `orbit_points`. |
| in | `transition_band_m` | Float64 | n/a | yes | Field `transition_band_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | GramTrackCacheConfig | n/a | — | Constructed `GramTrackCacheConfig`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `callers` · call · `src/simulation/callbacks/gram_track_cache/config.jl:145-145`
- [[simulation_a.gram_cache_config_gram_track_cache_config|_gram_track_cache_config]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/config.jl:145-145`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No invariants are enforced: a negative horizon or zero `entry_points` is accepted and will fail downstream. The `mode` symbol's legal values are defined only in the builder. Because the struct is immutable and captured at run start, tolerances cannot adapt to observed prediction error.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 590.
