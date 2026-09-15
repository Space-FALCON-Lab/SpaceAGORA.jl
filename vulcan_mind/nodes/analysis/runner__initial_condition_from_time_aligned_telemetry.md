---
id: analysis.runner__initial_condition_from_time_aligned_telemetry
label: _initial_condition_from_time_aligned_telemetry
kind: function
source:
  file: src/analysis/verification/telemetry_verification/runner.jl
  symbol: _initial_condition_from_time_aligned_telemetry
  lines:
  - 100
  - 100
inputs:
- id: cfg
  type: TimeAlignedScenarioConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: telemetry
  type: Any
  units: n/a
  required: true
  description: Positional argument `telemetry`.
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
  type: Union{CartesianInitialCondition, InitialCondition}
  units: n/a
  description: Return value of `_initial_condition_from_time_aligned_telemetry`. Returns
    `CartesianInitialCondition(` or `InitialCondition(`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- analysis
origin: agent
---

# _initial_condition_from_time_aligned_telemetry

## Purpose
`_initial_condition_from_time_aligned_telemetry` builds the simulation initial condition for a time-aligned telemetry scenario, preferring a Cartesian state from the telemetry's `*_ic_km`/`*_ic_kmps` columns and falling back to Keplerian elements. It converts units, frames, and applies scenario-configured offsets.

## Design & Implementation
Signature `(cfg::TimeAlignedScenarioConfig, telemetry)`. If `telemetry.x_ic_km !== nothing` it assembles `r_m` and `v_mps` as `SVector{3,Float64}` scaled by `1e3` (km to m, km/s to m/s); when `cfg.cartesian_ic_frame == :planet_fixed` it computes `et0 = _initial_time_et(cfg.initial_time)` and rotates the state to J2000 via `_planet_fixed_to_j2000_state(cfg.planet_name, et0, r_m, v_mps)`. The offsets `cfg.ic_offset_m` and `cfg.ic_offset_mps` are added and a `CartesianInitialCondition` with plain `Vector` fields is returned. Without Cartesian columns, any non-zero offset throws `ArgumentError`; then all six Keplerian fields (`sma_km`, `ecc`, `inc_deg`, `aop_deg`, `raan_deg`, `ta_deg`) must be finite or another `ArgumentError` is thrown. Apoapsis and periapsis radii are `sma_m * (1 ± ecc)` and an `InitialCondition(ra, rp, i, ω, Ω, ν)` in degrees is returned.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | TimeAlignedScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `telemetry` | Any | n/a | yes | Positional argument `telemetry`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{CartesianInitialCondition, InitialCondition} | n/a | — | Return value of `_initial_condition_from_time_aligned_telemetry`. Returns `CartesianInitialCondition(` or `InitialCondition(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_single_scenario|_run_single_scenario]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:220-220`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/runner.jl`

**Downstream**

- `callees` → [[analysis.telemetry_loading__initial_time_et|_initial_time_et]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:105-105`
- `callees` → [[analysis.telemetry_loading__planet_fixed_to_j2000_state|_planet_fixed_to_j2000_state]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:106-106`
- `callees` → [[vehicle.model_cartesianinitialcondition|CartesianInitialCondition]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:108-108`
- `callees` → [[vehicle.model_initialcondition|InitialCondition]] · `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:124-124`
<!-- vulcan:connections:end -->

## Limitations
Presence of Cartesian ICs is inferred from `x_ic_km` alone; if `y`/`z` or velocity columns are missing the `SVector` constructor fails with a `MethodError` rather than a descriptive message. Hyperbolic telemetry (`ecc >= 1`) yields a negative `rp_m` with no check. Keplerian angles are passed through in degrees and assumed to match the `InitialCondition` convention. Only `:planet_fixed` triggers a frame conversion; any other symbol is treated as inertial without validation.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/runner.jl` line 100.
