---
id: analysis.telemetry_loading__extract_extrema_series
label: _extract_extrema_series
kind: function
source:
  file: src/analysis/verification/telemetry_verification/telemetry_loading.jl
  symbol: _extract_extrema_series
  lines:
  - 72
  - 72
inputs:
- id: df
  type: DataFrame
  units: n/a
  required: true
  description: Positional argument `df`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: altitude_mode
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `altitude_mode`.
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
  description: Return value of `_extract_extrema_series`. Returns `(`.
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

# _extract_extrema_series

## Purpose
Derives per-orbit periapsis and apoapsis events (time, altitude in km, speed in km/s) from a simulated single-satellite state history stored in a `DataFrame`, using radial-velocity sign changes with a local-extrema fallback, so simulation output can be compared against telemetry-derived extrema.

## Theory & Math
Event time by linear interpolation of the radial-velocity zero: $$t_{evt} = (1-w)\,t_i + w\,t_{i+1}, \qquad w = \frac{|s_i|}{|s_i| + |s_{i+1}|}, \qquad s_i = \frac{\mathbf r_i \cdot \mathbf v_i}{|\mathbf r_i|}$$

## Design & Implementation
Reads `sc1_pos_*`/`sc1_position_*` and `sc1_vel_*`/`sc1_velocity_*` (metres, m/s) through `_require_column` and `df.time` (s); requires `n >= 3`. Altitude is `(r - planet.Rp_e)*1e-3` for `altitude_mode == :vacuum`, or geodetic via `r_intor_p!` and `rtolatlong` for `:oblate`; any other mode throws. Radial velocity `rdot = (r·v)/|r|` is scanned pairwise: a negative-to-non-negative crossing is a periapsis, positive-to-non-positive an apoapsis. The event is placed by weighting adjacent samples with `w = |s0|/(|s0|+|s1|)` (0.5 if both zero) and the same weight interpolates altitude and speed. Events closer than `min_sep_s = 500.0` s to the previous one of the same kind replace it only if more extreme. If either list is empty a second pass detects strict local minima/maxima of `alt_km` directly. Returns nested named tuples `peri` and `apo` with `orbit = 1.0:count`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `df` | DataFrame | n/a | yes | Positional argument `df`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `altitude_mode` | Symbol | n/a | yes | Positional argument `altitude_mode`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_extract_extrema_series`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.error_tables__time_aligned_rows_errors|_time_aligned_rows_errors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:91-91`
- [[envana.ana_error_tables_orbit_rows_errors|_orbit_rows_errors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:8-8`

**Downstream**

- `callees` → [[analysis.telemetry_loading__require_column|_require_column]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:73-73`
- `callees` → [[analysis.telemetry_loading__to_float_vector|_to_float_vector]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:79-79`
- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:95-95`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:96-96`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:128-128`
<!-- vulcan:connections:end -->

## Limitations
Only satellite 1 is processed; the column names are hard-coded with the `sc1_` prefix. The 500 s de-duplication window is fixed and unsuitable for orbits with periods below about 1000 s. The oblate branch calls `r_intor_p!` with no ephemeris time, so it relies on that method's default epoch handling. Non-finite `rdot` pairs are skipped, but a non-finite altitude still passes into the fallback comparison. The fallback uses non-strict comparisons on one side (`a1 <= a0`), so flat plateaus register an extremum at their last sample.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/telemetry_loading.jl` line 72.
