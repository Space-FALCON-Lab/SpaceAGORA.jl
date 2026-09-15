---
id: analysis.decay_diagnostics_flight_density_table
label: flight_density_table
kind: function
source:
  file: src/analysis/verification/telemetry_verification/decay_diagnostics.jl
  symbol: flight_density_table
  lines:
  - 119
  - 119
inputs:
- id: t_s
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `t_s`.
- id: sma_m
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `sma_m`.
- id: t_ref_s
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `t_ref_s`.
- id: sma_ref_m
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `sma_ref_m`.
- id: mu
  type: Real
  units: n/a
  required: true
  description: Keyword argument `mu`.
- id: cd_area_m2
  type: Real
  units: n/a
  required: true
  description: Keyword argument `cd_area_m2`.
- id: mass_kg
  type: Real
  units: n/a
  required: true
  description: Keyword argument `mass_kg`.
- id: period_s
  type: Real
  units: n/a
  required: true
  description: Keyword argument `period_s`.
- id: window_s
  type: Real
  units: n/a
  required: false
  description: Keyword argument `window_s` (default `21600.0`).
- id: step_s
  type: Real
  units: n/a
  required: false
  description: Keyword argument `step_s` (default `3600.0`).
- id: n_harmonics
  type: Int
  units: n/a
  required: false
  description: Keyword argument `n_harmonics` (default `3`).
- id: drifting_amplitudes
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `drifting_amplitudes` (default `false`).
- id: kwargs
  type: Vararg{Any}
  units: n/a
  required: false
  description: Keyword argument `kwargs` (variadic).
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
  type: DataFrame
  units: n/a
  description: Return value of `flight_density_table`.
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

# flight_density_table

## Purpose

`flight_density_table` infers an along-track atmospheric density history directly from flight telemetry. It slides a window over the measured semi-major-axis series, converts each window's zero-referenced decay rate into a density through the near-circular drag relation, and returns a `DataFrame(time_s, rho_kgm3)` in exactly the `tabulated_time` scenario-source format so the result can be replayed as a digital twin.

## Design & Implementation

After rejecting non-positive `cd_area_m2`, `mass_kg`, `window_s` and `step_s` with `ArgumentError`, the function walks window centres from `minimum(t) + window_s/2` in `step_s` increments while the centre stays within half a window of `maximum(t)`. For each centre it masks both the measurement and reference series to the half-open interval `[lo, hi)`, requires both counts to exceed `min_pts = 2 + 2*n_harmonics*(drifting_amplitudes ? 2 : 1) + 1`, and calls `zero_referenced_decay` to get `decay_m_per_day`, converted to metres per second by dividing by 86400. The window's mean SMA feeds the drag inversion, and only finite, strictly positive densities are pushed. If no window yields a valid density the function throws `ArgumentError`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `t_s` | AbstractVector{<:Real} | n/a | yes | Positional argument `t_s`. |
| in | `sma_m` | AbstractVector{<:Real} | n/a | yes | Positional argument `sma_m`. |
| in | `t_ref_s` | AbstractVector{<:Real} | n/a | yes | Positional argument `t_ref_s`. |
| in | `sma_ref_m` | AbstractVector{<:Real} | n/a | yes | Positional argument `sma_ref_m`. |
| in | `mu` | Real | n/a | yes | Keyword argument `mu`. |
| in | `cd_area_m2` | Real | n/a | yes | Keyword argument `cd_area_m2`. |
| in | `mass_kg` | Real | n/a | yes | Keyword argument `mass_kg`. |
| in | `period_s` | Real | n/a | yes | Keyword argument `period_s`. |
| in | `window_s` | Real | n/a | no | Keyword argument `window_s` (default `21600.0`). |
| in | `step_s` | Real | n/a | no | Keyword argument `step_s` (default `3600.0`). |
| in | `n_harmonics` | Int | n/a | no | Keyword argument `n_harmonics` (default `3`). |
| in | `drifting_amplitudes` | Bool | n/a | no | Keyword argument `drifting_amplitudes` (default `false`). |
| in | `kwargs` | Vararg{Any} | n/a | no | Keyword argument `kwargs` (variadic). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | DataFrame | n/a | — | Return value of `flight_density_table`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.decay_diagnostics_zero_referenced_decay|zero_referenced_decay]] · `callees` → `callers` · feedback · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl:104-104`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl`

**Downstream**

- `callees` → [[analysis.decay_diagnostics_zero_referenced_decay|zero_referenced_decay]] · `callers` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl:152-152`
- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl:158-158`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl:160-160`
<!-- vulcan:connections:end -->

## Theory & Math
The near-circular drag decay relation used is

$$\frac{da}{dt} = -\rho \frac{C_d A}{m}\sqrt{\mu a},$$

inverted per window as

$$\rho = \frac{-\dot a\, m}{C_d A \sqrt{\mu \bar a}},$$

with $\bar a$ the mean SMA over the window. The inferred $\rho$ therefore absorbs the supplied effective $C_d A$: a scale error in $C_d A$ rescales the whole profile uniformly without changing its shape.

## Limitations

The inversion assumes a near-circular orbit and that drag is the only non-conservative secular term, so solar radiation pressure and thrust events contaminate the result. Windows producing a negative or non-finite density are silently dropped, which can leave gaps or a heavily thinned table without warning. Density is only identifiable up to the assumed `cd_area_m2`, and overlapping windows (when `step_s < window_s`) make successive rows statistically correlated rather than independent samples.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/decay_diagnostics.jl` line 119.
