---
id: envana.ana_decay_diagnostics_secular_sma_slope
label: secular_sma_slope
kind: function
source:
  file: src/analysis/verification/telemetry_verification/decay_diagnostics.jl
  symbol: secular_sma_slope
  lines:
  - 41
  - 77
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: TelemetryVerification namespace providing the linear algebra backslash
    solver.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: slope
  type: Float64
  units: m/day
  description: Secular semi-major axis decay rate extracted from an harmonic-plus-trend
    least squares fit.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- envana
origin: agent
---
# secular_sma_slope

## Purpose
`secular_sma_slope` extracts the secular decay rate of semi-major axis from a sampled time history, separating the long-term trend from the periodic oscillations that orbital perturbations impose at multiples of the orbital frequency.

## Theory & Math
The design matrix models `a(t) = c_1 * tau + c_2 + sum_{k=1}^{K} (alpha_k sin(k omega tau) + beta_k cos(k omega tau))`, where `tau = t - mean(t)` is centred time in seconds, `omega = 2 pi / P` with orbital period `P` in seconds, and `K` is `n_harmonics`. When `drifting_amplitudes` is set, each harmonic gains the two extra columns `tau sin(k omega tau)` and `tau cos(k omega tau)`, letting the oscillation envelope grow or shrink linearly. The fit is the ordinary least squares solution `c = argmin |A c - a|_2`, and the returned value is `c_1 * 86400.0`, converting the slope from metres per second to metres per day.

## Model & Assumptions
Centring time about its mean decorrelates the trend column from the constant column, which markedly improves the conditioning of `A`. The routine requires `n > 2 + 2 K (1 or 2)` samples so the system is overdetermined, and validates that the time and semi-major axis vectors have equal length, that `period_s > 0`, and that `n_harmonics >= 1`; each violation raises an `ArgumentError` naming the function and the offending quantity. Uniform sampling is not required, since the solve is a general least squares problem.

## Design & Implementation
Columns are accumulated into a `Vector{Vector{Float64}}` and combined with `reduce(hcat, cols)`, placing the trend column first so that `coef[1]` is always the secular rate regardless of harmonic count. The solve uses Julia's `\` operator, which selects a QR factorisation for the rectangular system rather than forming the normal equations, avoiding the squaring of the condition number. Defaults are three harmonics with fixed amplitudes.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | TelemetryVerification namespace providing the linear algebra backslash solver. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `slope` | Float64 | m/day | — | Secular semi-major axis decay rate extracted from an harmonic-plus-trend least squares fit. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.decay_diagnostics_visviva_sma|visviva_sma]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl:27-27`
- [[analysis.decay_diagnostics_zero_referenced_decay|zero_referenced_decay]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl:98-98`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl:64-64`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl:61-61`
<!-- vulcan:connections:end -->

## Limitations
A single orbital period drives every harmonic column, so a decaying orbit whose period shrinks measurably over the arc is modelled with a stale frequency and leaks power into the trend. The fit is unweighted, so noisier late samples carry the same influence as clean early ones. Fewer than the minimum required samples raises rather than degrading gracefully to a plain linear fit.

## Provenance
Read directly from `src/analysis/verification/telemetry_verification/decay_diagnostics.jl:41-77`, including the argument validation block, the harmonic column construction, and the `86400.0` unit conversion.
