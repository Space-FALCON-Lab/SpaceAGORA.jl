---
id: analysis.comparison_metrics__apo_decay_diagnostic
label: _apo_decay_diagnostic
kind: function
source:
  file: src/analysis/verification/telemetry_verification/comparison_metrics.jl
  symbol: _apo_decay_diagnostic
  lines:
  - 46
  - 46
inputs:
- id: tele_orbit
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `tele_orbit`.
- id: tele_alt
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `tele_alt`.
- id: sim_axis
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `sim_axis`.
- id: sim_alt
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `sim_alt`.
- id: maneuver_orbits
  type: Vector{Float64}
  units: n/a
  required: true
  description: Positional argument `maneuver_orbits`.
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
  description: Return value of `_apo_decay_diagnostic`. Returns `(`.
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

# _apo_decay_diagnostic

## Purpose

Measures drag fidelity by comparing per-orbit apoapsis decay rates between telemetry and simulation, rather than comparing absolute apoapsis altitudes. Absolute apsis error accumulates any per-pass density bias over a whole campaign and inherits orbit-index misalignment once the period drifts; the rate ratio isolates per-pass drag directly.

## Design & Implementation

Given `tele_orbit`, `tele_alt`, `sim_axis`, `sim_alt` and `maneuver_orbits`, the function returns `_DECAY_DIAGNOSTIC_EMPTY` unless both series carry at least three samples. It forms telemetry and simulated decay rates with the local `_rates` closure, then builds a `keep` mask that drops telemetry intervals whose midpoint falls outside the simulated midpoint span and any interval bracketing a maneuver orbit, because the impulse rather than drag dominates those. The simulated rate is interpolated onto the kept telemetry midpoints with `_interp_linear`. Per-interval ratios are formed only where the telemetry rate magnitude exceeds `1e-9` km/orbit, and the aggregate weights each interval by its orbit span so that gaps in telemetry sampling do not count the same as single-orbit intervals. The named tuple returned carries `drag_decay_ratio_median`, `drag_decay_ratio_total` and `drag_decay_n`.

## Theory & Math

With telemetry rates $\dot{h}^{\text{tel}}_k$, interpolated simulated rates $\dot{h}^{\text{sim}}_k$ and orbit spans $s_k$ over the kept intervals $k \in K$, the diagnostic reports

$$\rho_{\text{med}} = \operatorname{median}_{k \in K,\; |\dot{h}^{\text{tel}}_k| > 10^{-9}} \frac{\dot{h}^{\text{sim}}_k}{\dot{h}^{\text{tel}}_k}$$

$$\rho_{\text{tot}} = \frac{\sum_{k \in K} \dot{h}^{\text{sim}}_k s_k}{\sum_{k \in K} \dot{h}^{\text{tel}}_k s_k}$$

A value of $1$ means the simulated drag removes energy at the same rate as flight; values above $1$ mean the model is too draggy. The span weighting makes $\rho_{\text{tot}}$ the ratio of total altitude lost over the overlapping campaign.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tele_orbit` | Vector{Float64} | n/a | yes | Positional argument `tele_orbit`. |
| in | `tele_alt` | Vector{Float64} | n/a | yes | Positional argument `tele_alt`. |
| in | `sim_axis` | Vector{Float64} | n/a | yes | Positional argument `sim_axis`. |
| in | `sim_alt` | Vector{Float64} | n/a | yes | Positional argument `sim_alt`. |
| in | `maneuver_orbits` | Vector{Float64} | n/a | yes | Positional argument `maneuver_orbits`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_apo_decay_diagnostic`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.error_tables__time_aligned_rows_errors|_time_aligned_rows_errors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:112-112`
- [[envana.ana_error_tables_orbit_rows_errors|_orbit_rows_errors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:50-50`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

Both series must be sampled on a common orbit-index abscissa; no cross-check confirms that orbit numbering is aligned between telemetry and simulation. Telemetry altitude quantisation makes individual ratios noisy, which is why the median is reported alongside the aggregate, but the median is still undefined when every telemetry rate falls under the `1e-9` threshold and `NaN` is returned instead. The maneuver exclusion test uses closed-interval containment on `tele_orbit[k] <= m <= tele_orbit[k+1]`, so a maneuver exactly on a sample boundary removes two intervals. A near-zero denominator in the aggregate produces `NaN` rather than a flagged error, and `drag_decay_n` counts only the ratio samples, not the intervals used in the aggregate.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/comparison_metrics.jl` line 46.
