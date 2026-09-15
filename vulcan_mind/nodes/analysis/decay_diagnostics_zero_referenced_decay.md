---
id: analysis.decay_diagnostics_zero_referenced_decay
label: zero_referenced_decay
kind: function
source:
  file: src/analysis/verification/telemetry_verification/decay_diagnostics.jl
  symbol: zero_referenced_decay
  lines:
  - 90
  - 90
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
- id: period_s
  type: Real
  units: n/a
  required: true
  description: Keyword argument `period_s`.
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
  type: Tuple
  units: n/a
  description: Return value of `zero_referenced_decay`. Returns `(decay_m_per_day=raw
    - ref, raw_m_per_day=raw, reference_m_per_day=ref)`.
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

# zero_referenced_decay

## Purpose

`zero_referenced_decay` produces the physically meaningful secular orbit decay rate of a measured semi-major-axis series by removing the estimator's own window leakage. It returns the named tuple `(decay_m_per_day, raw_m_per_day, reference_m_per_day)` in metres per day.

## Design & Implementation

The function applies `secular_sma_slope` twice with identical settings: once to the measurement pair `(t_s, sma_m)` and once to a drag-free reference propagation `(t_ref_s, sma_ref_m)` of the same arc, forwarding `period_s` and all remaining `kwargs...` to both calls so the harmonic order and amplitude-drift options match exactly. The reference slope is then subtracted from the raw slope. Since a conservative gravity field cannot change SMA secularly, the reference slope is by construction pure estimator and window leakage, and the subtraction cancels the common-mode J2 apsidal-precession term the module header measures at 11 to 12 metres per day on a drag-free 48-hour LEO arc.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `t_s` | AbstractVector{<:Real} | n/a | yes | Positional argument `t_s`. |
| in | `sma_m` | AbstractVector{<:Real} | n/a | yes | Positional argument `sma_m`. |
| in | `t_ref_s` | AbstractVector{<:Real} | n/a | yes | Positional argument `t_ref_s`. |
| in | `sma_ref_m` | AbstractVector{<:Real} | n/a | yes | Positional argument `sma_ref_m`. |
| in | `period_s` | Real | n/a | yes | Keyword argument `period_s`. |
| in | `kwargs` | Vararg{Any} | n/a | no | Keyword argument `kwargs` (variadic). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple | n/a | — | Return value of `zero_referenced_decay`. Returns `(decay_m_per_day=raw - ref, raw_m_per_day=raw, reference_m_per_day=ref)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.decay_diagnostics_flight_density_table|flight_density_table]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl:152-152`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl`

**Downstream**

- `callees` → [[analysis.decay_diagnostics_flight_density_table|flight_density_table]] · `callers` · feedback · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl:104-104`
- `callees` → [[envana.ana_decay_diagnostics_secular_sma_slope|secular_sma_slope]] · `callers` · call · `src/analysis/verification/telemetry_verification/decay_diagnostics.jl:98-98`
<!-- vulcan:connections:end -->

## Theory & Math
$$\dot a_{\text{decay}} = \dot a_{\text{raw}} - \dot a_{\text{ref}},$$

where each slope is the linear coefficient of the joint fit $a(t) \approx a_0 + \dot a\, t + \sum_{k=1}^{N} \left(s_k \sin k\omega t + c_k \cos k\omega t\right)$ with $\omega = 2\pi/T$ at the orbital period $T$.

## Limitations

The cancellation only holds when the reference propagation shares the measurement's gravity field, time window and gap structure; a mismatched reference leaves residual leakage in the answer. The two series are not required to be sampled on the same grid, and nothing checks that the reference is genuinely drag-free. All validation of vector lengths, sample count and positive `period_s` is delegated to `secular_sma_slope`, which throws `ArgumentError` on violation.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/decay_diagnostics.jl` line 90.
