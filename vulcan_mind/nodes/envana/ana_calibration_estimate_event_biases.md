---
id: envana.ana_calibration_estimate_event_biases
label: _estimate_event_biases
kind: function
source:
  file: src/analysis/verification/telemetry_verification/calibration.jl
  symbol: _estimate_event_biases
  lines:
  - 38
  - 56
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: TelemetryVerification namespace supplying the error tables and the
    float coercion utility.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: bias_map
  type: Dict{String,Float64}
  units: km
  description: Per-event additive altitude bias, or zero where the estimate saturates
    the configured cap.
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
# _estimate_event_biases

## Purpose
`_estimate_event_biases` derives a per-event additive altitude correction from the earliest part of each error table, so that a constant offset between simulation and telemetry does not contaminate the accuracy metrics that follow.

## Theory & Math
For each event the routine takes the first `k = min(_BIAS_ESTIMATE_POINTS, n)` error samples and computes `b = -median(e_1..e_k)`, in kilometres, where `e_j` is the signed simulated-minus-telemetry altitude error of sample `j`. The median is used rather than the mean because its breakdown point is 50 percent, so a single outlying early orbit cannot move the estimate; the mean has a breakdown point of zero. Applying `b` shifts the simulated curve so that the early-arc median residual is zero, which leaves the shape metrics unchanged while removing the constant term.

## Model & Assumptions
The estimate assumes the leading `k` samples are representative of a genuinely constant offset, which holds when the initial state is slightly wrong but the force model is right. If the absolute value of `b` reaches `bias_abs_max_km`, that assumption is judged false: the code prints a `calibration_bias_saturated` line naming the event, the rounded bias, and the cap, then sets `b = 0.0` so the uncorrected error is reported. Empty tables are skipped with `nrow(df) == 0 && continue`, leaving the event absent from the returned dictionary.

## Design & Implementation
Errors are pulled through `_to_float_vector(df.error_km, "error_km:$event")`, which both coerces the column and names the source in any failure message. The event key comes from the first row of the table, so each `DataFrame` in the input vector must hold exactly one event. The result is a plain `Dict{String,Float64}` keyed by event name, which downstream error tabulation adds to the simulated series before comparison.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | TelemetryVerification namespace supplying the error tables and the float coercion utility. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `bias_map` | Dict{String,Float64} | km | — | Per-event additive altitude bias, or zero where the estimate saturates the configured cap. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_single_scenario|_run_single_scenario]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:168-168`

**Downstream**

- `callees` → [[analysis.telemetry_loading__to_float_vector|_to_float_vector]] · `callers` · call · `src/analysis/verification/telemetry_verification/calibration.jl:43-43`
- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/analysis/verification/telemetry_verification/calibration.jl:47-47`
<!-- vulcan:connections:end -->

## Limitations
Only a constant offset is estimated; a linear or periodic bias passes through untouched and inflates the reported error. The saturation cap converts a large genuine offset into zero correction rather than raising, so a caller that ignores the printed line sees a silently uncalibrated event. Events absent from the returned dictionary receive no correction at all.

## Provenance
Read directly from `src/analysis/verification/telemetry_verification/calibration.jl:38-56`, including the median estimator, the saturation branch, and the diagnostic print.
