---
id: envana.ana_reporting_evaluate_thresholds
label: _evaluate_thresholds
kind: function
source:
  file: src/analysis/verification/telemetry_verification/reporting.jl
  symbol: _evaluate_thresholds
  lines:
  - 144
  - 171
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: TelemetryVerification namespace supplying per-profile tolerance maps
    and minimum sample counts.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: verdict
  type: NamedTuple
  units: km
  description: Overall pass flag, the four component flags, and the three numeric
    limits that were applied.
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
# _evaluate_thresholds

## Purpose
`_evaluate_thresholds` turns one summary row into a pass or fail verdict by comparing its accuracy metrics against the tolerances configured for the scenario and profile, and reports which individual criterion failed.

## Theory & Math
The verdict is the conjunction of four independent tests: `n_sim >= n_min`, `max_abs_km <= tol.max_abs_km`, `nmae <= tol.max_nmae`, and `rmse_km <= tol.max_rmse_km`. Two of these are absolute in kilometres and one, `nmae`, is dimensionless because it is normalised by the telemetry peak-to-peak range, so a scenario must satisfy both an absolute accuracy floor and a relative shape criterion. The sample-count test guards the other three, since a mean or maximum computed over too few points is not a meaningful estimate of accuracy.

## Model & Assumptions
Tolerances come from `_tolerances_for(cfg, profile)`, so the same scenario can be judged more loosely under the quick profile than under the strict one. Speed-derived events fall back to their base event's tolerance: when the event name ends in `_speed` and has no entry of its own, the suffix is stripped and the lookup repeated, which lets a manifest specify one tolerance for an altitude curve and inherit it for the matching rate curve. A completely missing tolerance is treated as a configuration error and raises an `ArgumentError` naming both the event and the scenario.

## Design & Implementation
The returned `NamedTuple` carries the overall `pass` flag, the four component flags `pass_points`, `pass_abs`, `pass_nmae`, and `pass_rmse`, and the three applied limits plus `min_eval_points`. Echoing the limits alongside the flags means the emitted report is self-describing: a reader can see the threshold that was applied without re-reading the manifest. Comparisons use non-strict inequalities, so a metric exactly at its limit passes.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | TelemetryVerification namespace supplying per-profile tolerance maps and minimum sample counts. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `verdict` | NamedTuple | km | — | Overall pass flag, the four component flags, and the three numeric limits that were applied. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:382-382`

**Downstream**

- `callees` → [[analysis.error_tables__tolerances_for|_tolerances_for]] · `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:145-145`
- `callees` → [[analysis.reporting__min_eval_points|_min_eval_points]] · `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:154-154`
- `callees` → [[analysis.reporting__tolerances_for|_tolerances_for]] · `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:145-145`
<!-- vulcan:connections:end -->

## Limitations
Inheriting a `_speed` event's tolerance from its base event applies a kilometre-scale altitude limit to a kilometre-per-second rate quantity, which is only sensible when the manifest was written with that in mind. The four criteria are weighted equally with no notion of margin, so a run that fails one criterion by a fraction of a percent is reported identically to one that fails by an order of magnitude. Metrics of `Inf` from an empty simulation fail every numeric test but produce no distinct diagnostic.

## Provenance
Read directly from `src/analysis/verification/telemetry_verification/reporting.jl:144-171`, including the `_speed` fallback lookup and the four-way conjunction.
