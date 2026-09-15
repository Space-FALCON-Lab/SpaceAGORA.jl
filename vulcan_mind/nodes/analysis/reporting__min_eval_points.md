---
id: analysis.reporting__min_eval_points
label: _min_eval_points
kind: function
source:
  file: src/analysis/verification/telemetry_verification/reporting.jl
  symbol: _min_eval_points
  lines:
  - 131
  - 131
inputs:
- id: cfg
  type: OrbitEventsScenarioConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  description: Return value of `_min_eval_points`. Returns `cfg.min_eval_points`.
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

# _min_eval_points

## Purpose
Returns the minimum number of simulated evaluation points an event needs before its accuracy metrics are considered trustworthy; `_evaluate_thresholds` compares `row.n_sim` against this value to set `pass_points`.

## Design & Implementation
Two `@inline` methods, for `OrbitEventsScenarioConfig` and `TimeAlignedScenarioConfig`, each returning the `min_eval_points` field of the config directly. The value is also copied into the threshold result tuple as `min_eval_points` so the summary CSV records the limit applied.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | OrbitEventsScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_min_eval_points`. Returns `cfg.min_eval_points`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[envana.ana_reporting_evaluate_thresholds|_evaluate_thresholds]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:154-154`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No lower bound is enforced; a config with `min_eval_points = 0` makes the point-count check vacuous. The type of the field is whatever the config parser produced, so a `Float64` value would still compare but would be serialised inconsistently.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/reporting.jl` line 131.
