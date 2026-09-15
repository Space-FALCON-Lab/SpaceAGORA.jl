---
id: analysis.reporting__tolerances_for
label: _tolerances_for
kind: function
source:
  file: src/analysis/verification/telemetry_verification/reporting.jl
  symbol: _tolerances_for
  lines:
  - 1
  - 1
inputs:
- id: cfg
  type: TimeAlignedScenarioConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: profile
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `profile`.
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
  description: 'Return value of `_tolerances_for`. Returns `profile == :quick ? cfg.tolerances_quick
    : cfg.tolerances_full`.'
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

# _tolerances_for

## Purpose
Selects which tolerance table applies to a telemetry-verification run: the `:quick` profile reads `cfg.tolerances_quick`, any other profile symbol reads `cfg.tolerances_full`. It is the single switch point that `_evaluate_thresholds` uses to look up per-event `max_abs_km`, `max_nmae` and `max_rmse_km` limits.

## Design & Implementation
An `@inline` one-liner taking `cfg::TimeAlignedScenarioConfig` and `profile::Symbol`. The comparison is `profile == :quick`; no other profile names are recognised, so `:full`, `:ci`, or a typo all resolve to the full tolerance map. Nothing is mutated and nothing is thrown; the returned object is whatever container the config stores (a map from event name to a tolerance tuple).

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | TimeAlignedScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `profile` | Symbol | n/a | yes | Positional argument `profile`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_tolerances_for`. Returns `profile == :quick ? cfg.tolerances_quick : cfg.tolerances_full`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[envana.ana_reporting_evaluate_thresholds|_evaluate_thresholds]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:145-145`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only a `TimeAlignedScenarioConfig` method exists even though `_evaluate_thresholds` accepts `AbstractScenarioConfig`; passing an `OrbitEventsScenarioConfig` hits a `MethodError` unless another file adds that method. Unknown profile symbols silently fall back to the full tolerances rather than raising, which can mask a misspelt CLI flag.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/reporting.jl` line 1.
