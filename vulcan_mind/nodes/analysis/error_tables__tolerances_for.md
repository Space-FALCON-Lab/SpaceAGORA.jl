---
id: analysis.error_tables__tolerances_for
label: _tolerances_for
kind: function
source:
  file: src/analysis/verification/telemetry_verification/error_tables.jl
  symbol: _tolerances_for
  lines:
  - 236
  - 236
inputs:
- id: cfg
  type: OrbitEventsScenarioConfig
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
Selects which tolerance set an orbit-events scenario is scored against, according to whether the verification run is a quick smoke check or a full comparison.

## Design & Implementation
A one-line `@inline` definition returning `cfg.tolerances_quick` when `profile` is `:quick` and `cfg.tolerances_full` for every other symbol. Keeping the choice in one function means the quick and full profiles cannot drift apart across the several call sites that score channels, and the scenario config carries both sets rather than being mutated between runs.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | OrbitEventsScenarioConfig | n/a | yes | Positional argument `cfg`. |
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
Any profile symbol other than `:quick` silently falls through to the full tolerances, so a misspelled profile name is scored strictly rather than rejected.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/error_tables.jl` line 236.
