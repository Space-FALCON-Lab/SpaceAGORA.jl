---
id: analysis.example_support__example_smoke_results_enabled
label: _example_smoke_results_enabled
kind: function
source:
  file: src/analysis/verification/telemetry_verification/example_support.jl
  symbol: _example_smoke_results_enabled
  lines:
  - 4
  - 4
inputs:
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
  description: Return value of `_example_smoke_results_enabled`. Returns `get(ENV,
    "SPACEAGORA_EXAMPLE_SMOKE_RESULTS", "0") == "1"`.
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

# _example_smoke_results_enabled

## Purpose
Decides whether a smoke-mode example run should still write its results files, which are normally suppressed to keep CI fast and clean.

## Design & Implementation
Mirrors the smoke-enabled predicate exactly: `@inline`, reads `SPACEAGORA_EXAMPLE_SMOKE_RESULTS` with default `"0"` and tests for `"1"`. Its value is consulted by `_example_smoke_args` to set both the `results` and `save_csv` flags of the trimmed simulation settings, so a single variable controls whether the smoke run leaves a CSV behind for a downstream verification step to read.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_example_smoke_results_enabled`. Returns `get(ENV, "SPACEAGORA_EXAMPLE_SMOKE_RESULTS", "0") == "1"`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.example_support__example_smoke_args|_example_smoke_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:26-26`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/example_support.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It is meaningful only when smoke mode is also enabled; outside smoke mode the ordinary settings apply and this variable is never consulted, which is easy to miss when a results file unexpectedly appears or does not.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/example_support.jl` line 4.
