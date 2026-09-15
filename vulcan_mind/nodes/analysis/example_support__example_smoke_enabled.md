---
id: analysis.example_support__example_smoke_enabled
label: _example_smoke_enabled
kind: function
source:
  file: src/analysis/verification/telemetry_verification/example_support.jl
  symbol: _example_smoke_enabled
  lines:
  - 3
  - 3
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
  description: Return value of `_example_smoke_enabled`. Returns `get(ENV, "SPACEAGORA_EXAMPLE_SMOKE",
    "0") == "1"`.
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

# _example_smoke_enabled

## Purpose
Reports whether the example scripts are running in smoke-test mode, where mission length and output are cut down so a CI job can exercise every example quickly.

## Design & Implementation
An `@inline` one-liner reading the `SPACEAGORA_EXAMPLE_SMOKE` environment variable through `get` with a default of `"0"` and comparing it to the literal string `"1"`. Reading the environment on every call rather than caching at load means a test harness can flip the mode between examples in one process without reloading the package.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_example_smoke_enabled`. Returns `get(ENV, "SPACEAGORA_EXAMPLE_SMOKE", "0") == "1"`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.example_support__example_smoke_args|_example_smoke_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:20-20`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/example_support.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the exact string `1` enables the mode; `true`, `yes` or `on` are treated as disabled with no warning, so a misconfigured CI variable silently runs the full-length examples.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/example_support.jl` line 3.
