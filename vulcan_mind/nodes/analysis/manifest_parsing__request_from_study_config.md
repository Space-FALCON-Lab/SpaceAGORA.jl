---
id: analysis.manifest_parsing__request_from_study_config
label: _request_from_study_config
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _request_from_study_config
  lines:
  - 772
  - 772
inputs:
- id: cfg
  type: StudyConfig
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
  type: VerificationRequest
  units: n/a
  description: Return value of `_request_from_study_config`.
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

# _request_from_study_config

## Purpose
The inverse conversion, rebuilding a `VerificationRequest` from a `StudyConfig` so a CLI-parsed run can be re-issued programmatically.

## Design & Implementation
Copies every field by name into the keyword constructor, copying the scenario list. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | StudyConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | VerificationRequest | n/a | — | Return value of `_request_from_study_config`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner_run_verification_cli|run_verification_cli]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:490-490`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.types_verificationrequest|VerificationRequest]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:773-773`
<!-- vulcan:connections:end -->

## Limitations
Round-tripping is lossy only in that absolute paths from the config are kept absolute, which is fine but means the request no longer reflects what the user typed.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 772.
