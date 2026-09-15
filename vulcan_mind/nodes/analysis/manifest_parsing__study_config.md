---
id: analysis.manifest_parsing__study_config
label: _study_config
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _study_config
  lines:
  - 758
  - 758
inputs:
- id: request
  type: VerificationRequest
  units: n/a
  required: true
  description: Positional argument `request`.
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
  type: StudyConfig
  units: n/a
  description: Return value of `_study_config`.
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

# _study_config

## Purpose
Converts a programmatic `VerificationRequest` into the `StudyConfig` the runner consumes, applying the same validation as the CLI path.

## Design & Implementation
Lowercases the profile into a symbol, requires `:quick` or `:full`, and copies the remaining fields with output and manifest paths made absolute and the scenario list copied. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `request` | VerificationRequest | n/a | yes | Positional argument `request`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | StudyConfig | n/a | — | Return value of `_study_config`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.run_verification|run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:479-479`

**Downstream**

- `callees` → [[envana.ana_types_studyconfig|StudyConfig]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:761-761`
<!-- vulcan:connections:end -->

## Limitations
It duplicates the profile validation of `parse_cli`; the two must be kept in step.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 758.
