---
id: envana.ana_manifest_parsing_parse_cli
label: parse_cli
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: parse_cli
  lines:
  - 712
  - 756
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: TelemetryVerification namespace supplying default paths and the StudyConfig
    type.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: study_config
  type: StudyConfig
  units: n/a
  description: Resolved study configuration naming profile, manifest, output paths,
    enforcement and scenario filter.
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
# parse_cli

## Purpose
`parse_cli` resolves the settings for a verification study from command-line arguments and environment variables, producing the `StudyConfig` that the runner consumes.

## Theory & Math
Resolution is a strict precedence fold over three layers: an explicit command-line argument overrides an environment variable, which in turn overrides a compiled-in default. Each setting is read through `get(ENV, key, default)`, so the environment layer is consulted exactly once per field and the default is evaluated as the fallback. The profile field is normalised by `Symbol(lowercase(...))`, making `Quick`, `QUICK`, and `quick` all resolve to the same `:quick` symbol and removing case as a source of silent misconfiguration.

## Model & Assumptions
Recognised environment variables are `SPACEAGORA_TELEMETRY_PROFILE`, which defaults to `"quick"`, `SPACEAGORA_TELEMETRY_OUT_SUMMARY`, defaulting to `telemetry_orbit_accuracy_summary.csv` under `DEFAULT_OUTPUT_DIR`, `SPACEAGORA_TELEMETRY_OUT_ERRORS`, defaulting to `telemetry_orbit_accuracy_errors.csv` in the same directory, and `SPACEAGORA_TELEMETRY_MANIFEST`, defaulting to `DEFAULT_MANIFEST_PATH`. The function assumes those module constants have already been computed from `REPO_ROOT`, so it inherits the module's assumption about its own location on disk.

## Design & Implementation
The declared return type `StudyConfig` makes the constructor call at the end of the body the single exit point, and the fields `enforce`, `generate_plots`, and `scenarios` are passed as keywords so the record reads the same as its definition in `types.jl`. Because the scenario filter is a `Vector{String}` that defaults to empty, an unfiltered study and a study filtered to zero scenarios are distinguished by that emptiness rather than by a nullable field. The `@inline _study_config(request)` helper immediately below builds the same record from a programmatic `VerificationRequest`, giving the CLI and library paths a common target type.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | TelemetryVerification namespace supplying default paths and the StudyConfig type. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `study_config` | StudyConfig | n/a | — | Resolved study configuration naming profile, manifest, output paths, enforcement and scenario filter. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner_run_verification_cli|run_verification_cli]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:489-489`

**Downstream**

- `callees` → [[analysis.manifest_parsing__safe_parse_bool|_safe_parse_bool]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:717-717`
- `callees` → [[analysis.types__parse_scenario_list|_parse_scenario_list]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:719-719`
- `callees` → [[envana.ana_types_studyconfig|StudyConfig]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:747-747`
<!-- vulcan:connections:end -->

## Limitations
Unknown environment variable values are not validated, so an unrecognised profile symbol propagates until a tolerance lookup fails deep in reporting. Output paths are accepted verbatim without checking that the parent directory exists or is writable. There is no help text or usage message emitted from this function, so an operator learns about a malformed invocation only through a later exception.

## Provenance
Read directly from `src/analysis/verification/telemetry_verification/manifest_parsing.jl:712-756`, including the four environment lookups and the terminating `StudyConfig` construction.
