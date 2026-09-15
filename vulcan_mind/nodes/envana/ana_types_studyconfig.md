---
id: envana.ana_types_studyconfig
label: StudyConfig
kind: struct
source:
  file: src/analysis/verification/telemetry_verification/types.jl
  symbol: StudyConfig
  lines:
  - 229
  - 237
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: TelemetryVerification namespace in which the configuration record types
    are defined.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: config_record
  type: StudyConfig
  units: n/a
  description: Immutable record of profile, output paths, manifest path, enforcement
    flag and scenario filter.
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
# StudyConfig

## Purpose
`StudyConfig` is the immutable record that fully describes one verification study run. It is the shared target type produced both by `parse_cli` from command-line and environment input and by `_study_config` from a programmatic `VerificationRequest`, so the runner has exactly one configuration shape to consume.

## Theory & Math
The record separates the two orthogonal axes that determine a run's cost and strictness. `profile` selects the numerical regime, choosing between the quick tolerances and the strict values of `1e-7` relative and `1e-9` absolute with orbit steps capped at 60.0 s and atmospheric steps at 0.2 s. `enforce` selects the consequence: when true, a failing threshold comparison aborts the study, and when false the metrics are still computed and written but the run completes. Together they span cheap exploratory runs, cheap gating runs, expensive exploratory runs, and expensive gating runs.

## Model & Assumptions
`profile` is a `Symbol` rather than an enumeration, so any symbol constructs successfully and an invalid value is caught only when tolerances are looked up. The three path fields `out_summary`, `out_errors`, and `manifest_path` are plain `String` values with no existence check at construction time. `scenarios` defaults to an empty `String[]`, and the inline comment states the convention explicitly: an empty vector means every scenario in the manifest, so an unfiltered study and a study filtered to nothing are the same value.

## Design & Implementation
The type is declared with `Base.@kwdef`, making every field keyword-constructible and letting the single defaulted field sit last without forcing positional ordering on the caller. Only `scenarios` carries a default; the remaining six fields are required, which prevents a study from silently running against an unintended manifest or profile. Being an immutable struct, it can be shared across the runner without any risk of mid-study mutation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | TelemetryVerification namespace in which the configuration record types are defined. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `config_record` | StudyConfig | n/a | — | Immutable record of profile, output paths, manifest path, enforcement flag and scenario filter. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__study_config|_study_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:761-761`
- [[envana.ana_manifest_parsing_parse_cli|parse_cli]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:747-747`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the empty vector means all scenarios, there is no way to express a deliberately empty study through this record. Path fields carry no validation, so a misspelled manifest path fails only when the runner opens it. The record carries no random seed or version stamp, so two studies with identical `StudyConfig` values are not guaranteed bit-identical if the underlying model changes.

## Provenance
Read directly from `src/analysis/verification/telemetry_verification/types.jl:229-237`, including the inline comment documenting the empty-scenarios convention.
