---
id: analysis.types_verificationrequest
label: VerificationRequest
kind: struct
source:
  file: src/analysis/verification/telemetry_verification/types.jl
  symbol: VerificationRequest
  lines:
  - 256
  - 256
inputs:
- id: profile
  type: Symbol
  units: n/a
  required: false
  description: Field `profile` (default `Symbol(lowercase(get(ENV, "SPACEAGORA_TELEMETRY_PROFILE",
    "quick")))`).
- id: out_summary
  type: String
  units: n/a
  required: false
  description: Field `out_summary` (default `abspath(get(ENV, "SPACEAGORA_TELEMETRY_OUT_SUMMARY",
    joinpath(DEFAULT_OUTPUT_DIR, "telemetry_orbit_accuracy_summary.csv")))`).
- id: out_errors
  type: String
  units: n/a
  required: false
  description: Field `out_errors` (default `abspath(get(ENV, "SPACEAGORA_TELEMETRY_OUT_ERRORS",
    joinpath(DEFAULT_OUTPUT_DIR, "telemetry_orbit_accuracy_errors.csv")))`).
- id: manifest_path
  type: String
  units: n/a
  required: false
  description: Field `manifest_path` (default `abspath(get(ENV, "SPACEAGORA_TELEMETRY_MANIFEST",
    DEFAULT_MANIFEST_PATH))`).
- id: enforce
  type: Bool
  units: n/a
  required: false
  description: Field `enforce` (default `false`).
- id: generate_plots
  type: Bool
  units: n/a
  required: false
  description: Field `generate_plots` (default `_safe_parse_bool(get(ENV, "SPACEAGORA_TELEMETRY_PLOTS",
    "1"), true)`).
- id: scenarios
  type: Vector{String}
  units: n/a
  required: false
  description: Field `scenarios` (default `String[]`).
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
  description: Constructed `VerificationRequest` (keyword constructor via @kwdef).
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

# VerificationRequest

## Purpose
Typed request object for the telemetry verification study entry point, capturing the study profile, output CSV paths, manifest path, enforcement flag, plot generation, and an optional scenario subset.

## Design & Implementation
`Base.@kwdef struct` whose defaults are resolved from environment variables at construction time: `profile = Symbol(lowercase(get(ENV, "SPACEAGORA_TELEMETRY_PROFILE", "quick")))`, `out_summary` and `out_errors` from `SPACEAGORA_TELEMETRY_OUT_SUMMARY`/`_OUT_ERRORS` (defaulting under `DEFAULT_OUTPUT_DIR`), `manifest_path` from `SPACEAGORA_TELEMETRY_MANIFEST` or `DEFAULT_MANIFEST_PATH`, all passed through `abspath`. `enforce=false`, `generate_plots` parsed with `_safe_parse_bool` from `SPACEAGORA_TELEMETRY_PLOTS` (default true), and `scenarios::Vector{String}=String[]` meaning every scenario. The `SPACEAGORA_TELEMETRY_SCENARIOS` filter is deliberately applied only by `parse_cli`, so requests built in code are never narrowed by an environment variable they did not request.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `profile` | Symbol | n/a | no | Field `profile` (default `Symbol(lowercase(get(ENV, "SPACEAGORA_TELEMETRY_PROFILE", "quick")))`). |
| in | `out_summary` | String | n/a | no | Field `out_summary` (default `abspath(get(ENV, "SPACEAGORA_TELEMETRY_OUT_SUMMARY", joinpath(DEFAULT_OUTPUT_DIR, "telemetry_orbit_accuracy_summary.csv")))`). |
| in | `out_errors` | String | n/a | no | Field `out_errors` (default `abspath(get(ENV, "SPACEAGORA_TELEMETRY_OUT_ERRORS", joinpath(DEFAULT_OUTPUT_DIR, "telemetry_orbit_accuracy_errors.csv")))`). |
| in | `manifest_path` | String | n/a | no | Field `manifest_path` (default `abspath(get(ENV, "SPACEAGORA_TELEMETRY_MANIFEST", DEFAULT_MANIFEST_PATH))`). |
| in | `enforce` | Bool | n/a | no | Field `enforce` (default `false`). |
| in | `generate_plots` | Bool | n/a | no | Field `generate_plots` (default `_safe_parse_bool(get(ENV, "SPACEAGORA_TELEMETRY_PLOTS", "1"), true)`). |
| in | `scenarios` | Vector{String} | n/a | no | Field `scenarios` (default `String[]`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | VerificationRequest | n/a | — | Constructed `VerificationRequest` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.ic_fit__ic_fit_run|_ic_fit_run]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/ic_fit.jl:50-50`
- [[analysis.manifest_parsing__request_from_study_config|_request_from_study_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:773-773`
- [[analysis.runner_run_verification_cli|run_verification_cli]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:491-491`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/types.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__safe_parse_bool|_safe_parse_bool]] · `callers` · call · `src/analysis/verification/telemetry_verification/types.jl:262-262`
<!-- vulcan:connections:end -->

## Limitations
Because defaults read `ENV` at construction, two requests built at different times in the same session can differ silently if the environment changes. `profile` is not validated against the known set (`:quick`, `:full`). `abspath` is resolved relative to the current working directory at construction, not at execution. Output directories are not created here. The struct is immutable, so overriding one field means constructing a new request.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/types.jl` line 256.
