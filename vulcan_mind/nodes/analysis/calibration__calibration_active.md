---
id: analysis.calibration__calibration_active
label: _calibration_active
kind: function
source:
  file: src/analysis/verification/telemetry_verification/calibration.jl
  symbol: _calibration_active
  lines:
  - 1
  - 1
inputs:
- id: cal
  type: CalibrationConfig
  units: n/a
  required: true
  description: Positional argument `cal`.
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
  type: Bool
  units: n/a
  description: Return value of `_calibration_active`.
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

# _calibration_active

## Purpose

`_calibration_active(cal::CalibrationConfig, profile::Symbol)` is the single predicate deciding whether drag and SRP calibration should run for a given solver profile. It returns `true` only when the configuration's master switch is on and the requested profile is one the configuration opted in.

## Design & Implementation

The body is one `@inline`d boolean expression, `cal.enabled && (profile in cal.profiles)`, returning `Bool`. Short-circuit evaluation means the membership test on `cal.profiles` is skipped entirely when `cal.enabled` is `false`, so a disabled configuration costs one field read. Because the profile is a `Symbol`, the containment check is by identity against the profiles collection carried on the config, keeping the predicate allocation-free and safe to call inside the verification runner's per-profile loop.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cal` | CalibrationConfig | n/a | yes | Positional argument `cal`. |
| in | `profile` | Symbol | n/a | yes | Positional argument `profile`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_calibration_active`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_single_scenario|_run_single_scenario]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:137-137`
- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:376-376`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/calibration.jl`

**Downstream**

- `callees` → [[analysis.calibration__single_point_calibration|_single_point_calibration]] · `callers` · call · `src/analysis/verification/telemetry_verification/calibration.jl:6-6`
<!-- vulcan:connections:end -->

## Limitations

The check is purely declarative: it does not verify that the candidate grids are non-empty, that telemetry is available, or that the objective string is supported, so a profile can pass this gate and still fail later in the calibration path. It reads `cal` without copying, so mutating the config concurrently during a run gives an inconsistent answer.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/calibration.jl` line 1.
