---
id: envana.ana_ic_fit_fit_initial_state
label: fit_initial_state
kind: function
source:
  file: src/analysis/verification/telemetry_verification/ic_fit.jl
  symbol: fit_initial_state
  lines:
  - 90
  - 194
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: TelemetryVerification namespace supplying manifest parsing and the
    perturbed scenario runner.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: ic_offsets
  type: NamedTuple
  units: m,m/s
  description: Fitted Cartesian position and velocity offsets applied to the manifest
    initial condition.
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
# fit_initial_state

## Purpose
`fit_initial_state` estimates the six-element Cartesian offset that best aligns a simulated trajectory with telemetry for one manifest scenario, correcting an initial condition that is close but not exact.

## Theory & Math
The method is a finite-difference Gauss-Newton step. Writing `d` for the six-vector of offsets, `r(d)` for the residual series between simulation and telemetry, and `r_0` for the residual at the manifest baseline, the sensitivity matrix is built column by column from `J[:, k] = (r(d_0 + h_k e_k) - r_0) / h_k`, where `e_k` is the k-th unit vector and `h_k` is the corresponding step. Position steps default to `pos_step_m = 100.0` metres and velocity steps to `vel_step_mps = 0.005` metres per second. The correction solves the linear least squares problem `min |J delta + r_0|_2`, and the reported offset is `d_0 + delta`.

## Model & Assumptions
The scenario must expose Cartesian initial-condition telemetry columns; the function checks for all six of `x_ic`, `y_ic`, `z_ic`, `vx_ic`, `vy_ic`, `vz_ic` and raises an `ArgumentError` naming the missing key otherwise. Baseline offsets are read from the manifest fields `ic_offset_m` and `ic_offset_mps`, defaulting to zero triples. Sensitivities are assumed locally linear over the chosen step sizes, which holds only when the trajectory response to a 100 metre position shift is small compared with the residual scale. Both step sizes must be strictly positive.

## Design & Implementation
Seven simulations are executed through `_ic_fit_run`: one baseline labelled `"base"` and six perturbations labelled `"p1"` through `"p6"`, each writing into `workdir`, which defaults to a fresh `mktempdir()`. The residual axes are intersected across the events in `_IC_FIT_EVENTS` using `intersect` over the key sets, so only sample points present in every run and every event contribute to the fit. The `validate` keyword controls whether a confirming run is executed after the correction is computed.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | TelemetryVerification namespace supplying manifest parsing and the perturbed scenario runner. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `ic_offsets` | NamedTuple | m,m/s | — | Fitted Cartesian position and velocity offsets applied to the manifest initial condition. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.ic_fit__ic_fit_mean_axis_rmse|_ic_fit_mean_axis_rmse]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/ic_fit.jl:73-73`

**Downstream**

- `callees` → [[analysis.ic_fit__ic_fit_mean_axis_rmse|_ic_fit_mean_axis_rmse]] · `callers` · feedback · `src/analysis/verification/telemetry_verification/ic_fit.jl:165-165`
- `callees` → [[analysis.ic_fit__ic_fit_run|_ic_fit_run]] · `callers` · call · `src/analysis/verification/telemetry_verification/ic_fit.jl:116-116`
- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/analysis/verification/telemetry_verification/ic_fit.jl:178-178`
<!-- vulcan:connections:end -->

## Limitations
Seven full simulations make each call expensive, and the cost scales linearly if the step is iterated. Only a single Gauss-Newton step is taken, so a strongly nonlinear response is not converged. The fixed step sizes are not scaled to the orbit regime, so the same 100 metre perturbation may be too small in a low orbit and too large in a high one. The scenario must be listed in the manifest under a matching `name` field or the lookup raises.

## Provenance
Read directly from `src/analysis/verification/telemetry_verification/ic_fit.jl:90-194`, including the telemetry column validation, the six-element perturbation table, and the shared-axis intersection.
