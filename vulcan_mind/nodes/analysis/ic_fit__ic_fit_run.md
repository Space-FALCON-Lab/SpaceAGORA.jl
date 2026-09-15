---
id: analysis.ic_fit__ic_fit_run
label: _ic_fit_run
kind: function
source:
  file: src/analysis/verification/telemetry_verification/ic_fit.jl
  symbol: _ic_fit_run
  lines:
  - 31
  - 31
inputs:
- id: manifest
  type: Dict{String, Any}
  units: n/a
  required: true
  description: Positional argument `manifest`.
- id: scenario_index
  type: Int
  units: n/a
  required: true
  description: Positional argument `scenario_index`.
- id: offsets_m
  type: NTuple{3, Float64}
  units: n/a
  required: true
  description: Positional argument `offsets_m`.
- id: offsets_mps
  type: NTuple{3, Float64}
  units: n/a
  required: true
  description: Positional argument `offsets_mps`.
- id: label
  type: String
  units: n/a
  required: true
  description: Positional argument `label`.
- id: workdir
  type: String
  units: n/a
  required: true
  description: Positional argument `workdir`.
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
  description: Return value of `_ic_fit_run`. Returns `(`.
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

# _ic_fit_run

## Purpose

Executes one propagation of the fit campaign at a specified Cartesian initial-condition offset and returns both its per-sample residual series and its summary table. It is the single primitive the differential correction calls for its baseline, its six finite-difference perturbations, and its validation run.

## Design & Implementation

Takes `manifest::Dict{String, Any}`, the `scenario_index` to perturb, `offsets_m` and `offsets_mps` as `NTuple{3, Float64}`, a `label` used in filenames, a `workdir`, and a run `profile::Symbol`. It `deepcopy`s the manifest so the caller's dictionary is never mutated, writes `ic_offset_m` and `ic_offset_mps` into the selected scenario, and serialises the copy to `icfit_manifest_<label>.toml` via `TOML.print`. It then builds a `VerificationRequest` with `enforce=false`, `generate_plots=false`, and `scenarios=[scen["name"]]` — restricting the run to the single fitted scenario so the other manifest scenarios are neither propagated nor able to filter it out — runs `run_verification`, and returns `(series=_ic_fit_series(req.out_errors, name), summary=result.summary)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `manifest` | Dict{String, Any} | n/a | yes | Positional argument `manifest`. |
| in | `scenario_index` | Int | n/a | yes | Positional argument `scenario_index`. |
| in | `offsets_m` | NTuple{3, Float64} | n/a | yes | Positional argument `offsets_m`. |
| in | `offsets_mps` | NTuple{3, Float64} | n/a | yes | Positional argument `offsets_mps`. |
| in | `label` | String | n/a | yes | Positional argument `label`. |
| in | `workdir` | String | n/a | yes | Positional argument `workdir`. |
| in | `profile` | Symbol | n/a | yes | Positional argument `profile`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_ic_fit_run`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[envana.ana_ic_fit_fit_initial_state|fit_initial_state]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/ic_fit.jl:116-116`

**Downstream**

- `callees` → [[analysis.ic_fit__ic_fit_series|_ic_fit_series]] · `callers` · call · `src/analysis/verification/telemetry_verification/ic_fit.jl:61-61`
- `callees` → [[analysis.run_verification|run_verification]] · `callers` · call · `src/analysis/verification/telemetry_verification/ic_fit.jl:59-59`
- `callees` → [[analysis.types_verificationrequest|VerificationRequest]] · `callers` · call · `src/analysis/verification/telemetry_verification/ic_fit.jl:50-50`
<!-- vulcan:connections:end -->

## Limitations

Filenames are derived solely from `label`, so two concurrent fits sharing a `workdir`, or a repeated label, overwrite each other's manifest, summary and errors CSVs; the function takes no lock and does nothing to make the writes atomic. Nothing is cleaned up, so a campaign leaves eight manifests plus sixteen CSVs behind. Offsets are written as absolute values into the scenario, so the caller — not this function — is responsible for adding perturbations on top of any offsets already in the manifest. `deepcopy` of a large manifest is repeated on every one of the eight runs.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/ic_fit.jl` line 31.
