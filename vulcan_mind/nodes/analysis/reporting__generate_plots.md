---
id: analysis.reporting__generate_plots
label: _generate_plots
kind: function
source:
  file: src/analysis/verification/telemetry_verification/reporting.jl
  symbol: _generate_plots
  lines:
  - 7
  - 7
inputs:
- id: summary_csv
  type: String
  units: n/a
  required: true
  description: Positional argument `summary_csv`.
- id: errors_csv
  type: String
  units: n/a
  required: true
  description: Positional argument `errors_csv`.
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
  type: String
  units: n/a
  description: Return value of `_generate_plots`.
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

# _generate_plots

## Purpose
Spawns a child Julia process that runs `scripts/plotting/telemetry_orbit_accuracy_plots.jl` against the summary and error CSVs produced by a verification run, returning the directory the figures were written to.

## Design & Implementation
Builds `plot_script` from `REPO_ROOT`, throwing `ErrorException("Missing plotting script: ...")` if it does not exist. The output directory comes from `_default_plots_outdir(summary_csv, profile)`. The project used for the child process is `REPO_ROOT/.AGORA` when that directory exists, otherwise `REPO_ROOT` itself. The command is `Base.julia_cmd()` with `--startup-file=no --project=<plot_project>` followed by the script and `--summary=`, `--errors=`, `--outdir=` arguments, executed with `run(cmd)`, which blocks until the child exits and throws `ProcessFailedException` on a non-zero status.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `summary_csv` | String | n/a | yes | Positional argument `summary_csv`. |
| in | `errors_csv` | String | n/a | yes | Positional argument `errors_csv`. |
| in | `profile` | Symbol | n/a | yes | Positional argument `profile`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_generate_plots`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:446-446`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`

**Downstream**

- `callees` → [[analysis.reporting__default_plots_outdir|_default_plots_outdir]] · `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:10-10`
<!-- vulcan:connections:end -->

## Limitations
Each call pays a full Julia startup plus package load in the child, so this is slow relative to the verification itself. Output of the child goes to the parent's stdio with no capture. Paths containing spaces are safe because `Cmd` interpolation quotes them, but the script itself must parse `--key=value` flags exactly as passed. There is no timeout on the child process.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/reporting.jl` line 7.
