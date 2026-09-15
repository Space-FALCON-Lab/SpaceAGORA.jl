---
id: analysis.reporting__default_plots_outdir
label: _default_plots_outdir
kind: function
source:
  file: src/analysis/verification/telemetry_verification/reporting.jl
  symbol: _default_plots_outdir
  lines:
  - 3
  - 3
inputs:
- id: out_summary
  type: String
  units: n/a
  required: true
  description: Positional argument `out_summary`.
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
  description: Return value of `_default_plots_outdir`.
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

# _default_plots_outdir

## Purpose
Derives the directory into which telemetry accuracy plots are written, placing it beside the summary CSV and suffixing the profile name so quick and full runs never overwrite each other's figures.

## Design & Implementation
Takes `out_summary::String` (path of the summary CSV) and `profile::Symbol`, and returns `normpath(joinpath(dirname(out_summary), "telemetry_plots_<profile>"))`. `String(profile)` converts the symbol without a leading colon. `normpath` collapses `..` and duplicate separators but does not make the path absolute, so a relative `out_summary` produces a relative output directory. The directory is not created here; `_generate_plots` passes it to the plotting script which is responsible for `mkpath`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `out_summary` | String | n/a | yes | Positional argument `out_summary`. |
| in | `profile` | Symbol | n/a | yes | Positional argument `profile`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_default_plots_outdir`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.reporting__generate_plots|_generate_plots]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/reporting.jl:10-10`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
If `out_summary` has no directory component, `dirname` returns an empty string and the plots land in `telemetry_plots_<profile>` under the current working directory, which may differ from where the CSV was written. Symbols containing path separators are not sanitised.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/reporting.jl` line 3.
