---
id: cli.spaceagora_cli__run_telemetry
label: _run_telemetry
kind: function
source:
  file: src/cli/spaceagora_cli.jl
  symbol: _run_telemetry
  lines:
  - 102
  - 102
inputs:
- id: args
  type: Vector{String}
  units: n/a
  required: true
  description: Positional argument `args`.
- id: io
  type: IO
  units: n/a
  required: false
  description: Keyword argument `io` (default `stdout`).
- id: errio
  type: IO
  units: n/a
  required: false
  description: Keyword argument `errio` (default `stderr`).
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
  type: Int
  units: n/a
  description: Return value of `_run_telemetry`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- cli
charts:
- cli
origin: agent
---

# _run_telemetry

## Purpose
Implements `spaceagora telemetry [quick|full|smoke] [flags]`: parses the profile and options, maps them onto the `SPACEAGORA_TELEMETRY_*` environment contract expected by the telemetry orbit-accuracy study, and launches that study script.

## Design & Implementation
Defaults are `profile="quick"`, `output_dir=<repo>/output/telemetry`, `enforce=false`, `generate_plots=false`, `scenarios=""`. Bare `quick`/`full` tokens set the profile; `smoke` maps to `quick`; `--profile=` is lowercased and stripped with the same smoke aliasing; `--output-dir=` is made absolute; `--enforce=` and `--plots=` are true for `1`, `true`, `yes`, `on`; `--scenarios=` is stripped; `--print-only` sets the dry-run flag; anything else throws `ArgumentError`. After `mkpath(output_dir)` it sets five env pairs (`SPACEAGORA_TELEMETRY_OUT_SUMMARY`, `_OUT_ERRORS` pointing at CSVs in the output directory, `_PLOTS`, `_ENFORCE` as "1"/"0", and `_SCENARIOS`) and returns `_run_subprocess(TELEMETRY_LAUNCHER, [profile]; ...)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Vector{String} | n/a | yes | Positional argument `args`. |
| in | `io` | IO | n/a | no | Keyword argument `io` (default `stdout`). |
| in | `errio` | IO | n/a | no | Keyword argument `errio` (default `stderr`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_run_telemetry`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.run_cli|run_cli]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:207-207`

**Downstream**

- `callees` → [[cli.spaceagora_cli__run_subprocess|_run_subprocess]] · `callers` · call · `src/cli/spaceagora_cli.jl:139-139`
- `callees` → [[cli.spaceagora_cli__starts_with|_starts_with]] · `callers` · call · `src/cli/spaceagora_cli.jl:114-114`
- `callees` → [[cli.spaceagora_cli__value_after_equals|_value_after_equals]] · `callers` · call · `src/cli/spaceagora_cli.jl:115-115`
<!-- vulcan:connections:end -->

## Limitations
Only the profile is passed as a script argument; the remaining settings ride on environment variables, so the launcher and this parser must stay in sync by convention. The output directory is created even for `--print-only`. An arbitrary `--profile=` string is forwarded unvalidated to the launcher. `SPACEAGORA_TELEMETRY_SCENARIOS` is always set, even to an empty string.

## Provenance
Mapped from `src/cli/spaceagora_cli.jl` line 102.
