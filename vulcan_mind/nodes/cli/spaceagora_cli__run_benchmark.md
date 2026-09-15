---
id: cli.spaceagora_cli__run_benchmark
label: _run_benchmark
kind: function
source:
  file: src/cli/spaceagora_cli.jl
  symbol: _run_benchmark
  lines:
  - 142
  - 142
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
  description: Return value of `_run_benchmark`.
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

# _run_benchmark

## Purpose
Implements `spaceagora benchmark <mode> ...`: selects one of the two benchmark study launchers, parses the output-directory and print-only flags, forwards remaining tokens (such as `quick`/`full`/`smoke`) to the study script, and runs it in a child Julia process.

## Design & Implementation
Takes `args::Vector{String}` with keyword `io` and `errio`. An empty vector throws `ArgumentError`. `mode = first(args)` must be `runtime-analysis` (launcher `PERF_RUNTIME_LAUNCHER`, default output `<repo>/output/performance`, env `SPACEAGORA_PERF_OUTDIR`) or `smart-parallel-ladder` (launcher `SMART_LADDER_LAUNCHER`, default `<repo>/output/performance/smart_parallel_ladder`, env `SPACEAGORA_SMART_LADDER_OUTDIR`); any other mode throws. For each remaining token, `--output-dir=` overrides the directory (made absolute), `--print-only` sets the dry-run flag, and everything else is appended verbatim to `script_args`. It calls `mkpath(output_dir)` and returns the exit code of `_run_subprocess(launcher, script_args; env_pairs, print_only, io, errio)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Vector{String} | n/a | yes | Positional argument `args`. |
| in | `io` | IO | n/a | no | Keyword argument `io` (default `stdout`). |
| in | `errio` | IO | n/a | no | Keyword argument `errio` (default `stderr`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_run_benchmark`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.run_cli|run_cli]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:209-209`

**Downstream**

- `callees` → [[cli.spaceagora_cli__run_subprocess|_run_subprocess]] · `callers` · call · `src/cli/spaceagora_cli.jl:181-181`
- `callees` → [[cli.spaceagora_cli__starts_with|_starts_with]] · `callers` · call · `src/cli/spaceagora_cli.jl:155-155`
- `callees` → [[cli.spaceagora_cli__value_after_equals|_value_after_equals]] · `callers` · call · `src/cli/spaceagora_cli.jl:156-156`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/cli/spaceagora_cli.jl:160-160`
<!-- vulcan:connections:end -->

## Limitations
Unrecognised flags are silently forwarded to the study script rather than rejected, unlike `_run_example` and `_run_telemetry`. The output directory is created even in `--print-only` mode. The two mode branches duplicate the parsing loop verbatim. Profile tokens are not validated here, so a typo surfaces only from the launched script.

## Provenance
Mapped from `src/cli/spaceagora_cli.jl` line 142.
