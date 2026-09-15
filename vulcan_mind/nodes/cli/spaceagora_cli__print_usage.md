---
id: cli.spaceagora_cli__print_usage
label: _print_usage
kind: function
source:
  file: src/cli/spaceagora_cli.jl
  symbol: _print_usage
  lines:
  - 40
  - 40
inputs:
- id: io
  type: IO
  units: n/a
  required: false
  description: Positional argument `io` (default `stdout`).
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
  description: Return value of `_print_usage`. Returns `0`.
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

# _print_usage

## Purpose
Writes the fixed command synopsis for the `spaceagora` CLI to an IO stream and returns exit code 0, serving both the no-argument and `help`/`--help`/`-h` cases in `run_cli`.

## Design & Implementation
`_print_usage(io::IO=stdout)` issues nine `println` calls listing the `run`, `telemetry`, `benchmark runtime-analysis`, `benchmark smart-parallel-ladder` and three `assets` subcommands with their accepted flags (`--example=`, `--output-dir=`, `--smoke`, `--print-only`, `--enforce=0|1`, `--plots=0|1`, profile tokens `quick|full|smoke`). The integer 0 is returned so `run_cli` can propagate it as a process exit status.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `io` | IO | n/a | no | Positional argument `io` (default `stdout`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_print_usage`. Returns `0`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.run_cli|run_cli]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:185-185`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/cli/spaceagora_cli.jl:41-41`
<!-- vulcan:connections:end -->

## Limitations
The usage text is a hand-maintained string literal, so flags accepted by the parsers (for example `--profile=` and `--scenarios=` in `_run_telemetry`) are not all listed, and additions to the parsers do not update it. No version or description line is printed.

## Provenance
Mapped from `src/cli/spaceagora_cli.jl` line 40.
