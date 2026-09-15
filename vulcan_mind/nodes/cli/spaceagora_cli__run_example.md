---
id: cli.spaceagora_cli__run_example
label: _run_example
kind: function
source:
  file: src/cli/spaceagora_cli.jl
  symbol: _run_example
  lines:
  - 73
  - 73
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
  description: Return value of `_run_example`.
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

# _run_example

## Purpose
Implements `spaceagora run --example=<file>`: parses the run flags, resolves the example script, and launches it in a fresh Julia process under the `.AGORA` project with optional smoke-mode and output-directory environment variables.

## Design & Implementation
Iterates over `args::Vector{String}`; `--example=` sets the script name, `--output-dir=` sets `output_dir` (made absolute), `--smoke` and `--print-only` set booleans, and any other token throws `ArgumentError("Unknown run argument ...")`. A missing example throws `ArgumentError("run requires --example=<file>.")`. The script is resolved via `_resolve_example_path`. Environment pairs are built as `SPACEAGORA_CLI_OUTPUT_DIR => output_dir` when given and, for smoke, `SPACEAGORA_EXAMPLE_SMOKE => "1"` plus `SPACEAGORA_EXAMPLE_SMOKE_RESULTS => "1"`. It returns the integer exit code from `_run_subprocess(script, String[]; env_pairs, print_only, io, errio)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Vector{String} | n/a | yes | Positional argument `args`. |
| in | `io` | IO | n/a | no | Keyword argument `io` (default `stdout`). |
| in | `errio` | IO | n/a | no | Keyword argument `errio` (default `stderr`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_run_example`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.run_cli|run_cli]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:205-205`

**Downstream**

- `callees` → [[cli.spaceagora_cli__resolve_example_path|_resolve_example_path]] · `callers` · call · `src/cli/spaceagora_cli.jl:92-92`
- `callees` → [[cli.spaceagora_cli__run_subprocess|_run_subprocess]] · `callers` · call · `src/cli/spaceagora_cli.jl:99-99`
- `callees` → [[cli.spaceagora_cli__starts_with|_starts_with]] · `callers` · call · `src/cli/spaceagora_cli.jl:79-79`
- `callees` → [[cli.spaceagora_cli__value_after_equals|_value_after_equals]] · `callers` · call · `src/cli/spaceagora_cli.jl:80-80`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/cli/spaceagora_cli.jl:94-94`
<!-- vulcan:connections:end -->

## Limitations
Flags must use the `--flag=value` form; `--example foo` (space separated) is rejected as two unknown tokens. No arguments can be passed through to the example script itself since `script_args` is always empty. Repeated flags silently overwrite earlier values.

## Provenance
Mapped from `src/cli/spaceagora_cli.jl` line 73.
