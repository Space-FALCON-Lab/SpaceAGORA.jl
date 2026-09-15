---
id: cli.spaceagora_cli__run_subprocess
label: _run_subprocess
kind: function
source:
  file: src/cli/spaceagora_cli.jl
  symbol: _run_subprocess
  lines:
  - 52
  - 52
inputs:
- id: script
  type: String
  units: n/a
  required: true
  description: Positional argument `script`.
- id: script_args
  type: Vector{String}
  units: n/a
  required: true
  description: Positional argument `script_args`.
- id: env_pairs
  type: Vector{Pair{String,String}}
  units: n/a
  required: false
  description: Keyword argument `env_pairs` (default `Pair{String,String}[]`).
- id: print_only
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `print_only` (default `false`).
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
  description: Return value of `_run_subprocess`.
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

# _run_subprocess

## Purpose
Common launcher used by every CLI subcommand that runs a script: builds a `julia --project=.AGORA <script> <args>` command, optionally prints it instead of executing, and otherwise runs it with extra environment variables and returns the child's exit code.

## Design & Implementation
Signature `_run_subprocess(script::String, script_args::Vector{String}; env_pairs::Vector{Pair{String,String}}=[], print_only::Bool=false, io::IO=stdout, errio::IO=stderr)::Int`. `Base.julia_cmd()` supplies the current Julia executable and flags; the command is interpolated as `` `$cmd --project=$DOT_AGORA_PROJECT $script $script_args` ``. In `print_only` mode it prints `project=`, `script=`, an `env:` block of `k=v` lines and `cmd=` to `io`, returning 0. Otherwise it creates `<repo>/output` with `mkpath`, applies `addenv(full, env_pairs...)` when pairs are present, and executes `run(pipeline(ignorestatus(cmd_env); stdout=io, stderr=errio), wait=true)`, returning `process.exitcode`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `script` | String | n/a | yes | Positional argument `script`. |
| in | `script_args` | Vector{String} | n/a | yes | Positional argument `script_args`. |
| in | `env_pairs` | Vector{Pair{String,String}} | n/a | no | Keyword argument `env_pairs` (default `Pair{String,String}[]`). |
| in | `print_only` | Bool | n/a | no | Keyword argument `print_only` (default `false`). |
| in | `io` | IO | n/a | no | Keyword argument `io` (default `stdout`). |
| in | `errio` | IO | n/a | no | Keyword argument `errio` (default `stderr`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Int | n/a | — | Return value of `_run_subprocess`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.spaceagora_cli__run_benchmark|_run_benchmark]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:181-181`
- [[cli.spaceagora_cli__run_example|_run_example]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:99-99`
- [[cli.spaceagora_cli__run_telemetry|_run_telemetry]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:139-139`
- [[module.cli|SpaceAGORACLI]] · `api` → `module_api` · call · `src/cli/spaceagora_cli.jl`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/cli/spaceagora_cli.jl:56-56`
<!-- vulcan:connections:end -->

## Limitations
`ignorestatus` means a failing child never raises; callers must inspect the returned code. The child inherits the parent's full environment plus overrides, so stale `SPACEAGORA_*` variables in the shell leak through. `io`/`errio` must be real streams accepted by `pipeline`; an `IOBuffer` works but a closed stream errors. Signals and the child's stdin are not managed.

## Provenance
Mapped from `src/cli/spaceagora_cli.jl` line 52.
