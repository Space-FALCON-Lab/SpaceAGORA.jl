---
id: misc.main_run_cli
label: run_cli
kind: function
source:
  file: src/cli/main.jl
  symbol: run_cli
  lines:
  - 4
  - 4
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: SpaceAGORACLI command dispatcher reached through the SpaceAGORA package
    namespace as SpaceAGORA.run_cli.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: exit_code
  type: Int
  units: n/a
  description: Process exit status passed to exit(); zero on success, non-zero when
    the dispatched subcommand reports failure.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- cli
- entrypoint
charts:
- misc
origin: agent
---

# run_cli

## Purpose
`src/cli/main.jl` is the executable shim that turns the SpaceAGORA package into a command-line program. It contains no simulation logic of its own: it locates the repository root, loads the package source, hands the raw process arguments to `run_cli`, and terminates the process with whatever status that dispatcher returns. Keeping the shim this thin means the exit-code contract and the argument contract both live in testable library code rather than in a script that can only be exercised by spawning a subprocess.

## Model & Assumptions
The script assumes it is executed from its checked-out location, because `REPO_ROOT` is derived structurally as `normpath(joinpath(@__DIR__, "..", ".."))` — two directories above `src/cli`. It further assumes the package entry file lives at `src/SpaceAGORA.jl` relative to that root, and it `include`s that file directly rather than relying on the Julia package manager, so the CLI runs against the working tree even when the package is not installed into an environment.

## Design & Implementation
The four lines execute in order: define `REPO_ROOT`, include the package, then `exit(SpaceAGORA.run_cli(copy(ARGS)))`. The defensive `copy(ARGS)` matters — `ARGS` is a global mutable vector, and the dispatcher consumes arguments as it walks subcommands, so copying prevents the parser from mutating process-global state that other code or a later interactive session might read. Because `exit` receives the dispatcher's return value directly, `run_cli` must return an integer status; shells and CI jobs branch on that value, so a subcommand that fails must report it as a non-zero return rather than by throwing past the shim.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | SpaceAGORACLI command dispatcher reached through the SpaceAGORA package namespace as SpaceAGORA.run_cli. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `exit_code` | Int | n/a | — | Process exit status passed to exit(); zero on success, non-zero when the dispatched subcommand reports failure. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.cli|SpaceAGORACLI]] · `api` → `module_api` · call · `src/cli/main.jl`
- [[spaceagora.spaceagora_render_asset_report|render_asset_report]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:562-562`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Including the package on every invocation pays full load and compilation latency unless a precompiled system image is in play, which is exactly the cost the precompile workload exists to reduce. Relocating `main.jl` within the tree silently breaks the computed root. An exception escaping `run_cli` bypasses the exit-code path and surfaces as a Julia stack trace with the interpreter's own failure status.

## Provenance
Mapped from `src/cli/main.jl:1-4`.
