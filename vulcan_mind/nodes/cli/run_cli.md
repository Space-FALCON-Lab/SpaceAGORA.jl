---
id: cli.run_cli
label: run_cli
kind: function
source:
  file: src/cli/spaceagora_cli.jl
  symbol: run_cli
  lines:
  - 184
  - 214
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: SpaceAGORACLI namespace containing command dispatch and reporting helpers.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: status
  type: Int
  units: n/a
  description: Process-style status code returned after command dispatch and child
    execution.
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

# run_cli

## Purpose
`run_cli` dispatches SpaceAGORA command-line arguments to the example, telemetry, benchmark, and asset-check handlers. It is the public adapter used by `src/cli/main.jl`, and its explicit IO arguments let tests and embedding callers capture messages without replacing global streams.

## Theory & Math
The function is a finite dispatcher from `Vector{String}` to an integer status. Valid commands map to handler status codes; malformed or unknown commands map to the usage/error branch. No spacecraft state or numerical integration occurs in this function itself.

## Model & Assumptions
Arguments follow the command grammar implemented by the four handlers. Example and benchmark names resolve inside the repository, and subprocess calls depend on the local Julia executable and environment. The caller owns the meaning of nonzero status values and may choose whether to terminate the process.

## Design & Implementation
`run_cli` inspects the first argument, routes to `_run_example`, `_run_telemetry`, `_run_benchmark`, or the asset-report path, and writes diagnostics to `io` or `errio`. `_print_usage` is used for missing and unknown commands. The wrapper `SpaceAGORA.run_cli` forwards keyword streams to this implementation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | SpaceAGORACLI namespace containing command dispatch and reporting helpers. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `status` | Int | n/a | — | Process-style status code returned after command dispatch and child execution. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[spaceagora.spaceagora_render_asset_report|render_asset_report]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:562-562`

**Downstream**

- `callees` → [[cli.assets_load_asset_manifest|load_asset_manifest]] · `callers` · call · `src/cli/spaceagora_cli.jl:197-197`
- `callees` → [[cli.assets_render_asset_manifest|render_asset_manifest]] · `callers` · call · `src/cli/spaceagora_cli.jl:197-197`
- `callees` → [[cli.assets_render_asset_report|render_asset_report]] · `callers` · call · `src/cli/spaceagora_cli.jl:194-194`
- `callees` → [[cli.assets_setup_open_assets|setup_open_assets]] · `callers` · call · `src/cli/spaceagora_cli.jl:200-200`
- `callees` → [[cli.spaceagora_cli__print_usage|_print_usage]] · `callers` · call · `src/cli/spaceagora_cli.jl:185-185`
- `callees` → [[cli.spaceagora_cli__run_benchmark|_run_benchmark]] · `callers` · call · `src/cli/spaceagora_cli.jl:209-209`
- `callees` → [[cli.spaceagora_cli__run_example|_run_example]] · `callers` · call · `src/cli/spaceagora_cli.jl:205-205`
- `callees` → [[cli.spaceagora_cli__run_telemetry|_run_telemetry]] · `callers` · call · `src/cli/spaceagora_cli.jl:207-207`
- `callees` → [[misc.assets_check_assets|check_assets]] · `callers` · feedback · `src/cli/spaceagora_cli.jl:194-194`
- `callees` → [[spaceagora.spaceagora_check_assets|check_assets]] · `callers` · feedback · `src/cli/spaceagora_cli.jl:194-194`
- `callees` → [[spaceagora.spaceagora_render_asset_report|render_asset_report]] · `callers` · feedback · `src/cli/spaceagora_cli.jl:194-194`
<!-- vulcan:connections:end -->

## Limitations
The dispatcher does not validate external solver installations or native-library compatibility before launching a child. A child can produce partial files before returning failure. Because command output is human-readable, integrations needing structured data should call the underlying analysis or asset APIs rather than parse text.

## Provenance
Mapped from `src/cli/spaceagora_cli.jl:184-214`.
