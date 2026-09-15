---
id: module.cli
label: SpaceAGORACLI
kind: module
source:
  file: src/cli/spaceagora_cli.jl
  symbol: SpaceAGORACLI
inputs: []
outputs:
- id: api
  type: Module
  units: n/a
  description: AssetCheckItem, AssetCheckReport, check_assets, render_asset_report,
    run_cli, and command dispatch helpers exported by SpaceAGORACLI.
tags:
- module
charts:
- master
origin: agent
---

# SpaceAGORACLI

## Purpose
`SpaceAGORACLI` is the thin command-line boundary for the package. It translates argument vectors into example, telemetry, benchmark, and asset-check actions, prints usage information, and returns integer process statuses. The module keeps subprocess execution and path resolution outside the simulation engine, allowing library callers to invoke the same Julia package without inheriting CLI output or process management.

## Theory & Math
The CLI has no physical model. Its observable contract is a finite dispatch function from an argument vector to an integer status code. Asset reporting adds a simple classification of expected files into present and missing sets, represented by `AssetCheckItem` and `AssetCheckReport` records.

## Model & Assumptions
Examples and benchmark scripts are resolved relative to the repository tree. Command names and option positions are interpreted exactly as implemented by `run_cli`; an unknown command follows the usage/error path. Subprocesses inherit selected environment pairs, and `print_only=true` suppresses execution while retaining the rendered command for diagnostics.

## Design & Implementation
`_resolve_example_path` maps a short example name to a checked-in Julia script. `_run_subprocess` builds the command, optionally prints environment assignments, and returns the child exit code. `_run_example`, `_run_telemetry`, and `_run_benchmark` validate their argument slices before delegating. `run_cli` is the public dispatcher and accepts explicit output streams, which makes the command path testable without redirecting global stdout.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `api` | Module | n/a | — | AssetCheckItem, AssetCheckReport, check_assets, render_asset_report, run_cli, and command dispatch helpers exported by SpaceAGORACLI. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `api` → [[cli.assets_assetcheckitem|AssetCheckItem]] · `module_api` · call · `src/cli/assets.jl`
- `api` → [[cli.assets_assetmanifestentry|AssetManifestEntry]] · `module_api` · call · `src/cli/assets.jl`
- `api` → [[cli.spaceagora_cli__normalize_example_name|_normalize_example_name]] · `module_api` · call · `src/cli/spaceagora_cli.jl`
- `api` → [[cli.spaceagora_cli__resolve_example_path|_resolve_example_path]] · `module_api` · call · `src/cli/spaceagora_cli.jl`
- `api` → [[cli.spaceagora_cli__run_subprocess|_run_subprocess]] · `module_api` · call · `src/cli/spaceagora_cli.jl`
- `api` → [[cli.spaceagora_cli__starts_with|_starts_with]] · `module_api` · call · `src/cli/spaceagora_cli.jl`
- `api` → [[cli.spaceagora_cli__value_after_equals|_value_after_equals]] · `module_api` · call · `src/cli/spaceagora_cli.jl`
- `api` → [[cli.spaceagora_cli_spaceagoracli|SpaceAGORACLI]] · `module_api` · call · `src/cli/spaceagora_cli.jl`
- `api` → [[misc.main_run_cli|run_cli]] · `module_api` · call · `src/cli/main.jl`
- `api` → [[module.spaceagora|SpaceAGORA]] · `cli` · call · `src/SpaceAGORA.jl:15-15`
<!-- vulcan:connections:end -->

## Limitations
The CLI relies on the local Julia executable and repository files. A missing example or unavailable benchmark dependency produces a nonzero status rather than a library exception. Subprocess output is intentionally textual, so consumers that need structured telemetry should use `TelemetryVerification` directly. Environment pairs are caller-controlled and therefore can change native-library behavior.

## Provenance
Mapped from `src/cli/spaceagora_cli.jl`.
