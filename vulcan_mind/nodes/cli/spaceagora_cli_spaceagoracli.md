---
id: cli.spaceagora_cli_spaceagoracli
label: SpaceAGORACLI
kind: module
source:
  file: src/cli/spaceagora_cli.jl
  symbol: SpaceAGORACLI
  lines:
  - 1
  - 1
inputs:
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
  type: Any
  units: n/a
  description: Value produced by this symbol.
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

# SpaceAGORACLI

## Purpose
Module implementing the `spaceagora` command-line front end: it resolves repository paths, includes the asset-management helpers, and dispatches `run`, `telemetry`, `benchmark`, `assets` and `help` commands to child Julia processes or in-process asset routines.

## Design & Implementation
`SpaceAGORACLI` exports `AssetCheckItem`, `AssetCheckReport`, `check_assets`, `render_asset_report` and `run_cli`. Constants fix `REPO_ROOT` (two levels above the file), `DOT_AGORA_PROJECT = <repo>/.AGORA`, `EXAMPLES_DIR`, and the three study launchers under `benchmarks/studies/`. `assets.jl` is included for manifest and asset checks. `run_cli(args=copy(ARGS); io, errio)` returns an `Int` exit code: no args or `help`/`--help`/`-h` print usage; `assets` requires `check`, `manifest` or `setup-open`; `run`, `telemetry` and `benchmark` delegate to `_run_example`, `_run_telemetry` and `_run_benchmark`; unknown commands throw `ArgumentError`. Sub-parsers share the `_starts_with`/`_value_after_equals` idiom and `_run_subprocess` for execution.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.cli|SpaceAGORACLI]] · `api` → `module_api` · call · `src/cli/spaceagora_cli.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Errors are signalled by throwing `ArgumentError` from `run_cli` rather than returning a non-zero code, so a wrapper script must catch them to produce a clean exit status. Paths are derived from `@__DIR__` at load time, so a relocated or precompiled copy of the module still points at the original checkout. All heavy work runs in a separate Julia process, incurring a fresh startup and compile cost per command.

## Provenance
Mapped from `src/cli/spaceagora_cli.jl` line 1.
