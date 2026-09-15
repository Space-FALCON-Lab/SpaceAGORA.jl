---
id: spaceagora.spaceagora_run_cli
label: run_cli
kind: function
source:
  file: src/SpaceAGORA.jl
  symbol: run_cli
  lines:
  - 573
  - 573
inputs:
- id: args
  type: Vararg{Any}
  units: n/a
  required: false
  description: Positional argument `args` (variadic).
- id: kwargs
  type: Vararg{Any}
  units: n/a
  required: false
  description: Keyword argument `kwargs` (variadic).
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
  type: SpaceAGORACLI.run_cli
  units: n/a
  description: Return value of `run_cli`. Returns `SpaceAGORACLI.run_cli(args...;
    kwargs...)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- spaceagora
charts:
- spaceagora
origin: agent
---

# run_cli

## Purpose
`run_cli` is the stable, package-owned command-line entrypoint used by the `bin/spaceagora` wrapper. From the top-level module it simply forwards to `SpaceAGORACLI.run_cli`, which dispatches the `run`, `telemetry`, `benchmark`, `help` and `assets {check|manifest|setup-open}` commands and returns a process exit code.

## Design & Implementation
Signature `run_cli(args::Vector{String}=copy(ARGS); io::IO=stdout, errio::IO=stderr)::Int`. The wrapper in `SpaceAGORA.jl` is `run_cli(args...; kwargs...) = SpaceAGORACLI.run_cli(args...; kwargs...)`. Inside `SpaceAGORACLI`, an empty vector or `help`/`--help`/`-h` prints usage and returns its code; `assets` requires a subcommand and routes `check` to `render_asset_report(check_assets(); io)`, `manifest` to `render_asset_manifest(load_asset_manifest())`, and `setup-open` to `setup_open_assets`; `run`, `telemetry` and `benchmark` delegate to `_run_example`, `_run_telemetry` and `_run_benchmark` with the remaining `tail` arguments. All output goes to the injected `io`/`errio` streams so tests can capture it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Vararg{Any} | n/a | no | Positional argument `args` (variadic). |
| in | `kwargs` | Vararg{Any} | n/a | no | Keyword argument `kwargs` (variadic). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SpaceAGORACLI.run_cli | n/a | — | Return value of `run_cli`. Returns `SpaceAGORACLI.run_cli(args...; kwargs...)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.spaceagora|SpaceAGORA]] · `api` → `module_api` · call · `src/SpaceAGORA.jl`
- [[spaceagora.spaceagora_render_asset_report|render_asset_report]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:562-562`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Unknown commands and missing or unknown `assets` subcommands are reported by throwing `ArgumentError` rather than by returning a non-zero exit code, so a wrapper must catch exceptions to produce clean shell exits. `copy(ARGS)` is evaluated at call time, meaning a default call inside a running REPL picks up the REPL's own argument vector. Global flags are not parsed at this level; every command parses its own `tail` independently and the forwarding method offers no argument typing.

## Provenance
Mapped from `src/SpaceAGORA.jl` line 573.
