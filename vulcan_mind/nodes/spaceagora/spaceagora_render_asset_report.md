---
id: spaceagora.spaceagora_render_asset_report
label: render_asset_report
kind: function
source:
  file: src/SpaceAGORA.jl
  symbol: render_asset_report
  lines:
  - 559
  - 559
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
  type: SpaceAGORACLI.render_asset_report
  units: n/a
  description: Return value of `render_asset_report`. Returns `SpaceAGORACLI.render_asset_report(args...;
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

# render_asset_report

## Purpose
`render_asset_report` prints an `AssetCheckReport` as a human-readable text block, one entry per asset root, so that `spaceagora assets check` and interactive users get the same formatted status listing. The package-level symbol forwards to `SpaceAGORACLI.render_asset_report`.

## Design & Implementation
Signature `render_asset_report(report::AssetCheckReport; io::IO=stdout)`. It writes a header `SpaceAGORA asset check` and `repo_root=<path>`, then for each `AssetCheckItem` computes a status word: `available` when `item.available`, otherwise `missing-required` if `item.required` else `missing-optional`. Four lines follow per item (`- name: status`, then indented `scope`, `path`, and `detail`). The function returns `nothing` and performs no filesystem access; all state comes from the report built by `check_assets`. Output is plain `println` to the supplied stream so tests can capture it with an `IOBuffer`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Vararg{Any} | n/a | no | Positional argument `args` (variadic). |
| in | `kwargs` | Vararg{Any} | n/a | no | Keyword argument `kwargs` (variadic). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SpaceAGORACLI.render_asset_report | n/a | — | Return value of `render_asset_report`. Returns `SpaceAGORACLI.render_asset_report(args...; kwargs...)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.assets_setup_open_assets|setup_open_assets]] · `callees` → `callers` · feedback · `src/cli/assets.jl:154-154`
- [[cli.run_cli|run_cli]] · `callees` → `callers` · feedback · `src/cli/spaceagora_cli.jl:194-194`
- [[spaceagora.spaceagora_check_assets|check_assets]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:555-555`

**Downstream**

- `callees` → [[cli.run_cli|run_cli]] · `callers` · call · `src/SpaceAGORA.jl:562-562`
- `callees` → [[misc.main_run_cli|run_cli]] · `callers` · call · `src/SpaceAGORA.jl:562-562`
- `callees` → [[spaceagora.spaceagora_run_cli|run_cli]] · `callers` · call · `src/SpaceAGORA.jl:562-562`
<!-- vulcan:connections:end -->

## Limitations
The format is fixed text with no machine-readable option (no JSON or TOML rendering), so consumers must parse indentation and prefix strings. Nothing is done with the status beyond printing: the function does not set a return code or throw when required assets are missing, leaving that decision to the caller. Very long `detail` strings are printed unwrapped. The forwarding method at package level accepts `args...` so a wrong argument type fails only inside `SpaceAGORACLI`.

## Provenance
Mapped from `src/SpaceAGORA.jl` line 559.
