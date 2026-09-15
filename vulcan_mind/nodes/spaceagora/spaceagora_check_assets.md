---
id: spaceagora.spaceagora_check_assets
label: check_assets
kind: function
source:
  file: src/SpaceAGORA.jl
  symbol: check_assets
  lines:
  - 552
  - 552
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
  type: SpaceAGORACLI.check_assets
  units: n/a
  description: Return value of `check_assets`. Returns `SpaceAGORACLI.check_assets(args...;
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

# check_assets

## Purpose
`check_assets` inspects the repository's asset layout and returns an `AssetCheckReport` describing which baseline, optional, and high-fidelity asset roots (SPICE kernels, station geometry, atmosphere tables) are present on disk. The top-level module re-exports it so users and the `assets check` CLI command share one implementation.

## Design & Implementation
The package-level definition is a forwarder, `check_assets(args...; kwargs...) = SpaceAGORACLI.check_assets(args...; kwargs...)`. The implementation signature is `check_assets(; repo_root::String=REPO_ROOT, manifest_path::String=joinpath(repo_root, "data", "assets_manifest.toml"))::AssetCheckReport`. It calls `load_asset_manifest` to parse the TOML manifest into `AssetManifestEntry` records, then for each entry builds an `AssetCheckItem` with `name`, `scope`, the resolved absolute `path` from `_manifest_entry_path`, `available` from `_manifest_entry_available` (a filesystem existence test), the manifest's `required` flag, and a `detail` string that always embeds `licensing=...`. The items are wrapped with `repo_root` in the returned report; nothing is mutated and no I/O beyond reads occurs.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | Vararg{Any} | n/a | no | Positional argument `args` (variadic). |
| in | `kwargs` | Vararg{Any} | n/a | no | Keyword argument `kwargs` (variadic). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SpaceAGORACLI.check_assets | n/a | — | Return value of `check_assets`. Returns `SpaceAGORACLI.check_assets(args...; kwargs...)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.assets__manifest_entry_available|_manifest_entry_available]] · `callees` → `callers` · call · `src/cli/assets.jl:79-79`
- [[cli.assets_assetcheckitem|AssetCheckItem]] · `callees` → `callers` · call · `src/cli/assets.jl:20-20`
- [[cli.assets_setup_open_assets|setup_open_assets]] · `callees` → `callers` · feedback · `src/cli/assets.jl:136-136`
- [[cli.run_cli|run_cli]] · `callees` → `callers` · feedback · `src/cli/spaceagora_cli.jl:194-194`
- [[spaceagora.spaceagora_load_nbody_ephemeris_cache_bang|load_nbody_ephemeris_cache!]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:547-547`

**Downstream**

- `callees` → [[cli.assets_render_asset_report|render_asset_report]] · `callers` · call · `src/SpaceAGORA.jl:555-555`
- `callees` → [[spaceagora.spaceagora_render_asset_report|render_asset_report]] · `callers` · call · `src/SpaceAGORA.jl:555-555`
<!-- vulcan:connections:end -->

## Limitations
Availability is a pure existence check: a present but truncated or corrupt kernel file is reported as available. `REPO_ROOT` is captured at package load, so the default resolves relative to the installed package rather than the current working directory despite the docstring's `repo_root=pwd()` wording. A missing or malformed `assets_manifest.toml` propagates as a TOML or file error rather than an empty report. The `Any`-typed forwarding wrapper hides the keyword-only signature from method introspection.

## Provenance
Mapped from `src/SpaceAGORA.jl` line 552.
