---
id: cli.assets_render_asset_manifest
label: render_asset_manifest
kind: function
source:
  file: src/cli/assets.jl
  symbol: render_asset_manifest
  lines:
  - 119
  - 119
inputs:
- id: entries
  type: Vector{AssetManifestEntry}
  units: n/a
  required: true
  description: Positional argument `entries`.
- id: io
  type: IO
  units: n/a
  required: false
  description: Keyword argument `io` (default `stdout`).
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
  type: Nothing
  units: n/a
  description: Return value of `render_asset_manifest`. Returns `nothing`.
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

# render_asset_manifest

## Purpose
Prints a human-readable listing of every `AssetManifestEntry` so users can see what asset roots the repository declares, independent of whether they are installed. It backs the CLI's manifest inspection command.

## Design & Implementation
Signature `render_asset_manifest(entries::Vector{AssetManifestEntry}; io::IO=stdout)`. It writes a header `SpaceAGORA asset manifest` and `entries=<count>`, then for each entry prints `- <name>` followed by indented lines for `scope`, `kind`, `relative_path`, `required` and `licensing`, and a `detail` line only when `entry.detail` is non-empty. Output is plain text via `println(io, ...)`; the function returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `entries` | Vector{AssetManifestEntry} | n/a | yes | Positional argument `entries`. |
| in | `io` | IO | n/a | no | Keyword argument `io` (default `stdout`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `render_asset_manifest`. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.run_cli|run_cli]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:197-197`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/cli/assets.jl:120-120`
<!-- vulcan:connections:end -->

## Limitations
The format is fixed and not machine-parseable (no JSON or TOML option), and there is no column alignment or sorting, so long manifests are read in declaration order only. Nothing is escaped, so newlines inside `detail` break the two-space indentation convention. The function does not resolve paths, so it cannot show where an entry would live on this machine; use `render_asset_report` for that.

## Provenance
Mapped from `src/cli/assets.jl` line 119.
