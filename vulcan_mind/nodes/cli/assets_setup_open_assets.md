---
id: cli.assets_setup_open_assets
label: setup_open_assets
kind: function
source:
  file: src/cli/assets.jl
  symbol: setup_open_assets
  lines:
  - 134
  - 134
inputs:
- id: repo_root
  type: String
  units: n/a
  required: false
  description: Keyword argument `repo_root` (default `REPO_ROOT`).
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
  type: Any
  units: n/a
  description: Return value of `setup_open_assets`. Returns `report`.
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

# setup_open_assets

## Purpose
Guides a new user through the asset situation: it explains that the baseline no-GRAM mode needs no downloads, lists which manifest entries are open versus licensed and user-provided, and finishes with the live availability report. It is the CLI entry for an open-source-first setup flow.

## Design & Implementation
Keyword arguments `repo_root::String=REPO_ROOT` and `io::IO=stdout`. It loads `entries = load_asset_manifest(; repo_root)` and `report = check_assets(; repo_root)`, recomputes `manifest_path` as `joinpath(repo_root, "data", "assets_manifest.toml")`, then prints a banner, the sentence about no downloads for baseline mode, and the manifest path. Two passes over `entries` follow: entries whose `licensing != "licensed-external"` are printed under `Baseline/open entries:` and those equal to `"licensed-external"` under `Licensed external entries remain user-provided:`, each as `- <name>: <relative_path>`. Finally `render_asset_report(report; io)` is called and the `AssetCheckReport` is returned.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `repo_root` | String | n/a | no | Keyword argument `repo_root` (default `REPO_ROOT`). |
| in | `io` | IO | n/a | no | Keyword argument `io` (default `stdout`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `setup_open_assets`. Returns `report`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.run_cli|run_cli]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:200-200`

**Downstream**

- `callees` → [[cli.assets_load_asset_manifest|load_asset_manifest]] · `callers` · call · `src/cli/assets.jl:135-135`
- `callees` → [[cli.assets_render_asset_report|render_asset_report]] · `callers` · call · `src/cli/assets.jl:154-154`
- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/cli/assets.jl:138-138`
- `callees` → [[misc.assets_check_assets|check_assets]] · `callers` · feedback · `src/cli/assets.jl:136-136`
- `callees` → [[spaceagora.spaceagora_check_assets|check_assets]] · `callers` · feedback · `src/cli/assets.jl:136-136`
- `callees` → [[spaceagora.spaceagora_render_asset_report|render_asset_report]] · `callers` · feedback · `src/cli/assets.jl:154-154`
<!-- vulcan:connections:end -->

## Limitations
The manifest is parsed twice (once directly and once inside `check_assets`), and the `manifest_path` default is duplicated rather than passed through, so a custom manifest location cannot be used. Despite its name the function performs no setup actions: it creates no directories and downloads nothing. The licensing split relies on the exact string `"licensed-external"`; any other licensing label is treated as open.

## Provenance
Mapped from `src/cli/assets.jl` line 134.
