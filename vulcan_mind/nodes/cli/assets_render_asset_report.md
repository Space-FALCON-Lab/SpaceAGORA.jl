---
id: cli.assets_render_asset_report
label: render_asset_report
kind: function
source:
  file: src/cli/assets.jl
  symbol: render_asset_report
  lines:
  - 106
  - 106
inputs:
- id: report
  type: AssetCheckReport
  units: n/a
  required: true
  description: Positional argument `report`.
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
  description: Return value of `render_asset_report`. Returns `nothing`.
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

# render_asset_report

## Purpose
Prints an `AssetCheckReport` as a plain-text status listing, classifying each asset as `available`, `missing-required` or `missing-optional` so a user can immediately see whether the repository is ready for baseline runs or needs licensed data installed.

## Design & Implementation
Signature `render_asset_report(report::AssetCheckReport; io::IO=stdout)`. After the header lines `SpaceAGORA asset check` and `repo_root=<path>`, it loops over `report.items` and computes `status` as `"available"` when `item.available`, otherwise `"missing-required"` if `item.required` else `"missing-optional"`. Each item produces four lines: `- <name>: <status>` plus indented `scope`, `path` and `detail`. All output goes through `println(io, ...)` and the function returns `nothing`. `setup_open_assets` calls it as the final step of its summary.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `report` | AssetCheckReport | n/a | yes | Positional argument `report`. |
| in | `io` | IO | n/a | no | Keyword argument `io` (default `stdout`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `render_asset_report`. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[cli.assets_setup_open_assets|setup_open_assets]] · `callees` → `callers` · call · `src/cli/assets.jl:154-154`
- [[cli.run_cli|run_cli]] · `callees` → `callers` · call · `src/cli/spaceagora_cli.jl:194-194`
- [[spaceagora.spaceagora_check_assets|check_assets]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:555-555`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/cli/assets.jl:107-107`
<!-- vulcan:connections:end -->

## Limitations
There is no aggregate summary or non-zero return signal when required assets are missing, so scripts must parse the text or inspect the report themselves. The `detail` string is printed even when empty, yielding a bare `detail:` line. Output ordering and wording are not versioned, so any tooling that scrapes it is fragile to formatting changes.

## Provenance
Mapped from `src/cli/assets.jl` line 106.
