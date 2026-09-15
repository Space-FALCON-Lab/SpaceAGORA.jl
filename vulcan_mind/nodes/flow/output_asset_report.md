---
id: output.asset_report
label: Asset check report
kind: external
inputs:
- id: asset_report
  type: AssetCheckReport
  units: n/a
  description: Rendered by the CLI.
outputs: []
tags:
- master-flow
charts:
- master
origin: agent
---

# Asset check report

## Purpose
A rendered table stating, for every asset in the data manifest — kernels, GRAM data, telemetry, station geometry — whether it is present, its size and whether its digest matches, so a new machine can be checked before any simulation is attempted.

## Design & Implementation
Produced by `check_assets` walking `data/assets_manifest.toml` and `render_asset_report` printing to the CLI's output stream; invoked by the `check-assets` subcommand.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `asset_report` | AssetCheckReport | n/a | — | Rendered by the CLI. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.configure|Configure a run]] · `asset_report` → `asset_report` · dataflow · `src/cli/assets.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It reports presence and integrity, not sufficiency — a scenario that needs a kernel not listed in the manifest still fails at run time.
