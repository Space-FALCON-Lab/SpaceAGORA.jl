---
id: output.results_bundle
label: Results bundle (feather + manifest)
kind: external
inputs:
- id: bundle
  type: Arrow + TOML
  units: n/a
  description: Written by the results writer.
outputs: []
tags:
- master-flow
charts:
- master
origin: agent
---

# Results bundle (feather + manifest)

## Purpose
The machine-readable product: the same results table as an Arrow (feather) file alongside a TOML manifest recording schema version, creation time, mission time, step count, satellite count, orientation mode, and the size and SHA-256 of every file in the bundle.

## Design & Implementation
Written by `_write_results_bundle!` under the bundle prefix from `IOConfig`; the CSV is referenced in the manifest when `save_csv` is on. The manifest is what downstream tooling should trust when deciding whether a bundle is complete and unmodified.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `bundle` | Arrow + TOML | n/a | — | Written by the results writer. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.write_results|Write results & checkpoints]] · `bundle` → `bundle` · dataflow · `src/io/outputs/io_outputs.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Bundle writing is gated by `SPACEAGORA_SAVE_BUNDLE` and defaults on, so disabling it silently leaves only the CSV; the manifest's digests are recorded but not verified by any reader in this repository.
