---
id: flow.write_results
label: Write results & checkpoints
kind: group
inputs:
- id: saved_values
  type: SavedValues / DataFrame
  units: n/a
  description: Collected time series.
- id: checkpoint_state
  type: (t, u, solver_mode)
  units: n/a
  description: Integrator state at a segment boundary.
  required: false
outputs:
- id: csv
  type: simulation_results.csv
  units: n/a
  description: Row-per-save-instant table.
- id: bundle
  type: feather + manifest
  units: n/a
  description: Arrow table with a TOML manifest.
- id: checkpoint
  type: .jls + manifest
  units: n/a
  description: Resumable state.
tags:
- master-flow
charts:
- master
origin: agent
opens: io
---

# Write results & checkpoints

## Purpose
Persists what a run produced: the results table as CSV and as an Arrow bundle with an integrity manifest, and — when checkpointing is on — the integrator state at each segment boundary so a killed run can resume.

## Design & Implementation
`io_outputs.jl` builds the results `DataFrame` from `SaveField` columns, writes the CSV with collision-safe naming, and writes the feather file plus a TOML manifest recording size, SHA-256 and run metadata; `io_serialization.jl` writes checkpoints atomically as a serialized payload plus manifest and loads them on resume. All writes go through `_atomic_write_file` so a crash never leaves a partial file.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `saved_values` | SavedValues / DataFrame | n/a | — | Collected time series. |
| in | `checkpoint_state` | (t, u, solver_mode) | n/a | no | Integrator state at a segment boundary. |
| out | `csv` | simulation_results.csv | n/a | — | Row-per-save-instant table. |
| out | `bundle` | feather + manifest | n/a | — | Arrow table with a TOML manifest. |
| out | `checkpoint` | .jls + manifest | n/a | — | Resumable state. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.callbacks|Integration callbacks]] · `saved_rows` → `checkpoint_state` · dataflow · `src/io/serialization/io_serialization.jl`
- [[flow.solve_loop|Solve loop]] · `saved_values` → `saved_values` · dataflow · `src/io/outputs/io_outputs.jl`

**Downstream**

- `bundle` → [[output.results_bundle|Results bundle (feather + manifest)]] · `bundle` · dataflow · `src/io/outputs/io_outputs.jl`
- `checkpoint` → [[output.checkpoint|Checkpoint (.jls + manifest)]] · `checkpoint` · dataflow · `src/io/serialization/io_serialization.jl`
- `csv` → [[output.results_csv|simulation_results.csv]] · `csv` · dataflow · `src/io/outputs/io_outputs.jl`
<!-- vulcan:connections:end -->

## Limitations
Checkpoints are trusted on load without verifying the manifest's digest; Julia's serialization format ties a checkpoint to the Julia version that wrote it.
