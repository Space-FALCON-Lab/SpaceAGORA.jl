---
id: output.results_csv
label: simulation_results.csv
kind: external
inputs:
- id: csv
  type: CSV file
  units: n/a
  description: Written by the results writer.
outputs: []
tags:
- master-flow
charts:
- master
origin: agent
---

# simulation_results.csv

## Purpose
The primary human-readable product of a run: one row per saved instant with time and every configured `SaveField` column, per satellite — position, velocity, mass, density, heat rate, forces, attitude and any custom fields.

## Design & Implementation
Written to `results_directory/simulation_results.csv` by `_write_results_csv!`, with a collision-safe alternate name when a file already exists and `generate_filenames` is off. Column names follow `sc<i>_<field>_<component>`, which the telemetry error tables read back by name.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `csv` | CSV file | n/a | — | Written by the results writer. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.write_results|Write results & checkpoints]] · `csv` → `csv` · dataflow · `src/io/outputs/io_outputs.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Wide constellations produce very wide tables; the CSV carries no schema or units beyond the column names, so consumers rely on the conventions in `save_fields.jl`.
