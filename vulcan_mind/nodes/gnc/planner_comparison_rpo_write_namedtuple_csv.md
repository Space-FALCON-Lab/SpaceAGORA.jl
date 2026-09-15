---
id: gnc.planner_comparison_rpo_write_namedtuple_csv
label: rpo_write_namedtuple_csv
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_write_namedtuple_csv
  lines:
  - 658
  - 658
inputs:
- id: path
  type: AbstractString
  units: n/a
  required: true
  description: Positional argument `path`.
- id: rows
  type: Any
  units: n/a
  required: true
  description: Positional argument `rows`.
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
  description: Return value of `rpo_write_namedtuple_csv`. Returns `path`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# rpo_write_namedtuple_csv

## Purpose
Writes a vector of NamedTuple rows to a CSV file using the first row's field order as the header, without depending on CSV.jl.

## Design & Implementation
`rpo_write_namedtuple_csv(path::AbstractString, rows)` returns `path` immediately for empty input; otherwise creates the parent directory with `mkpath`, takes `names = collect(keys(first(rows)))`, and writes the header followed by one line per row where each value is `string(value)` with commas replaced by semicolons. Returns `path`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | AbstractString | n/a | yes | Positional argument `path`. |
| in | `rows` | Any | n/a | yes | Positional argument `rows`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_write_namedtuple_csv`. Returns `path`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_write_planner_comparison_outputs|rpo_write_planner_comparison_outputs]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1157-1157`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:663-663`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:668-668`
<!-- vulcan:connections:end -->

## Limitations
Values are not quoted, so embedded newlines or double quotes corrupt the file; only commas are sanitised. Rows with fields missing from the first row raise an error via `getproperty`, and extra fields in later rows are silently dropped. Floats are written with Julia's default `string` formatting (for example `1.0e-5`), which some spreadsheet importers mis-parse.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 658.
