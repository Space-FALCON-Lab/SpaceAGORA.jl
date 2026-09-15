---
id: gnc.planner_comparison__rpo_comparison_progress_line_bang
label: _rpo_comparison_progress_line!
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: _rpo_comparison_progress_line!
  lines:
  - 75
  - 75
inputs:
- id: completed
  type: Integer
  units: n/a
  required: true
  description: Positional argument `completed`.
- id: total
  type: Integer
  units: n/a
  required: true
  description: Positional argument `total`.
- id: planner
  type: Symbol
  units: n/a
  required: false
  description: Keyword argument `planner` (default `:unknown`).
- id: case_label
  type: AbstractString
  units: n/a
  required: false
  description: Keyword argument `case_label` (default `""`).
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
  description: Return value of `_rpo_comparison_progress_line!`; mutates `completed`
    in place. Returns `nothing`.
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

# _rpo_comparison_progress_line!

## Purpose
Prints or overwrites a single-line terminal progress indicator for the comparison batch, showing the bar, percentage, completed/total counts, and the current planner label and case label.

## Design & Implementation
`_rpo_comparison_progress_line!(completed, total; planner = :unknown, case_label = "")` computes `frac = total <= 0 ? 1.0 : completed / total`, formats the percentage to one decimal padded to width 5, appends `rpo_comparison_planner_label(planner)` and `/ case_label` when supplied, and writes `"\r" * rpad(line, 120)` to stdout followed by `flush(stdout)`. When `completed >= total` it emits a newline so the final state persists. Returns `nothing`. The `!` suffix denotes the stdout side effect rather than argument mutation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `completed` | Integer | n/a | yes | Positional argument `completed`. |
| in | `total` | Integer | n/a | yes | Positional argument `total`. |
| in | `planner` | Symbol | n/a | no | Keyword argument `planner` (default `:unknown`). |
| in | `case_label` | AbstractString | n/a | no | Keyword argument `case_label` (default `""`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_rpo_comparison_progress_line!`; mutates `completed` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.planner_comparison_rpo_run_planner_comparison_batch|rpo_run_planner_comparison_batch]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:557-557`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:88-88`
- `callees` → [[gnc.planner_comparison__rpo_comparison_progress_bar|_rpo_comparison_progress_bar]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:85-85`
- `callees` → [[gnc.trajectory_optimizers_rpo_comparison_planner_label|rpo_comparison_planner_label]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:83-83`
<!-- vulcan:connections:end -->

## Limitations
Lines longer than 120 characters (long case labels) are not truncated and will wrap, breaking the carriage-return overwrite. Output goes to `stdout` unconditionally; there is no logger integration, and in non-TTY contexts the `\r` sequences accumulate in captured logs. `planner` not equal to `:unknown` triggers normalisation, which throws for unsupported symbols.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 75.
