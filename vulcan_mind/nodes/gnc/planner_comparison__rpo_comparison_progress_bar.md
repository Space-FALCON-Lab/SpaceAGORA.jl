---
id: gnc.planner_comparison__rpo_comparison_progress_bar
label: _rpo_comparison_progress_bar
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: _rpo_comparison_progress_bar
  lines:
  - 68
  - 68
inputs:
- id: frac
  type: Real
  units: n/a
  required: true
  description: Positional argument `frac`.
- id: width
  type: Integer
  units: n/a
  required: false
  description: Keyword argument `width` (default `30`).
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
  type: String
  units: n/a
  description: Return value of `_rpo_comparison_progress_bar`. Returns `"[" * repeat("=",
    filled) * repeat(" ", width - filled) * "]"`.
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

# _rpo_comparison_progress_bar

## Purpose
Renders a fixed-width ASCII progress bar string such as `[=========           ]` for the terminal progress line printed during comparison batches.

## Design & Implementation
`_rpo_comparison_progress_bar(frac::Real; width::Integer = 30)` clamps `frac` to [0, 1], computes `filled = Int(floor(f * width))`, and returns `"[" * repeat("=", filled) * repeat(" ", width - filled) * "]"`. Pure and allocation-light; used only by `_rpo_comparison_progress_line!`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `frac` | Real | n/a | yes | Positional argument `frac`. |
| in | `width` | Integer | n/a | no | Keyword argument `width` (default `30`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_rpo_comparison_progress_bar`. Returns `"[" * repeat("=", filled) * repeat(" ", width - filled) * "]"`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison__rpo_comparison_progress_line_bang|_rpo_comparison_progress_line!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:85-85`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:69-69`
<!-- vulcan:connections:end -->

## Limitations
`NaN` input passes through `clamp` as `NaN`, and `Int(floor(NaN))` throws `InexactError`; the caller guards `total <= 0` but not a NaN ratio. `width <= 0` yields `repeat` with a negative count error for the space padding when `filled` is 0.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 68.
