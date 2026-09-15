---
id: gnc.planner_comparison_rpo_comparison_metric_specs
label: rpo_comparison_metric_specs
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_comparison_metric_specs
  lines:
  - 690
  - 690
inputs:
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
  type: AbstractArray
  units: n/a
  description: Return value of `rpo_comparison_metric_specs`. Returns `[`.
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

# rpo_comparison_metric_specs

## Purpose
Declares the ordered list of nine metrics rendered in the summary plot grid, each as a `(field, display label, units)` triple, giving the plotting and export code one place to look up names and units.

## Design & Implementation
Returns a literal `Vector` of tuples: `(:success_pct, "Success", "success (%)")`, `(:planner_compute_time, "Planner Runtime", "s")`, `(:fuel_used_pct, "Fuel", "%")`, `(:control_effort_total, "Control Effort", "m/s")`, `(:thrust_saturation_fraction, "Thruster Sat.", "fraction")`, `(:min_clearance, "Min Clearance", "m")`, `(:final_pos_error, "Goal Error", "m")`, `(:actual_travel_duration, "Travel Duration", "s")`, `(:keepout_violations, "Keepout Violations", "count")`. `rpo_comparison_metric_summary_plot` lays these out in a fixed 3 x 3 grid.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractArray | n/a | — | Return value of `rpo_comparison_metric_specs`. Returns `[`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_comparison_metric_summary_plot|rpo_comparison_metric_summary_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:729-729`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The plot grid is hard-coded to 3 rows by 3 columns, so adding a tenth spec would place a subplot outside the layout domain. `:success_pct` is a synthetic metric that exists only through `rpo_comparison_metric_value`; it is not a column of the CSV rows.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 690.
