---
id: output.planner_comparison_plots
label: RPO planner comparison figures
kind: external
inputs:
- id: planner_plots
  type: PNG / HTML
  units: n/a
  description: Saved by the planner comparison batch.
outputs: []
tags:
- master-flow
charts:
- master
origin: agent
---

# RPO planner comparison figures

## Purpose
Figures from the RPO planner comparison study: failed-path families, metric summaries, path families per planner and PSO cost-versus-iteration curves, saved when a comparison batch is run with plotting on.

## Design & Implementation
Written with `PlotlyJS.savefig` from `planner_comparison.jl` after a batch of PSO, RRT and optimiser runs against the station geometry, alongside a serialized batch record.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planner_plots` | PNG / HTML | n/a | — | Saved by the planner comparison batch. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.gnc|Guidance, navigation & control]] · `planner_plots` → `planner_plots` · dataflow · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
These are produced by a study entrypoint in the GNC module rather than by a simulation run, so they do not go through the results writer or its atomic-write discipline.
