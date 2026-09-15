---
id: analysis.rpo_visualization_rpo_tracking_plot
label: rpo_tracking_plot
kind: function
source:
  file: src/analysis/visualization/rpo/rpo_visualization.jl
  symbol: rpo_tracking_plot
  lines:
  - 31
  - 31
inputs:
- id: t
  type: Any
  units: n/a
  required: true
  description: Positional argument `t`.
- id: err_norm
  type: Any
  units: n/a
  required: true
  description: Positional argument `err_norm`.
- id: title
  type: AbstractString
  units: n/a
  required: false
  description: Keyword argument `title` (default `"RPO Tracking Error"`).
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
  type: Plot
  units: n/a
  description: Return value of `rpo_tracking_plot`. Returns `Plot(`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- analysis
origin: agent
---

# rpo_tracking_plot

## Purpose

`rpo_tracking_plot(t, err_norm; title)` builds a single two-dimensional PlotlyJS figure of relative-navigation tracking error against time. It is the quick-look diagnostic for an RPO run: `t` is the time vector in seconds and `err_norm` the matching per-sample error norm in metres.

## Design & Implementation

The function is a one-expression constructor. It creates one `scatter` trace in `lines` mode named `tracking error` with `x=t` and `y=err_norm`, wraps it in a `Layout` carrying the caller-supplied `title` (default `"RPO Tracking Error"`), `xaxis_title="time (s)"` and `yaxis_title="error norm (m)"`, and returns the resulting `Plot`. Nothing is mutated, no file is written and no display is triggered, so the caller owns rendering and export.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `err_norm` | Any | n/a | yes | Positional argument `err_norm`. |
| in | `title` | AbstractString | n/a | no | Keyword argument `title` (default `"RPO Tracking Error"`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Plot | n/a | — | Return value of `rpo_tracking_plot`. Returns `Plot(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/visualization/rpo/rpo_visualization.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

The two arguments are passed straight through to PlotlyJS with no length check, so mismatched `t` and `err_norm` produce a malformed trace instead of an error. Axis labels hard-code seconds and metres regardless of the units actually supplied. Only a single series is supported, so comparing several runs requires composing traces outside this function.

## Provenance
Mapped from `src/analysis/visualization/rpo/rpo_visualization.jl` line 31.
