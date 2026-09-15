---
id: output.verification_reports
label: Verification summary, errors & plots
kind: external
inputs:
- id: reports
  type: CSV + PNG
  units: n/a
  description: Written by the verification study.
outputs: []
tags:
- master-flow
charts:
- master
origin: agent
---

# Verification summary, errors & plots

## Purpose
What the telemetry study produces: a summary CSV with one row per scenario and channel — pass/fail against tolerance, NMAE, RMSE, bias — an errors CSV with per-sample residuals, and, when enabled, a directory of comparison plots.

## Design & Implementation
Written by `runner.jl` to the `out_summary` and `out_errors` paths (defaults under the output directory) and by `reporting.jl`'s `_generate_plots`; the calibration pass adds the fitted drag scale, reflectivity and bias to the summary. Under `--enforce` a failed tolerance exits non-zero.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `reports` | CSV + PNG | n/a | — | Written by the verification study. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.verification|Telemetry verification study]] · `reports` → `reports` · dataflow · `src/analysis/verification/telemetry_verification/runner.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Plot generation requires the plotting stack and is on by default, so a headless run should disable it; summary rows for the derived speed channels are only present when the manifest declares their tolerances.
