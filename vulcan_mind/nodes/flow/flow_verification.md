---
id: flow.verification
label: Telemetry verification study
kind: group
inputs:
- id: scenarios
  type: Vector{AbstractScenarioConfig}
  units: n/a
  description: Parsed manifest scenarios.
- id: telemetry_tables
  type: Arrow / CSV
  units: n/a
  description: Flight truth.
outputs:
- id: scenario_config
  type: SimulationConfiguration
  units: n/a
  description: Each scenario's run for the solve loop.
- id: reports
  type: CSV + plots
  units: n/a
  description: Summary and error tables and figures.
tags:
- master-flow
charts:
- master
origin: agent
opens: analysis
---

# Telemetry verification study

## Purpose
Scores the simulator against real missions: builds each scenario's configuration, runs it in a quick or full profile, compares apsis histories or state channels against telemetry within declared tolerances, optionally calibrates drag scale, reflectivity and bias, and writes the results.

## Design & Implementation
`scenario_builders.jl` assembles the planet, vehicle, effectors and atmosphere truth per scenario; `runner.jl` runs each through `run_simulation` with study tolerances, handling solver retries; `error_tables.jl` and `comparison_metrics.jl` compute NMAE, RMSE and decay diagnostics; `calibration.jl` grid-searches the fit parameters; `reporting.jl` writes CSVs and generates plots.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `scenarios` | Vector{AbstractScenarioConfig} | n/a | — | Parsed manifest scenarios. |
| in | `telemetry_tables` | Arrow / CSV | n/a | — | Flight truth. |
| out | `scenario_config` | SimulationConfiguration | n/a | — | Each scenario's run for the solve loop. |
| out | `reports` | CSV + plots | n/a | — | Summary and error tables and figures. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[flow.configure|Configure a run]] · `scenarios` → `scenarios` · dataflow · `src/analysis/verification/telemetry_verification/runner.jl`
- [[input.telemetry_truth|Flight telemetry (truth)]] · `telemetry_tables` → `telemetry_tables` · dataflow · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`

**Downstream**

- `reports` → [[output.verification_reports|Verification summary, errors & plots]] · `reports` · dataflow · `src/analysis/verification/telemetry_verification/runner.jl`
- `scenario_config` → [[flow.solve_loop|Solve loop]] · `scenario_config` · dataflow · `src/analysis/verification/telemetry_verification/runner.jl`
<!-- vulcan:connections:end -->

## Limitations
Telemetry is treated as truth; the entry-interface and tolerance settings interact with six environment variables, so the effective tolerances of a run are not echoed anywhere.
