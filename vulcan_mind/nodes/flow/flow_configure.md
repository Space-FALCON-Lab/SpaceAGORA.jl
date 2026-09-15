---
id: flow.configure
label: Configure a run
kind: group
inputs:
- id: argv
  type: Vector{String}
  units: n/a
  description: Command-line tokens.
  required: false
- id: manifest
  type: TOML
  units: n/a
  description: Scenario manifest.
  required: false
- id: api_args
  type: SimulationConfiguration
  units: n/a
  description: Programmatic configuration.
  required: false
outputs:
- id: run_config
  type: SimulationConfiguration
  units: n/a
  description: A validated configuration ready for setup.
- id: scenarios
  type: Vector{AbstractScenarioConfig}
  units: n/a
  description: Typed verification scenarios.
- id: campaign_spec
  type: MonteCarloSpec / closure
  units: n/a
  description: Seeds and a per-seed configuration builder.
- id: asset_report
  type: AssetCheckReport
  units: n/a
  description: Presence and integrity of data assets.
tags:
- master-flow
charts:
- master
origin: agent
opens: cli
---

# Configure a run

## Purpose
Turns whatever the user provided — a command line, a scenario manifest, or a configuration built in a script — into validated, typed configuration objects, and checks that the data assets those configurations will need are present.

## Design & Implementation
`run_cli` parses subcommands and options; `_load_scenarios_from_manifest` parses TOML into scenario configs with every field validated; the example builders assemble a `SimulationConfiguration` from a planet, a three-body vehicle and force models; and `check_assets` walks the asset manifest and renders a report. Environment variables layered on top become typed snapshots at setup, never read on the hot path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `argv` | Vector{String} | n/a | no | Command-line tokens. |
| in | `manifest` | TOML | n/a | no | Scenario manifest. |
| in | `api_args` | SimulationConfiguration | n/a | no | Programmatic configuration. |
| out | `run_config` | SimulationConfiguration | n/a | — | A validated configuration ready for setup. |
| out | `scenarios` | Vector{AbstractScenarioConfig} | n/a | — | Typed verification scenarios. |
| out | `campaign_spec` | MonteCarloSpec / closure | n/a | — | Seeds and a per-seed configuration builder. |
| out | `asset_report` | AssetCheckReport | n/a | — | Presence and integrity of data assets. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[input.cli_args|CLI arguments]] · `argv` → `argv` · dataflow · `src/cli/spaceagora_cli.jl`
- [[input.scenario_manifest|Scenario manifest (TOML)]] · `manifest_toml` → `manifest` · dataflow · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
- [[input.simulation_configuration|SimulationConfiguration (script / API)]] · `args` → `api_args` · dataflow · `src/core/state/simulation_configuration.jl`

**Downstream**

- `asset_report` → [[output.asset_report|Asset check report]] · `asset_report` · dataflow · `src/cli/assets.jl`
- `campaign_spec` → [[flow.campaigns|Campaigns]] · `campaign_spec` · dataflow · `src/simulation/campaigns/monte_carlo.jl`
- `run_config` → [[flow.setup_run|Set up the run]] · `run_config` · dataflow · `src/simulation/engine/execution.jl`
- `run_config` → [[flow.vehicle|Spacecraft model]] · `run_config` · dataflow · `src/vehicle/spacecraft/model.jl`
- `scenarios` → [[flow.verification|Telemetry verification study]] · `scenarios` · dataflow · `src/analysis/verification/telemetry_verification/runner.jl`
<!-- vulcan:connections:end -->

## Limitations
Configuration validation is spread across constructors and parsers rather than one schema, so some inconsistent combinations — attitude quaternions under the wrong incidence mode, for instance — are caught only where the code happens to check.
