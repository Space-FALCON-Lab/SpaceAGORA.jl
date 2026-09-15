---
id: input.scenario_manifest
label: Scenario manifest (TOML)
kind: external
inputs: []
outputs:
- id: manifest_toml
  type: TOML file
  units: n/a
  description: Scenario table listing planet, spacecraft, events, tolerances, atmosphere
    truth and telemetry paths.
tags:
- master-flow
charts:
- master
origin: agent
---

# Scenario manifest (TOML)

## Purpose
The TOML manifest that defines a telemetry verification study: one table per scenario naming the planet, vehicle geometry, initial orbit, force models, atmosphere truth source, comparison tolerances and the telemetry files to score against.

## Design & Implementation
Parsed by `_load_scenarios_from_manifest` in `manifest_parsing.jl` into typed `OrbitEventsScenarioConfig` and `TimeAlignedScenarioConfig` records. Relative paths are resolved against the repository root, and every field is validated at parse time so a typo fails before any simulation runs.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `manifest_toml` | TOML file | n/a | — | Scenario table listing planet, spacecraft, events, tolerances, atmosphere truth and telemetry paths. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `manifest_toml` → [[flow.configure|Configure a run]] · `manifest` · dataflow · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`
<!-- vulcan:connections:end -->

## Limitations
The manifest is the only route into the verification study; there is no programmatic constructor for scenario configs, so a one-off scenario still has to be written as TOML.
