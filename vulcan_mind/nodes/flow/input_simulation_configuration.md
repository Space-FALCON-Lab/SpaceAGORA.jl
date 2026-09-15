---
id: input.simulation_configuration
label: SimulationConfiguration (script / API)
kind: external
inputs: []
outputs:
- id: args
  type: SimulationConfiguration
  units: n/a
  description: 'Fully specified run: settings, mission, environment, dynamics, GNC
    models, epoch, tolerances.'
tags:
- master-flow
charts:
- master
origin: agent
---

# SimulationConfiguration (script / API)

## Purpose
The programmatic entry: a `SimulationConfiguration` built in a Julia script or example, bundling simulation settings, mission definition, environment model, dynamics model, guidance/navigation/control models, initial epoch and integrator tolerances.

## Design & Implementation
Constructed with keyword constructors in `src/core/state/simulation_configuration.jl`, typically via `make_example_config` in the examples. Passed directly to `run_simulation`, or to `run_monte_carlo` inside a per-seed closure for campaigns. Nothing on disk is required beyond the data assets the chosen models reference.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| out | `args` | SimulationConfiguration | n/a | — | Fully specified run: settings, mission, environment, dynamics, GNC models, epoch, tolerances. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `args` → [[flow.configure|Configure a run]] · `api_args` · dataflow · `src/core/state/simulation_configuration.jl`
<!-- vulcan:connections:end -->

## Limitations
The struct is immutable, so every variation — a different atmosphere, a different mission length — is a rebuilt configuration; the example helpers copy field-by-field and silently drop `solver_config` in some paths.
