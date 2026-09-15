---
id: analysis.scenario_builders__make_orbit_args
label: _make_orbit_args
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _make_orbit_args
  lines:
  - 549
  - 549
inputs:
- id: cfg
  type: OrbitEventsScenarioConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: target_orbits
  type: Int
  units: n/a
  required: true
  description: Positional argument `target_orbits`.
- id: cd_scale
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `cd_scale` (default `1.0`).
- id: cr_override
  type: Union{Nothing, Float64}
  units: n/a
  required: false
  description: Keyword argument `cr_override` (default `nothing`).
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
  type: SimulationConfiguration
  units: n/a
  description: Return value of `_make_orbit_args`.
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

# _make_orbit_args

## Purpose
The top-level builder for an orbit-events scenario: resolves planet, initial state, vehicle, effectors and atmosphere into a runnable configuration for a target orbit count.

## Design & Implementation
Resolves the planet, the initial condition, the spacecraft and the effector tuple with optional drag and reflectivity overrides, and the density model. It computes mission time as `target_orbits` times the two-body period, calls `make_example_config` with orientation off, non-Keplerian propagation and the scenario's entry interface, then layers on orbit-counted termination, the wind flag and campaign manoeuvres.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | OrbitEventsScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `target_orbits` | Int | n/a | yes | Positional argument `target_orbits`. |
| in | `cd_scale` | Float64 | n/a | no | Keyword argument `cd_scale` (default `1.0`). |
| in | `cr_override` | Union{Nothing, Float64} | n/a | no | Keyword argument `cr_override` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationConfiguration | n/a | — | Return value of `_make_orbit_args`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_single_scenario|_run_single_scenario]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:161-161`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[analysis.scenario_builders__make_spacecraft|_make_spacecraft]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:559-559`
- `callees` → [[analysis.scenario_builders__period_seconds|_period_seconds]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:569-569`
- `callees` → [[analysis.scenario_builders__planet_from_name|_planet_from_name]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:555-555`
- `callees` → [[analysis.scenario_builders__scenario_density_model|_scenario_density_model]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:567-567`
- `callees` → [[analysis.scenario_builders__scenario_dynamic_effectors|_scenario_dynamic_effectors]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:560-560`
- `callees` → [[analysis.scenario_builders__scenario_initial_condition|_scenario_initial_condition]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:557-557`
- `callees` → [[analysis.scenario_builders__with_campaign_maneuvers|_with_campaign_maneuvers]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:586-586`
- `callees` → [[analysis.scenario_builders__with_environment_wind|_with_environment_wind]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:585-585`
- `callees` → [[analysis.scenario_builders__with_orbit_mission|_with_orbit_mission]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:584-584`
- `callees` → [[envana.ana_example_support_make_example_config|make_example_config]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:572-572`
<!-- vulcan:connections:end -->

## Limitations
Each layer rebuilds the immutable configuration, so the returned object is the fourth full copy; `solver_config` is lost at the wind and manoeuvre layers because those rebuilders do not forward it.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 549.
