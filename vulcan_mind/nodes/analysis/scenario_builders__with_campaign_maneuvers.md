---
id: analysis.scenario_builders__with_campaign_maneuvers
label: _with_campaign_maneuvers
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _with_campaign_maneuvers
  lines:
  - 398
  - 398
inputs:
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
- id: cfg
  type: OrbitEventsScenarioConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  description: Return value of `_with_campaign_maneuvers`.
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

# _with_campaign_maneuvers

## Purpose
Attaches the thruster and the campaign manoeuvre guidance to a scenario configuration so replayed apoapsis burns fire at the flight orbit numbers.

## Design & Implementation
Returns the input unchanged when there are no manoeuvres. It requires exactly one spacecraft and positive thrust, Isp, guidance rate and control rate, raising `ArgumentError` naming the scenario otherwise. A `BaseThrusterModel` is built with per-satellite vectors and sentinel burn times of minus one. When `maneuver_replay_scale_mode` is `flight_apoapsis_ratio` the flight apoapsis altitudes are converted to radii with `Rp_e` for the guidance model's diagnostic scaling. The result replaces the guidance and control models with a single `AerobrakingCampaignPropulsiveManeuverGuidanceModel` and the thruster at their configured rates.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `cfg` | OrbitEventsScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationConfiguration | n/a | — | Return value of `_with_campaign_maneuvers`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__make_orbit_args|_make_orbit_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:586-586`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[analysis.scenario_builders__has_campaign_maneuvers|_has_campaign_maneuvers]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:399-399`
- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:423-423`
- `callees` → [[gncz.thruster_guidance_models_aerobrakingcampaignpropulsivemaneuverguidancemodel|AerobrakingCampaignPropulsiveManeuverGuidanceModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:428-428`
- `callees` → [[parcore.simulation_configuration_simulationconfiguration|SimulationConfiguration]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:433-433`
- `callees` → [[vehicle.model_controlmodel|ControlModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:444-444`
- `callees` → [[vehicle.model_guidancemodel|GuidanceModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:439-439`
- `callees` → [[vehx.actuators_thruster_models_basethrustermodel|BaseThrusterModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:410-410`
<!-- vulcan:connections:end -->

## Limitations
Single-spacecraft only, enforced by an explicit check; the equatorial-radius conversion is acknowledged in the code as an approximation whose error cancels in the flight-to-sim ratio at the 0.2 percent level.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 398.
