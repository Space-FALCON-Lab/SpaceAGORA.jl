---
id: core.simulation_configuration_missionconfiguration
label: MissionConfiguration
kind: struct
source:
  file: src/core/state/simulation_configuration.jl
  symbol: MissionConfiguration
  lines:
  - 139
  - 139
inputs:
- id: mission_type
  type: MissionType
  units: n/a
  required: true
  description: Field `mission_type`.
- id: keplerian
  type: Bool
  units: n/a
  required: true
  description: Field `keplerian`.
- id: number_of_orbits
  type: Int
  units: n/a
  required: true
  description: Field `number_of_orbits`.
- id: mission_time
  type: Float64
  units: n/a
  required: true
  description: Field `mission_time`.
- id: orientation_sim
  type: Bool
  units: n/a
  required: true
  description: Field `orientation_sim`.
- id: num_steps_to_save
  type: Int
  units: n/a
  required: true
  description: Field `num_steps_to_save`.
- id: data_rate
  type: Float64
  units: n/a
  required: true
  description: Field `data_rate`.
- id: mission_type_2
  type: MissionType,
  units: n/a
  required: true
  description: Field `mission_type`.
- id: keplerian_2
  type: Bool,
  units: n/a
  required: true
  description: Field `keplerian`.
- id: number_of_orbits_2
  type: Integer,
  units: n/a
  required: true
  description: Field `number_of_orbits`.
- id: mission_time_2
  type: Real,
  units: n/a
  required: true
  description: Field `mission_time`.
- id: orientation_sim_2
  type: Bool,
  units: n/a
  required: true
  description: Field `orientation_sim`.
- id: num_steps_to_save_2
  type: Integer,
  units: n/a
  required: true
  description: Field `num_steps_to_save`.
- id: data_rate_2
  type: Float64
  units: n/a
  required: false
  description: Field `data_rate` (default `10.0`).
- id: data_rate_3
  type: Any
  units: n/a
  required: true
  description: Field `data_rate`.
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
  type: MissionConfiguration
  units: n/a
  description: Constructed `MissionConfiguration`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# MissionConfiguration

## Purpose
`MissionConfiguration` defines how long a simulation runs and what is propagated: the termination criterion (`MissionTime` in seconds or `MissionOrbits` by orbit count), whether the drag passage is integrated as a separate Keplerian step, whether attitude dynamics are simulated, and the output sampling cadence. It is validated on construction so the engine can rely on positive values.

## Design & Implementation
An immutable struct with fields `mission_type::MissionType`, `keplerian::Bool`, `number_of_orbits::Int`, `mission_time::Float64` (s), `orientation_sim::Bool`, `num_steps_to_save::Int` and `data_rate::Float64` (s, used as `saveat`). The positional inner constructor throws `ArgumentError` when `number_of_orbits <= 0`, `mission_time <= 0`, `num_steps_to_save <= 0` or `data_rate <= 0.0`, then converts integers with `Int(...)` and time with `Float64(...)`. A keyword outer constructor supplies defaults (`mission_type = MissionTime`, `keplerian = true`, `number_of_orbits = 1`, `mission_time = 90*60*20*10` = 1,080,000 s, `orientation_sim = false`, `num_steps_to_save = 1000`, `data_rate = 10.0`) and routes `mission_type` through `_parse_mission_type` so strings and symbols are accepted.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `mission_type` | MissionType | n/a | yes | Field `mission_type`. |
| in | `keplerian` | Bool | n/a | yes | Field `keplerian`. |
| in | `number_of_orbits` | Int | n/a | yes | Field `number_of_orbits`. |
| in | `mission_time` | Float64 | n/a | yes | Field `mission_time`. |
| in | `orientation_sim` | Bool | n/a | yes | Field `orientation_sim`. |
| in | `num_steps_to_save` | Int | n/a | yes | Field `num_steps_to_save`. |
| in | `data_rate` | Float64 | n/a | yes | Field `data_rate`. |
| in | `mission_type_2` | MissionType, | n/a | yes | Field `mission_type`. |
| in | `keplerian_2` | Bool, | n/a | yes | Field `keplerian`. |
| in | `number_of_orbits_2` | Integer, | n/a | yes | Field `number_of_orbits`. |
| in | `mission_time_2` | Real, | n/a | yes | Field `mission_time`. |
| in | `orientation_sim_2` | Bool, | n/a | yes | Field `orientation_sim`. |
| in | `num_steps_to_save_2` | Integer, | n/a | yes | Field `num_steps_to_save`. |
| in | `data_rate_2` | Float64 | n/a | no | Field `data_rate` (default `10.0`). |
| in | `data_rate_3` | Any | n/a | yes | Field `data_rate`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | MissionConfiguration | n/a | — | Constructed `MissionConfiguration`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.example_support__example_smoke_args|_example_smoke_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:27-27`
- [[analysis.scenario_builders__with_orbit_mission|_with_orbit_mission]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:459-459`
- [[analysis.scenario_builders__with_study_settings|_with_study_settings]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:672-672`
- [[envana.ana_example_support_make_example_config|make_example_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/example_support.jl:143-143`
- [[parcore.simulation_configuration_simulationconfiguration|SimulationConfiguration]] · `callees` → `callers` · call · `src/core/state/simulation_configuration.jl:238-238`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/core/state/simulation_configuration.jl:166-166`
- `callees` → [[core.simulation_configuration__parse_mission_type|_parse_mission_type]] · `callers` · call · `src/core/state/simulation_configuration.jl:184-184`
<!-- vulcan:connections:end -->

## Limitations
Both `number_of_orbits` and `mission_time` must be positive even though only one is meaningful for a given `mission_type`, so a caller cannot signal 'unused' with zero. The 1,080,000 s default (12.5 days) is an arbitrary literal expressed as a product. `data_rate` is typed `Float64` in the keyword constructor, so passing an `Int` fails with `MethodError` rather than converting. There is no upper bound on `num_steps_to_save`, which governs in-memory buffering.

## Provenance
Mapped from `src/core/state/simulation_configuration.jl` line 139.
