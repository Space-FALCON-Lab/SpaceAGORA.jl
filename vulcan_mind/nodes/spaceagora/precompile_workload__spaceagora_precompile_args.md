---
id: spaceagora.precompile_workload__spaceagora_precompile_args
label: _spaceagora_precompile_args
kind: function
source:
  file: src/precompile_workload.jl
  symbol: _spaceagora_precompile_args
  lines:
  - 1
  - 1
inputs:
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
  type: TelemetryVerification.make_example_config
  units: n/a
  description: Return value of `_spaceagora_precompile_args`. Returns `TelemetryVerification.make_example_config(`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- spaceagora
charts:
- spaceagora
origin: agent
---

# _spaceagora_precompile_args

## Purpose
Builds the concrete argument set the precompile workload runs against, so the package's hot paths are compiled at build time rather than on a user's first simulation.

## Design & Implementation
Constructs a Mars planet model, then a three-body spacecraft via `TelemetryVerification.make_three_body_spacecraft` with explicit bus and panel dimensions, masses and a 0.7 m panel offset. The initial condition is an elliptical Mars orbit with apoapsis at `Rp_e + 220 km`, periapsis at `Rp_e + 150 km`, inclination 28 degrees, and fixed argument of periapsis, RAAN and true anomaly. The result is handed to `make_example_config`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | TelemetryVerification.make_example_config | n/a | — | Return value of `_spaceagora_precompile_args`. Returns `TelemetryVerification.make_example_config(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[misc.precompile_workload_run_spaceagora_precompile_workload|_run_spaceagora_precompile_workload]] · `callees` → `callers` · call · `src/precompile_workload.jl:44-44`

**Downstream**

- `callees` → [[analysis.example_support_make_three_body_spacecraft|make_three_body_spacecraft]] · `callers` · call · `src/precompile_workload.jl:4-4`
- `callees` → [[core.simulation_configuration_initialtime|InitialTime]] · `callers` · call · `src/precompile_workload.jl:26-26`
- `callees` → [[envana.ana_example_support_make_example_config|make_example_config]] · `callers` · call · `src/precompile_workload.jl:22-22`
- `callees` → [[envana.env_simple_ephemerides_simpleephemeridesmodel|SimpleEphemeridesModel]] · `callers` · call · `src/precompile_workload.jl:29-29`
- `callees` → [[environment.density_models_exponentialatmospheremodel|ExponentialAtmosphereModel]] · `callers` · call · `src/precompile_workload.jl:28-28`
- `callees` → [[environment.gravity_models_inversesquaredgravitymodel|InverseSquaredGravityModel]] · `callers` · call · `src/precompile_workload.jl:27-27`
- `callees` → [[environment.planets_mars|Mars]] · `callers` · call · `src/precompile_workload.jl:3-3`
- `callees` → [[vehicle.model_initialcondition|InitialCondition]] · `callers` · call · `src/precompile_workload.jl:10-10`
<!-- vulcan:connections:end -->

## Limitations
Every value is hard-coded, so the workload exercises one Mars aerobraking geometry only; code paths reached solely by other planets or vehicle configurations stay uncompiled.

## Provenance
Mapped from `src/precompile_workload.jl` line 1.
