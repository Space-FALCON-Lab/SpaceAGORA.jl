---
id: misc.precompile_workload_run_spaceagora_precompile_workload
label: _run_spaceagora_precompile_workload
kind: function
source:
  file: src/precompile_workload.jl
  symbol: _run_spaceagora_precompile_workload
  lines:
  - 41
  - 50
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: SpaceAGORA package scope supplying parse_parallel_profile, simulation_engine_config_from_env
    and run_simulation at precompile time.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: compiled_methods
  type: PrecompileCache
  units: n/a
  description: Native code cached into the package image for the orbit-only Mars aerobraking
    configuration exercised by the workload.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- precompile
- startup
- simulation
charts:
- misc
origin: agent
---

# _run_spaceagora_precompile_workload

## Purpose
`_run_spaceagora_precompile_workload` is the representative simulation executed while the package precompiles, so that the expensive type inference and code generation for the hot solver path is paid once at build time instead of on every CLI invocation. Without it, the first `run_simulation` call in a fresh process spends a long stretch compiling before any integration begins, which is painful both for interactive work and for short CI jobs.

## Model & Assumptions
The workload deliberately configures the cheapest run that still touches the real code paths: a Mars planet model, a three-body spacecraft built by `TelemetryVerification.make_three_body_spacecraft` with a 150 kg bus, two 2 kg panels and 15 kg of propellant, an initial orbit of 220 km apoapsis by 150 km periapsis at 28 degrees inclination, and a mission time of 5.0 seconds. Attitude propagation is disabled (`orientation_sim=false`), the Keplerian path is selected, gravity is the inverse-squared model, the atmosphere is the exponential model, and ephemerides are the simple model. The assumption is that these choices instantiate the same method signatures a production run needs, differing only in how long the integration runs.

## Design & Implementation
`_spaceagora_precompile_args` builds the configuration through `TelemetryVerification.make_example_config`, with `results=false`, `verbose=false` and a results directory under `tempdir()`. The environment the engine reads is not the process environment but the fixed dictionary `_SPACEAGORA_PRECOMPILE_ENV`, which pins `SPACEAGORA_PARALLEL_PROFILE` to `R2`, disables bundle saving with `SPACEAGORA_SAVE_BUNDLE=0`, and silences deprecated-config warnings — so precompilation is reproducible on any machine. The body warms `parse_parallel_profile("R2")`, builds the engine config with `simulation_engine_config_from_env`, then runs inside `mktempdir` and `cd` so any stray file lands in a directory that is deleted on exit. `run_simulation(..., return_solution=true)` forces the solution-returning method to be inferred. `@setup_workload` and `@compile_workload` wrap the call so PrecompileTools records the generated code into the package image.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | SpaceAGORA package scope supplying parse_parallel_profile, simulation_engine_config_from_env and run_simulation at precompile time. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `compiled_methods` | PrecompileCache | n/a | — | Native code cached into the package image for the orbit-only Mars aerobraking configuration exercised by the workload. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[parallel.profile_definitions_parse_parallel_profile|parse_parallel_profile]] · `callers` · call · `src/precompile_workload.jl:42-42`
- `callees` → [[simulation.run_simulation|run_simulation]] · `callers` · call · `src/precompile_workload.jl:47-47`
- `callees` → [[simx.engine_adapters_from_env_simulation_engine_config_from_env|simulation_engine_config_from_env]] · `callers` · call · `src/precompile_workload.jl:43-43`
- `callees` → [[simx.engine_execution_run_simulation|run_simulation]] · `callers` · call · `src/precompile_workload.jl:47-47`
- `callees` → [[spaceagora.precompile_workload__spaceagora_precompile_args|_spaceagora_precompile_args]] · `callers` · call · `src/precompile_workload.jl:44-44`
- `callees` → [[spaceagora.run_simulation|run_simulation]] · `callers` · call · `src/precompile_workload.jl:47-47`
<!-- vulcan:connections:end -->

## Limitations
Only the signatures this configuration reaches are cached: attitude dynamics, GRAM atmospheres, alternative parallel profiles, output-bundle writing and the plotting stack are all excluded by the chosen flags, so runs using them still compile on first call. Any failure inside the workload becomes a package precompilation failure rather than a runtime error, which makes the build sensitive to changes in the example-config helper. The workload writes into `tempdir()` during build, so a build environment with no writable temporary directory cannot precompile.

## Provenance
Mapped from `src/precompile_workload.jl:41-50`, with the configuration builder at lines 1-37 and the PrecompileTools wrapper at lines 52-54.
