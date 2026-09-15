---
id: analysis.scenario_builders__make_required_gram_density_model
label: _make_required_gram_density_model
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _make_required_gram_density_model
  lines:
  - 316
  - 316
inputs:
- id: planet_name
  type: String
  units: n/a
  required: true
  description: Positional argument `planet_name`.
- id: initial_time
  type: InitialTime
  units: n/a
  required: true
  description: Positional argument `initial_time`.
- id: truth
  type: AtmosphereTruthConfig
  units: n/a
  required: true
  description: Positional argument `truth`.
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
  type: Base.invokelatest
  units: n/a
  description: Return value of `_make_required_gram_density_model`. Returns `Base.invokelatest(`
    or `offline_model`.
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

# _make_required_gram_density_model

## Purpose
Constructs the native GRAM atmosphere for a scenario, forwarding every truth-configuration knob, with a library-missing fallback to the offline surrogate.

## Design & Implementation
Calls `GRAMAtmosphereModel` through `Base.invokelatest` with planet, initial time, seed, minimum relative step, perturbation scales and the seven Mars-specific settings. On any exception it checks `_is_gram_library_missing_error`; if so and `_try_libraryless_gram_surrogate` returns a model, it warns and returns that. Every other failure is rethrown as an `ErrorException` naming the planet, initial time and original message.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet_name` | String | n/a | yes | Positional argument `planet_name`. |
| in | `initial_time` | InitialTime | n/a | yes | Positional argument `initial_time`. |
| in | `truth` | AtmosphereTruthConfig | n/a | yes | Positional argument `truth`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Base.invokelatest | n/a | — | Return value of `_make_required_gram_density_model`. Returns `Base.invokelatest(` or `offline_model`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__scenario_density_model|_scenario_density_model]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:191-191`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[analysis.scenario_builders__is_gram_library_missing_error|_is_gram_library_missing_error]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:339-339`
- `callees` → [[analysis.scenario_builders__try_libraryless_gram_surrogate|_try_libraryless_gram_surrogate]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:340-340`
<!-- vulcan:connections:end -->

## Limitations
The Mars-specific keywords are always forwarded regardless of planet, relying on the constructor to ignore them; the fallback warning is the only record that a run used the surrogate, and the returned configuration carries no flag.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 316.
