---
id: analysis.scenario_builders__try_libraryless_gram_surrogate
label: _try_libraryless_gram_surrogate
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _try_libraryless_gram_surrogate
  lines:
  - 298
  - 298
inputs:
- id: planet_name
  type: String
  units: n/a
  required: true
  description: Positional argument `planet_name`.
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
  type: GRAMAtmosphereModelSurrogate
  units: n/a
  description: Return value of `_try_libraryless_gram_surrogate`. Returns `GRAMAtmosphereModelSurrogate(base_model,
    surrogate_file, nothing)`.
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

# _try_libraryless_gram_surrogate

## Purpose
Attempts to construct the offline GRAM surrogate for a planet when the native library cannot be loaded.

## Design & Implementation
Returns `nothing` if the scenario pins `gram_offline_surrogate` to `off` or the environment switch is disabled. It resolves the default surrogate file through `Base.invokelatest` on the extension-installed resolver, swallowing any error into an empty string, and returns `nothing` if the file is missing. Otherwise it wraps a `_GRAMOfflineSurrogateFallbackBase` in `GRAMAtmosphereModelSurrogate` with no point-fallback altitude.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet_name` | String | n/a | yes | Positional argument `planet_name`. |
| in | `truth` | AtmosphereTruthConfig | n/a | yes | Positional argument `truth`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | GRAMAtmosphereModelSurrogate | n/a | — | Return value of `_try_libraryless_gram_surrogate`. Returns `GRAMAtmosphereModelSurrogate(base_model, surrogate_file, nothing)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__make_required_gram_density_model|_make_required_gram_density_model]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:340-340`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[analysis.scenario_builders__gramofflinesurrogatefallbackbase|_GRAMOfflineSurrogateFallbackBase]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:312-312`
- `callees` → [[analysis.scenario_builders__libraryless_gram_surrogate_enabled|_libraryless_gram_surrogate_enabled]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:303-303`
- `callees` → [[environment.density_models_gramatmospheremodelsurrogate|GRAMAtmosphereModelSurrogate]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:313-313`
<!-- vulcan:connections:end -->

## Limitations
`invokelatest` is needed because the resolver slot is filled at extension load time; if the extension never loaded the call throws and is swallowed, so the function cannot distinguish no file from no extension.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 298.
