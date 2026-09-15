---
id: environment.density_models_timetabulatedatmospheremodel
label: TimeTabulatedAtmosphereModel
kind: struct
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: TimeTabulatedAtmosphereModel
  lines:
  - 748
  - 748
inputs:
- id: t_el_s
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `t_el_s`.
- id: log_rho
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `log_rho`.
- id: scale
  type: Float64
  units: n/a
  required: true
  description: Field `scale`.
- id: temperature_k
  type: Float64
  units: n/a
  required: true
  description: Field `temperature_k`.
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
  type: TimeTabulatedAtmosphereModel
  units: n/a
  description: Constructed `TimeTabulatedAtmosphereModel`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- environment
origin: agent
---

# TimeTabulatedAtmosphereModel

## Purpose
Density as a pure function of scenario elapsed time, for replaying an along-track density history or an assimilated product sampled along the trajectory.

## Design & Implementation
Immutable with sorted elapsed times, log-density at each node, a multiplicative `scale` and a constant `temperature_k`. The constructor requires equal lengths, at least two nodes, sorted finite times, positive finite densities, and positive finite scale and temperature. `getDensity` holds the end values beyond the table and interpolates log-linearly inside.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `t_el_s` | Vector{Float64} | n/a | yes | Field `t_el_s`. |
| in | `log_rho` | Vector{Float64} | n/a | yes | Field `log_rho`. |
| in | `scale` | Float64 | n/a | yes | Field `scale`. |
| in | `temperature_k` | Float64 | n/a | yes | Field `temperature_k`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | TimeTabulatedAtmosphereModel | n/a | — | Constructed `TimeTabulatedAtmosphereModel`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__make_time_tabulated_density_model|_make_time_tabulated_density_model]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:265-265`
- [[environment.get_density|getDensity]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:755-755`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/environment/atmosphere/density_models.jl:773-773`
- `callees` → [[environment.density_models__exponential_density|_exponential_density]] · `callers` · call · `src/environment/atmosphere/density_models.jl:813-813`
- `callees` → [[environment.density_models__piecewise_layer_index|_piecewise_layer_index]] · `callers` · call · `src/environment/atmosphere/density_models.jl:824-824`
- `callees` → [[environment.get_density|getDensity]] · `callers` · feedback · `src/environment/atmosphere/density_models.jl:784-784`
<!-- vulcan:connections:end -->

## Limitations
It ignores altitude entirely, so it is only meaningful when the replayed trajectory matches the one the table was derived from; the table epoch must equal the scenario epoch and nothing here can verify that.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl` line 748.
