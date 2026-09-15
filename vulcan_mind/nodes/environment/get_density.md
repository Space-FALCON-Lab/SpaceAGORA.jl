---
id: environment.get_density
label: getDensity
kind: function
source:
  file: src/environment/atmosphere/density_models.jl
  symbol: getDensity
  lines:
  - 714
  - 813
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: EnvironmentModels namespace registering atmosphere-model density methods.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: atmosphere
  type: Tuple{Float64,Float64,SVector{3,Float64}}
  units: kg/m^3,K,m/s
  description: Density, temperature, and wind vector produced by a concrete atmosphere
    model dispatch.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
- atmosphere
charts:
- environment
origin: agent
---

# getDensity

## Purpose
`getDensity` is the atmosphere-model dispatch surface used by aerodynamic effectors and callbacks. The density-model file defines concrete methods for vacuum, exponential, tabulated, polynomial, NRLMSISE, constant-density, and optional GRAM-backed models. Each method returns the environmental quantities needed to compute aerodynamic force and thermal or telemetry diagnostics.

## Theory & Math
The basic exponential model follows `ρ(h)=ρ_ref exp((h_ref-h)/H)`, with density `ρ`, altitude `h`, reference density and altitude, and scale height `H`. Other methods interpolate tables or call external atmosphere models. The return tuple is `(ρ,T,w)`, where `T` is temperature and `w` is the wind vector.

## Model & Assumptions
Inputs use the package’s altitude, latitude, longitude, epoch, and wind-request conventions. Model validity ranges are advisory for some analytic models, while tabulated models depend on coverage and interpolation rules. Aerodynamic callers assume density units are kg/m³ and wind is expressed in the expected planet-relative frame.

## Design & Implementation
The methods near lines 805 onward implement no-atmosphere, exponential, and piecewise-exponential behavior; later methods cover tabulated and external models. `physical_models.jl` exports `getDensity`, and `aerodynamic_wrench_models.jl` calls it while building a wrench. The GRAM extension adds `EM.getDensity` methods when the optional dependency is loaded.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | EnvironmentModels namespace registering atmosphere-model density methods. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `atmosphere` | Tuple{Float64,Float64,SVector{3,Float64}} | kg/m^3,K,m/s | — | Density, temperature, and wind vector produced by a concrete atmosphere model dispatch. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models__aero_link_atmosphere_query|_aero_link_atmosphere_query]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:509-509`
- [[environment.density_models__density_scalar_for_batch|_density_scalar_for_batch]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:873-873`
- [[environment.density_models__gram_point_density|_gram_point_density]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:1058-1058`
- [[environment.density_models_density_polyfit|density_polyfit]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:1083-1083`
- [[environment.density_models_getdensitybatch_bang|getDensityBatch!]] · `callees` → `callers` · call · `src/environment/atmosphere/density_models.jl:1019-1019`
- [[environment.density_models_timetabulatedatmospheremodel|TimeTabulatedAtmosphereModel]] · `callees` → `callers` · feedback · `src/environment/atmosphere/density_models.jl:784-784`
- [[ext.spaceagoragramsuiteext__gram_spice_ephemeris_state|_gram_spice_ephemeris_state]] · `callees` → `callers` · call · `ext/SpaceAGORAGRAMSuiteExt.jl:271-271`
- [[gnc.heat_load_control__edg_sample_prediction_atmosphere|_edg_sample_prediction_atmosphere]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:112-112`
- [[gnc.targeting_control__edg_environment_state|_edg_environment_state]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:72-72`
- [[gnc.targeting_control__edg_targeting_prediction_environment|_edg_targeting_prediction_environment]] · `callees` → `callers` · call · `src/gnc/control/targeting_control.jl:335-335`
- [[grp.src_simulation_callbacks|simulation/callbacks/]] · `members_out` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:201-201`
- [[simulation.runtime__density_state_from_kinematics_bang|_density_state_from_kinematics!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/runtime.jl:201-201`
- [[simulation.vacuum_predicted_gram__query_vacuum_gram_cache__query_vacuum_gram_cache_bang|_query_vacuum_gram_cache!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:273-273`
- [[simulation_a.refresh_gram_track_cache_refresh__gram_track_cache_refresh_bang|_gram_track_cache_refresh!]] · `callees` → `callers` · call · `src/simulation/callbacks/gram_track_cache/refresh.jl:291-291`
- [[simulation_a.vacuum_predicted_gram_build_vacuum_gram_cache__build_vacuum_gram_cache_bang|_build_vacuum_gram_cache!]] · `callees` → `callers` · call · `src/simulation/callbacks/density_callbacks/vacuum_predicted_gram.jl:218-218`
- [[spaceagora.spaceagora_spaceagora|SpaceAGORA]] · `callees` → `callers` · call · `src/SpaceAGORA.jl:354-354`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/environment/atmosphere/density_models.jl:773-773`
- `callees` → [[environment.density_models__exponential_density|_exponential_density]] · `callers` · call · `src/environment/atmosphere/density_models.jl:813-813`
- `callees` → [[environment.density_models__tab_flight_interp|_tab_flight_interp]] · `callers` · call · `src/environment/atmosphere/density_models.jl:727-727`
- `callees` → [[environment.density_models_timetabulatedatmospheremodel|TimeTabulatedAtmosphereModel]] · `callers` · call · `src/environment/atmosphere/density_models.jl:755-755`
<!-- vulcan:connections:end -->

## Limitations
Analytic models extrapolate outside their documented altitude bands, while tables can fail or become inaccurate outside coverage. External GRAM evaluation depends on native data and shared locking. The function does not decide whether density is negligible for a particular force model; that thresholding belongs to aerodynamic callers.

## Provenance
Mapped from `src/environment/atmosphere/density_models.jl:805-813` and the model dispatches included by `src/environment/physical_models.jl`.
