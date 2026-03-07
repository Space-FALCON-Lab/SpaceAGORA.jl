__precompile__(true)

module SpaceAGORA

using PrecompileTools: @compile_workload, @setup_workload

include(joinpath(@__DIR__, "parallel", "routing", "parallel_profiles.jl"))
include(joinpath(@__DIR__, "core", "simulation_model.jl"))
include(joinpath(@__DIR__, "simulation", "engine", "simulation_engine.jl"))
include(joinpath(@__DIR__, "analysis", "verification", "telemetry_verification.jl"))
include(joinpath(@__DIR__, "cli", "spaceagora_cli.jl"))

using .ParallelProfiles: ParallelProfile, ParallelProfileConfig
using .ParallelProfiles: parse_parallel_profile, parallel_profile_name, profile_config, profile_env_pairs, with_parallel_profile
using .ParallelProfiles: OuterRouteFeatures, OuterRouteTuning, OuterRouteState
using .ParallelProfiles: reset_outer_route_state!, outer_route_signature, outer_route_stats_snapshot
using .ParallelProfiles: default_outer_route, outer_route_candidates, select_outer_route!, record_outer_route_feedback!
using .SimulationEngine: ParallelConfig, SolverConfig, RuntimePolicyConfig, ArtifactConfig, SimulationEngineConfig
using .SimulationEngine: simulation_engine_config_from_env
using .SimulationModel.AbstractTypes: AbstractForceTorqueModel, AbstractPlanet, AbstractDensityModel
using .SimulationModel.AbstractTypes: AbstractControlEffectorModel, AbstractEphemeridesModel
using .SimulationModel.AbstractTypes: AbstractThermalModel, AbstractThrusterModel, AbstractGuidanceModel
using .SimulationModel: NoAtmosphereModel, ExponentialAtmosphereModel, SimpleEphemeridesModel
using .SimulationModel: make_no_gram_planet, make_no_gram_density_model, make_no_gram_environment
using .SimulationModel: calcForceTorque, getDensity, getDensityBatch!
using .SimulationModel: calcControlEffect!, calcControlForceTorque, calcControlMassFlowRate
using .TelemetryVerification: VerificationRequest, VerificationResult
using .TelemetryVerification: run_verification, run_verification_cli, run_study
using .SpaceAGORACLI: AssetCheckItem, AssetCheckReport

@doc (@doc SimulationEngine.ParallelConfig) ParallelConfig
@doc (@doc SimulationEngine.SolverConfig) SolverConfig
@doc (@doc SimulationEngine.RuntimePolicyConfig) RuntimePolicyConfig
@doc (@doc SimulationEngine.ArtifactConfig) ArtifactConfig
@doc (@doc SimulationEngine.SimulationEngineConfig) SimulationEngineConfig
@doc (@doc SimulationEngine.simulation_engine_config_from_env) simulation_engine_config_from_env

@doc (@doc SimulationModel.AbstractTypes.AbstractForceTorqueModel) AbstractForceTorqueModel
@doc (@doc SimulationModel.AbstractTypes.AbstractPlanet) AbstractPlanet
@doc (@doc SimulationModel.AbstractTypes.AbstractDensityModel) AbstractDensityModel
@doc (@doc SimulationModel.AbstractTypes.AbstractControlEffectorModel) AbstractControlEffectorModel
@doc (@doc SimulationModel.AbstractTypes.AbstractEphemeridesModel) AbstractEphemeridesModel
@doc (@doc SimulationModel.AbstractTypes.AbstractThermalModel) AbstractThermalModel
@doc (@doc SimulationModel.AbstractTypes.AbstractThrusterModel) AbstractThrusterModel
@doc (@doc SimulationModel.AbstractTypes.AbstractGuidanceModel) AbstractGuidanceModel

"""
    NoAtmosphereModel()

Density-model constructor for no-atmosphere baseline runs and no-GRAM onboarding
scenarios.
"""
NoAtmosphereModel

"""
    ExponentialAtmosphereModel(planet)
    ExponentialAtmosphereModel(rho_ref, h_ref, H)

Analytic density-model constructor for baseline runs that should not depend on
GRAM assets. The `planet` convenience form uses the built-in reference density
and scale-height constants for the chosen body.
"""
ExponentialAtmosphereModel

"""
    SimpleEphemeridesModel(; reference_epoch_seconds=0.0, prime_meridian_at_reference_rad=0.0)

Analytic ephemerides/frame backend for onboarding and open-data runs that should
not depend on local SPICE kernels.
"""
SimpleEphemeridesModel

@doc (@doc SimulationModel.make_no_gram_planet) make_no_gram_planet
@doc (@doc SimulationModel.make_no_gram_density_model) make_no_gram_density_model
@doc (@doc SimulationModel.make_no_gram_environment) make_no_gram_environment

"""
    calcForceTorque(model, x, p, i) -> (force_n, torque_n_m)

Stable extension hook for custom [`AbstractForceTorqueModel`](@ref)
implementations. Extend this method for package or user models that contribute
translational and rotational wrench terms to the simulation RHS.
"""
calcForceTorque

"""
    getDensity(model, h, lat, lon, el_time, wind[, p]) -> (rho, temperature, wind_vec)

Stable extension hook for custom [`AbstractDensityModel`](@ref)
implementations. The scalar form returns density, temperature, and wind for a
single atmosphere query.
"""
getDensity

"""
    getDensityBatch!(rhos, Ts, winds, model, hs, lats, lons, el_time, wind, p)

Optional batch extension hook for [`AbstractDensityModel`](@ref)
implementations that can answer many atmosphere queries more efficiently than
repeated scalar `getDensity` dispatch.
"""
getDensityBatch!

@doc (@doc SimulationModel.calcControlEffect!) calcControlEffect!
@doc (@doc SimulationModel.calcControlForceTorque) calcControlForceTorque
@doc (@doc SimulationModel.calcControlMassFlowRate) calcControlMassFlowRate

@doc (@doc ParallelProfiles.ParallelProfile) ParallelProfile
@doc (@doc ParallelProfiles.ParallelProfileConfig) ParallelProfileConfig
@doc (@doc ParallelProfiles.parse_parallel_profile) parse_parallel_profile
@doc (@doc ParallelProfiles.parallel_profile_name) parallel_profile_name
@doc (@doc ParallelProfiles.profile_config) profile_config
@doc (@doc ParallelProfiles.profile_env_pairs) profile_env_pairs
@doc (@doc ParallelProfiles.with_parallel_profile) with_parallel_profile
@doc (@doc ParallelProfiles.OuterRouteFeatures) OuterRouteFeatures
@doc (@doc ParallelProfiles.OuterRouteTuning) OuterRouteTuning
@doc (@doc ParallelProfiles.OuterRouteState) OuterRouteState
@doc (@doc ParallelProfiles.reset_outer_route_state!) reset_outer_route_state!
@doc (@doc ParallelProfiles.outer_route_signature) outer_route_signature
@doc (@doc ParallelProfiles.outer_route_stats_snapshot) outer_route_stats_snapshot
@doc (@doc ParallelProfiles.default_outer_route) default_outer_route
@doc (@doc ParallelProfiles.outer_route_candidates) outer_route_candidates
@doc (@doc ParallelProfiles.select_outer_route!) select_outer_route!
@doc (@doc ParallelProfiles.record_outer_route_feedback!) record_outer_route_feedback!

@doc (@doc TelemetryVerification.VerificationRequest) VerificationRequest
@doc (@doc TelemetryVerification.VerificationResult) VerificationResult
@doc (@doc TelemetryVerification.run_verification) run_verification
@doc (@doc TelemetryVerification.run_verification_cli) run_verification_cli
@doc (@doc TelemetryVerification.run_study) run_study

@doc (@doc SpaceAGORACLI.AssetCheckItem) AssetCheckItem
@doc (@doc SpaceAGORACLI.AssetCheckReport) AssetCheckReport

export ParallelProfile, ParallelProfileConfig
export parse_parallel_profile, parallel_profile_name, profile_config, profile_env_pairs, with_parallel_profile
export OuterRouteFeatures, OuterRouteTuning, OuterRouteState
export reset_outer_route_state!, outer_route_signature, outer_route_stats_snapshot
export default_outer_route, outer_route_candidates, select_outer_route!, record_outer_route_feedback!
export ParallelConfig, SolverConfig, RuntimePolicyConfig, ArtifactConfig, SimulationEngineConfig
export simulation_engine_config_from_env
export AbstractForceTorqueModel, AbstractPlanet, AbstractDensityModel, AbstractControlEffectorModel
export AbstractEphemeridesModel, AbstractThermalModel, AbstractThrusterModel, AbstractGuidanceModel
export NoAtmosphereModel, ExponentialAtmosphereModel, SimpleEphemeridesModel
export make_no_gram_planet, make_no_gram_density_model, make_no_gram_environment
export calcForceTorque, getDensity, getDensityBatch!
export calcControlEffect!, calcControlForceTorque, calcControlMassFlowRate
export VerificationRequest, VerificationResult
export run_verification, run_verification_cli, run_study, run_simulation
export AssetCheckItem, AssetCheckReport, check_assets, render_asset_report, run_cli

"""
    run_simulation(args...; isolate_state=true, kwargs...)

Stable package entrypoint for simulation execution used by calibration integrations.

By default, `isolate_state=true` deep-copies the simulation configuration before execution.
This preserves correctness and reentrancy across repeated runs and concurrent callers by
preventing one run from mutating shared campaign or model state that another run still
references.

Set `isolate_state=false` only as an advanced performance lever when the caller owns the
configuration instance and will not reuse it concurrently or across runs that may mutate
shared state. This can reduce setup cost for large mission definitions or many short runs,
but it trades away the default isolation guarantee.

# Examples
```jldoctest
julia> args = spaceagora_no_gram_example_args();

julia> sol = run_simulation(args; return_solution=true);

julia> length(sol.t) > 1
true
```
"""
run_simulation(args...; kwargs...) = SimulationEngine.run_simulation(args...; kwargs...)
run_simulation(config::SimulationEngineConfig, args...; kwargs...) = SimulationEngine.run_simulation(config, args...; kwargs...)

"""
    check_assets(; repo_root=pwd()) -> AssetCheckReport

Inspect the current repository asset layout and report which baseline, optional,
and high-fidelity asset roots are available.
"""
check_assets(args...; kwargs...) = SpaceAGORACLI.check_assets(args...; kwargs...)

"""
    render_asset_report(report; io=stdout)

Render a human-readable asset status report.
"""
render_asset_report(args...; kwargs...) = SpaceAGORACLI.render_asset_report(args...; kwargs...)

"""
    run_cli([args=ARGS]; io=stdout, errio=stderr) -> Int

Stable CLI entrypoint for SpaceAGORA operational commands:

- `run`
- `telemetry`
- `benchmark`
- `assets check`

This is the package-owned command surface used by the `bin/spaceagora` wrapper.
"""
run_cli(args...; kwargs...) = SpaceAGORACLI.run_cli(args...; kwargs...)

include(joinpath(@__DIR__, "precompile_workload.jl"))

end # module SpaceAGORA
