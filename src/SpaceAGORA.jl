__precompile__(true)

module SpaceAGORA

using PrecompileTools: @compile_workload, @setup_workload

include(joinpath(@__DIR__, "parallel", "routing", "parallel_profiles.jl"))
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
using .SimulationEngine.SimulationModel.AbstractTypes: AbstractForceTorqueModel, AbstractDensityModel
using .SimulationEngine.SimulationModel.AbstractTypes: AbstractControlEffectorModel, AbstractEphemeridesModel
using .SimulationEngine.SimulationModel.AbstractTypes: AbstractThermalModel, AbstractThrusterModel, AbstractGuidanceModel
using .TelemetryVerification: VerificationRequest, VerificationResult
using .TelemetryVerification: run_verification, run_verification_cli, run_study
using .SpaceAGORACLI: AssetCheckItem, AssetCheckReport

@doc (@doc SimulationEngine.ParallelConfig) ParallelConfig
@doc (@doc SimulationEngine.SolverConfig) SolverConfig
@doc (@doc SimulationEngine.RuntimePolicyConfig) RuntimePolicyConfig
@doc (@doc SimulationEngine.ArtifactConfig) ArtifactConfig
@doc (@doc SimulationEngine.SimulationEngineConfig) SimulationEngineConfig
@doc (@doc SimulationEngine.simulation_engine_config_from_env) simulation_engine_config_from_env

@doc (@doc SimulationEngine.SimulationModel.AbstractTypes.AbstractForceTorqueModel) AbstractForceTorqueModel
@doc (@doc SimulationEngine.SimulationModel.AbstractTypes.AbstractDensityModel) AbstractDensityModel
@doc (@doc SimulationEngine.SimulationModel.AbstractTypes.AbstractControlEffectorModel) AbstractControlEffectorModel
@doc (@doc SimulationEngine.SimulationModel.AbstractTypes.AbstractEphemeridesModel) AbstractEphemeridesModel
@doc (@doc SimulationEngine.SimulationModel.AbstractTypes.AbstractThermalModel) AbstractThermalModel
@doc (@doc SimulationEngine.SimulationModel.AbstractTypes.AbstractThrusterModel) AbstractThrusterModel
@doc (@doc SimulationEngine.SimulationModel.AbstractTypes.AbstractGuidanceModel) AbstractGuidanceModel

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
export AbstractForceTorqueModel, AbstractDensityModel, AbstractControlEffectorModel
export AbstractEphemeridesModel, AbstractThermalModel, AbstractThrusterModel, AbstractGuidanceModel
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
