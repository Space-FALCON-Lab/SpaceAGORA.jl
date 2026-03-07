__precompile__(false)

module SpaceAGORA

include(joinpath(@__DIR__, "parallel", "routing", "parallel_profiles.jl"))
include(joinpath(@__DIR__, "simulation", "engine", "simulation_engine.jl"))
include(joinpath(@__DIR__, "analysis", "verification", "telemetry_verification.jl"))

using .ParallelProfiles: ParallelProfile, ParallelProfileConfig
using .ParallelProfiles: parse_parallel_profile, parallel_profile_name, profile_config, profile_env_pairs, with_parallel_profile
using .ParallelProfiles: OuterRouteFeatures, OuterRouteTuning, OuterRouteState
using .ParallelProfiles: reset_outer_route_state!, outer_route_signature, outer_route_stats_snapshot
using .ParallelProfiles: default_outer_route, outer_route_candidates, select_outer_route!, record_outer_route_feedback!
using .SimulationEngine: ParallelConfig, SolverConfig, RuntimePolicyConfig, ArtifactConfig, SimulationEngineConfig
using .SimulationEngine: simulation_engine_config_from_env
using .TelemetryVerification: VerificationRequest, VerificationResult
using .TelemetryVerification: run_verification, run_verification_cli, run_study

export ParallelProfile, ParallelProfileConfig
export parse_parallel_profile, parallel_profile_name, profile_config, profile_env_pairs, with_parallel_profile
export OuterRouteFeatures, OuterRouteTuning, OuterRouteState
export reset_outer_route_state!, outer_route_signature, outer_route_stats_snapshot
export default_outer_route, outer_route_candidates, select_outer_route!, record_outer_route_feedback!
export ParallelConfig, SolverConfig, RuntimePolicyConfig, ArtifactConfig, SimulationEngineConfig
export simulation_engine_config_from_env
export VerificationRequest, VerificationResult
export run_verification, run_verification_cli, run_study, run_simulation

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
"""
run_simulation(args...; kwargs...) = SimulationEngine.run_simulation(args...; kwargs...)
run_simulation(config::SimulationEngineConfig, args...; kwargs...) = SimulationEngine.run_simulation(config, args...; kwargs...)

end # module SpaceAGORA
