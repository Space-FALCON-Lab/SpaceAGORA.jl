__precompile__(false)

module SpaceAGORA

include(joinpath(@__DIR__, "parallel", "ParallelProfiles.jl"))
include(joinpath(@__DIR__, "analysis", "verification", "TelemetryVerification.jl"))

using .ParallelProfiles: ParallelProfile, ParallelProfileConfig
using .ParallelProfiles: parse_parallel_profile, parallel_profile_name, profile_config, profile_env_pairs, with_parallel_profile
using .ParallelProfiles: OuterRouteFeatures, OuterRouteTuning, OuterRouteState
using .ParallelProfiles: reset_outer_route_state!, outer_route_signature, outer_route_stats_snapshot
using .ParallelProfiles: default_outer_route, outer_route_candidates, select_outer_route!, record_outer_route_feedback!
using .TelemetryVerification: VerificationRequest, VerificationResult
using .TelemetryVerification: run_verification, run_verification_cli, run_study

export ParallelProfile, ParallelProfileConfig
export parse_parallel_profile, parallel_profile_name, profile_config, profile_env_pairs, with_parallel_profile
export OuterRouteFeatures, OuterRouteTuning, OuterRouteState
export reset_outer_route_state!, outer_route_signature, outer_route_stats_snapshot
export default_outer_route, outer_route_candidates, select_outer_route!, record_outer_route_feedback!
export VerificationRequest, VerificationResult
export run_verification, run_verification_cli, run_study, run_simulation

"""
    run_simulation(args...; kwargs...)

Stable package entrypoint for simulation execution used by calibration integrations.
"""
run_simulation(args...; kwargs...) = TelemetryVerification.run_simulation(args...; kwargs...)

end # module SpaceAGORA
