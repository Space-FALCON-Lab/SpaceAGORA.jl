# Canonical aggregator: no behavior ownership.
module ParallelProfiles

using Dates
using TOML

include(joinpath(@__DIR__, "machine_topology.jl"))
include(joinpath(@__DIR__, "profile_definitions.jl"))
include(joinpath(@__DIR__, "env_mapping.jl"))
include(joinpath(@__DIR__, "outer_route_state.jl"))
include(joinpath(@__DIR__, "outer_route_selection.jl"))
include(joinpath(@__DIR__, "outer_route_metrics.jl"))

export machine_topology, refresh_machine_topology!, usable_core_budget
export physical_core_count, allowed_cpus, cgroup_cpu_quota
export ParallelProfile, ParallelProfileConfig
export parse_parallel_profile, parallel_profile_name, profile_config, profile_env_pairs, with_parallel_profile
export OuterRouteFeatures, OuterRouteTuning, OuterRouteState
export reset_outer_route_state!, outer_route_signature, outer_route_stats_snapshot
export default_outer_route, outer_route_candidates, select_outer_route!, record_outer_route_feedback!
export load_outer_route_state!, save_outer_route_state
export outer_split_candidates, select_outer_split!, record_outer_split_feedback!

end # module ParallelProfiles
