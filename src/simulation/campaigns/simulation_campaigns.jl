module SimulationCampaigns

using ..ParallelProfiles
using ..ParallelProfiles: OuterRouteFeatures, OuterRouteTuning, OuterRouteState
using ..ParallelProfiles: select_outer_route!, record_outer_route_feedback!
using ..ParallelProcess: ProcessPool, campaign_process_pool, ensure_process_workers!
using Distributed: remotecall_fetch

include(joinpath(@__DIR__, "monte_carlo.jl"))
include(joinpath(@__DIR__, "adaptive_routing.jl"))
include(joinpath(@__DIR__, "constellation_ensemble.jl"))

export MonteCarloSpec
export MonteCarloSampleResult
export MonteCarloResult
export run_monte_carlo
export run_constellation_ensemble
export campaign_route_features
export campaign_outer_route_state

end # module SimulationCampaigns
