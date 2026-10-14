"""Opt-in contract adapter for existing HYPR planning and configured retiming."""
module HYPRRPOPlanning
using LinearAlgebra: norm
using ..RPOPlannerInterfaces
using ..SimulationModel
const P = RPOPlannerInterfaces
const S = SimulationModel
const G = S.GuidanceHooks
import ..RPOPlannerInterfaces: planner_capabilities, plan_rpo!

"""
Wrap a copied HYPR configuration. The physical request owns clearance and the
reference interval. Enabled request limits cap the configured retiming limits
using prospective headroom. Other settings and optimizer-returned adaptations
are retained. `rrt_on_replan` explicitly enables the existing warm start for
`:replan` requests. Retiming/restart are not advertised by this first adapter.
"""
struct HYPRRPOPlanner <: P.AbstractRPOPlanner
    config::S.RPOPSOConfig
    headroom::P.RPOPlanningHeadroom
    max_reference_samples::Int
    rrt_on_replan::Bool
    function HYPRRPOPlanner(config::S.RPOPSOConfig; headroom=P.RPOPlanningHeadroom(),
                            max_reference_samples::Integer=100_000, rrt_on_replan::Bool=false)
        2 <= max_reference_samples < typemax(Int) || throw(ArgumentError("Invalid reference sample budget."))
        new(deepcopy(S.rpo_pso_config(config)), headroom, Int(max_reference_samples), rrt_on_replan)
    end
end
planner_capabilities(::HYPRRPOPlanner) = P.RPOPlannerCapabilities(state_sources=(:truth,), frames=(:target_rtn,))

function _plan_hypr_rpo! end
function plan_rpo!(state::Nothing, planner::HYPRRPOPlanner, request::P.RPOPlanningRequest, rng::P.AbstractRNG)
    S.HYPRSupport.require_hypr()
    return _plan_hypr_rpo!(state, planner, request, rng)
end
P.require_planner_support(::HYPRRPOPlanner) = S.HYPRSupport.require_hypr()
P.initialize_planner(::HYPRRPOPlanner, context) = (S.HYPRSupport.require_hypr(); nothing)
import ..SimulationLifecycle
function SimulationLifecycle.preflight_guidance(g::S.RPOGuidanceModel, args; isolate_state)
    replanning = G._rpo_replanning_config(g)
    if !g.plan_buffer.valid || g.force_replan || (replanning !== nothing && replanning.enabled)
        S.HYPRSupport.require_hypr()
    end
    return nothing
end
end
