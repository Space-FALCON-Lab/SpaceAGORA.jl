"""Guidance model state for RPO planning, tracking, and optional replanning."""
Base.@kwdef mutable struct RPOGuidanceModel <: AbstractGuidanceModel
    chaser_idx::Int = 1
    target_idx::Int = 2
    goal_rtn::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 0.0)
    geometry::Any = nothing
    plan_buffer::RPOPlanBuffer = RPOPlanBuffer()
    pso_config::Any = nothing
    safe_distance_m::Float64 = 0.0
    force_replan::Bool = false
    replanning_config::Any = nothing
    replanning_events::Vector{NamedTuple} = NamedTuple[]
    replan_count::Int = 0
    retime_count::Int = 0
    safe_hold_count::Int = 0
    replan_failure_count::Int = 0
    last_replanning_time_s::Float64 = -Inf
    last_replanning_signature::UInt = UInt(0)
    replanning_persistence_count::Int = 0
end
