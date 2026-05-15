Base.@kwdef mutable struct RPOGuidanceModel <: AbstractGuidanceModel
    chaser_idx::Int = 1
    target_idx::Int = 2
    goal_rtn::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 0.0)
    geometry::Any = nothing
    plan_buffer::RPOPlanBuffer = RPOPlanBuffer()
    pso_config::Any = nothing
    safe_distance_m::Float64 = 0.0
    force_replan::Bool = false
end

