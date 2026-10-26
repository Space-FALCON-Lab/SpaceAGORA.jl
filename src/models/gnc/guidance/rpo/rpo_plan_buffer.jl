"""Retimed RPO trajectory reference and diagnostics produced by a planner."""
Base.@kwdef mutable struct RPOPlan
    valid::Bool = false
    t_ref_s::Vector{Float64} = Float64[]
    r_ref_rtn::Matrix{Float64} = zeros(3, 0)
    v_ref_rtn::Matrix{Float64} = zeros(3, 0)
    path_rtn::Matrix{Float64} = zeros(3, 0)
    cost::Float64 = Inf
    diagnostics::NamedTuple = NamedTuple()
end

"""Mutable buffer that stores the active and previous RPO plans across guidance updates."""
Base.@kwdef mutable struct RPOPlanBuffer
    valid::Bool = false
    plan::RPOPlan = RPOPlan()
    updated_at_s::Float64 = NaN
end

"""Install a new RPO plan while preserving the previous plan for diagnostics."""
function update_rpo_plan_buffer!(buffer::RPOPlanBuffer, plan::RPOPlan, t::Real)
    buffer.plan = plan
    buffer.valid = plan.valid
    buffer.updated_at_s = Float64(t)
    return buffer
end

