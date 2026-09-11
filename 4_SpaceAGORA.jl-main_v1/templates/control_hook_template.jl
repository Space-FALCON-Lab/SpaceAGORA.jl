module ExampleControlHookTemplate

using StaticArrays
using ComponentArrays
using SpaceAGORA

struct ExampleControlHook <: SpaceAGORA.AbstractControlEffectorModel
    force_n::SVector{3, Float64}
    torque_n_m::SVector{3, Float64}
    mass_flow_kg_s::Float64
end

function SpaceAGORA.calcControlEffect!(
    model::ExampleControlHook,
    u::ComponentVector,
    p,
    t::Float64,
    i::Int64,
)
    return nothing
end

function SpaceAGORA.calcControlForceTorque(
    model::ExampleControlHook,
    u::AbstractVector,
    p,
    i::Int64,
    t::Float64,
)
    return model.force_n, model.torque_n_m
end

function SpaceAGORA.calcControlMassFlowRate(
    model::ExampleControlHook,
    u::AbstractVector,
    p,
    i::Int64,
    t::Float64,
)
    return model.mass_flow_kg_s
end

end # module ExampleControlHookTemplate
