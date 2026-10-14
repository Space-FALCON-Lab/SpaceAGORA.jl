module ExampleForceTorqueModelTemplate

using StaticArrays
using SpaceAGORA

struct ExampleForceTorqueModel <: SpaceAGORA.AbstractForceTorqueModel
    force_n::SVector{3, Float64}
    torque_n_m::SVector{3, Float64}
end

function SpaceAGORA.calcForceTorque(
    model::ExampleForceTorqueModel,
    x::AbstractVector{Float64},
    p,
    i::Int64,
)
    return model.force_n, model.torque_n_m
end

end # module ExampleForceTorqueModelTemplate
