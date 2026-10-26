module ExampleForceTorqueModelTemplate

using StaticArrays
using SpaceAGORA

struct ExampleForceTorqueModel <: SpaceAGORA.AbstractForceTorqueModel
    force_n::SVector{3, Float64}
    torque_n_m::SVector{3, Float64}
end

# The force in the inertial frame (N) and the torque in the bus frame (N·m) on
# the spacecraft `x` describes; `x.spacecraft_index` is its position in the run.
function SpaceAGORA.wrench(
    model::ExampleForceTorqueModel,
    x::SpaceAGORA.StateSample,
    env::SpaceAGORA.EnvironmentSample,
    t::Float64,
)
    return model.force_n, model.torque_n_m
end

end # module ExampleForceTorqueModelTemplate
