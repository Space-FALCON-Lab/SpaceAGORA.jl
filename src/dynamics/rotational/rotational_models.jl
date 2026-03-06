module DynamicsRotational

using LinearAlgebra
using StaticArrays

include(joinpath(@__DIR__, "attitude_kinematics.jl"))
include(joinpath(@__DIR__, "rigid_body_dynamics.jl"))
include(joinpath(@__DIR__, "torque_models.jl"))

export quaternion_derivative
export angular_acceleration
export body_torque
export body_angular_velocity
export combine_torques

end # module DynamicsRotational
