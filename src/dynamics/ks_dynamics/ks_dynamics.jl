module DynamicsKS

using LinearAlgebra
using StaticArrays

export KSPropagationParams
export ks_energy_parameter, specific_energy_from_ks
export ks_position, ks_velocity, cartesian_to_ks_state, ks_state_to_cartesian
export ks_j2_acceleration_si, ks_drag_acceleration_si, ks_rhs, ks_rk4_step
export ks_kinematics_jacobians, ks_j2_acceleration_jacobian_si
export ks_rhs_jacobian, ks_step_jacobian

include(joinpath(@__DIR__, "transform.jl"))
include(joinpath(@__DIR__, "dynamics.jl"))
include(joinpath(@__DIR__, "linearization.jl"))

end # module DynamicsKS
