module DynamicsKS

using LinearAlgebra
using StaticArrays

export KSPropagationParams
export ks_energy_parameter, specific_energy_from_ks
export ks_position, ks_velocity, cartesian_to_ks_state, ks_state_to_cartesian
export ks_rotation_cross_matrix
export ks_j2_acceleration_si, ks_drag_acceleration_si, ks_rhs!, ks_rhs
export nonlinear_ks_rhs!, nonlinear_ks_rhs
export ks_rk4_step
export ks_kinematics_jacobians, ks_j2_acceleration_jacobian_si
export ks_density_value_gradient
export ks_rhs_jacobians, ks_rhs_jacobian
export evaluate_ks_dynamics_and_jacobians
export ks_linear_implicit_midpoint_step, ks_first_order_tangent_map

include(joinpath(@__DIR__, "transform.jl"))
include(joinpath(@__DIR__, "nonlinear_dynamics.jl"))
include(joinpath(@__DIR__, "linearized_dynamics.jl"))

end # module DynamicsKS
