"""
Kustaanheimo-Stiefel (KS) regularized orbital dynamics: state conversion, the
nonlinear right-hand side under J2 and drag, an RK4 propagator, and the
analytic Jacobians and discrete tangent maps used by the aerobraking MPC.

The state is `[u (4), u′ (4), h_KS, t]`, with `u′ = du/ds` in the Sundman
fictitious time `s` (`dt/ds = |r|`) and `h_KS = -ε` the negative specific
orbital energy.
"""
module DynamicsKS

using LinearAlgebra
using StaticArrays

export KSPropagationParams
export ks_energy_parameter, specific_energy_from_ks
export ks_position, ks_velocity, cartesian_to_ks_state, ks_state_to_cartesian
export ks_rotation_cross_matrix, ks_relative_velocity
export ks_j2_acceleration_si, ks_drag_acceleration_si, ks_density_value
export nonlinear_ks_rhs!, nonlinear_ks_rhs
export ks_rk4_step
export ks_kinematics_jacobians, ks_j2_acceleration_jacobian_si
export ks_density_value_gradient, ks_density_difference
export evaluate_ks_dynamics_and_jacobians
export ks_linear_implicit_midpoint_step, ks_first_order_tangent_map

include(joinpath(@__DIR__, "transform.jl"))
include(joinpath(@__DIR__, "nonlinear_dynamics.jl"))
include(joinpath(@__DIR__, "linearized_dynamics.jl"))

end # module DynamicsKS
