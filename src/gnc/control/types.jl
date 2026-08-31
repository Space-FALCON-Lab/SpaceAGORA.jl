#=
"""
    Main data types for the aerobraking MPC.

    The controller configuration does not carry mission limits or tuning
    defaults. Those values must be supplied by the example or user script so
    there is no hidden behind-the-scenes control case.
"""
=#
abstract type AerobrakingMPCMode end

struct TargetEnergyMode <: AerobrakingMPCMode end
struct MaxEnergyDepletionMode <: AerobrakingMPCMode end

Base.@kwdef struct AerobrakingMPCParams
    Re::Float64
    μ::Float64
    J2::Float64
    Ω::Float64
end

Base.@kwdef struct AerobrakingMPCConfig{M <: AerobrakingMPCMode}
    mode::M
    bus_reference_area_m2::Float64
    controllable_area_m2::Float64
    mass_kg::Float64
    drag_coefficient::Float64
    qdot_max_w_cm2::Float64
    heat_load_max_j_cm2::Float64
    drag_max_n::Float64
    area_slew_max_m2_s::Float64
    use_constraints::Bool
    use_slew_constraint::Bool
    use_qdot_constraint::Bool
    use_heat_load_constraint::Bool
    use_drag_constraint::Bool
    target_energy_mj_kg::Float64
    area_weight::Float64
    area_slew_weight::Float64
    slack_weight::Float64
    target_energy_weight::Float64
    max_depletion_energy_weight::Float64
    osqp_eps_abs::Float64
    osqp_eps_rel::Float64
    osqp_max_iter::Int
end

Base.@kwdef struct AerobrakingMPCProblem
    params::AerobrakingMPCParams
    H::Matrix{Float64}
    Mx::Matrix{Float64}
    δX0::Vector{Float64}
    N::Int
    ny::Int
    t::Vector{Float64}
    Ybar::Matrix{Float64}
    Xbar::Matrix{Float64}
    Abar_m2::Float64
end

Base.@kwdef struct AerobrakingMPCSolution
    mode::Symbol
    ok::Bool
    ΔU::Vector{Float64}
    slacks::Vector{Float64}
    objective::Float64
    terminal_energy::Float64
    commanded_area_m2::Vector{Float64}
    predicted_outputs::Matrix{Float64}
    solver_status::Symbol = :not_solved
end

Base.@kwdef mutable struct AerobrakingMPCState
    selected_mode::Symbol = :inactive
    last_solve_time_s::Float64 = -Inf
    solver_status::Symbol = :not_solved
    held_commanded_area_m2::Float64 = NaN
    held_alpha_rad::Float64 = NaN
    last_command_time_s::Float64 = NaN
    estimated_heat_load_j_cm2::Float64 = 0.0
    last_heat_rate_w_cm2::Float64 = NaN
    last_drag_n::Float64 = NaN
    last_density_kg_m3::Float64 = NaN
    active_limit::Symbol = :none
    plan_epoch_s::Float64 = NaN
    plan_time_s::Vector{Float64} = Float64[]
    plan_area_m2::Vector{Float64} = Float64[]
    predicted_terminal_energy::Float64 = NaN
    last_solution::Union{Nothing, AerobrakingMPCSolution} = nothing
end
