# Core aerobraking MPC modes, configuration, problem, solution, and runtime state.
"""
Objective selected by an aerobraking MPC configuration: [`TargetEnergyMode`](@ref)
or [`MaxEnergyDepletionMode`](@ref).
"""
abstract type AerobrakingMPCMode end

"""
Drive the predicted end-of-pass specific energy to `target_energy_mj_kg`.

The QP minimises `target_energy_weight * (1e6 * (E_N - E*))^2` with `E_N` in
MJ/kg, plus the area and slew penalties; heat-rate, heat-load and drag limits
are soft (penalised slacks), as in the energy-targeting benchmark formulation.
"""
struct TargetEnergyMode <: AerobrakingMPCMode end
"""
Maximise the energy removed during one atmospheric pass.

The QP minimises the predicted end-of-pass specific energy (J/kg), scaled by
`max_depletion_energy_weight / (max|E| * controllable_area_m2)`, plus the area
and slew penalties. The heat-rate and heat-load limits are hard and the drag
limit is soft, as in the maximum-energy-depletion benchmark formulation.
"""
struct MaxEnergyDepletionMode <: AerobrakingMPCMode end

"""
Settings of one aerobraking MPC (keyword constructor; every field is required).

- `mode`: [`TargetEnergyMode`](@ref) or [`MaxEnergyDepletionMode`](@ref).
- `bus_reference_area_m2`, `controllable_area_m2`: the exposed area ranges over
  `[bus, bus + controllable]` m²; panels give `bus + controllable * sin(α)`.
- `mass_kg`, `drag_coefficient`: the drag model `0.5 ρ V² C_D A`.
- `qdot_max_w_cm2`, `heat_load_max_j_cm2`, `drag_max_n`, `area_slew_max_m2_s`:
  path and actuator limits; `use_constraints` and the `use_*_constraint` flags
  select the active ones (see [`mpc_constraints`](@ref)).
- `target_energy_mj_kg`, `target_energy_weight`: energy-targeting objective.
- `max_depletion_energy_weight`: maximum-depletion objective weight.
- `area_weight`, `area_slew_weight`, `slack_weight`: quadratic penalties on the
  area increment, its node-to-node change and the constraint slacks.
- `osqp_eps_abs`, `osqp_eps_rel`, `osqp_max_iter`: OSQP settings.

Build one from a configuration with [`mpc_config_from_spaceagora`](@ref).
"""
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

"""
Condensed linear prediction of one atmospheric pass: `δY = H ΔU + Mx δX₀`
about the nominal outputs `Ybar` (`N × 4`, columns altitude m, drag N, heat rate
W/m² and specific energy, MJ/kg for [`TargetEnergyMode`](@ref) and J/kg for
[`MaxEnergyDepletionMode`](@ref)), at node times `t` (s) of the nominal KS
states `Xbar`, linearised about the exposed area `Abar_m2`. Built by
[`build_mpc_problem`](@ref) and solved by [`solve_mpc_qp`](@ref).
"""
Base.@kwdef struct AerobrakingMPCProblem
    params::KSPropagationParams
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

"""
Result of [`solve_mpc_qp`](@ref): the area increments `ΔU` (m²), the commanded
area plan `commanded_area_m2`, the constraint `slacks`, the OSQP `objective` and
`solver_status`, the predicted outputs (as [`AerobrakingMPCProblem`](@ref)
`Ybar`) and their `terminal_energy`. `ok` is `false` when OSQP did not solve the
problem; then `ΔU` is zero.
"""
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

"""
Run state of an [`AerobrakingMPCControlModel`](@ref): the held area and panel
angle, the current plan (relative times from `plan_epoch_s`), the estimated
heat load of the current pass, the last solution, solver status and error, and
the limit that last constrained the command (`active_limit`).
"""
Base.@kwdef mutable struct AerobrakingMPCState
    last_solve_time_s::Float64 = -Inf
    solver_status::Symbol = :not_solved
    held_commanded_area_m2::Float64 = NaN
    held_alpha_rad::Float64 = NaN
    last_command_time_s::Float64 = NaN
    estimated_heat_load_j_cm2::Float64 = 0.0
    last_heat_rate_w_cm2::Float64 = NaN
    last_drag_n::Float64 = NaN
    active_limit::Symbol = :none
    plan_epoch_s::Float64 = NaN
    plan_time_s::Vector{Float64} = Float64[]
    plan_area_m2::Vector{Float64} = Float64[]
    predicted_terminal_energy::Float64 = NaN
    solve_armed::Bool = true
    solve_count::Int = 0
    last_solution::Union{Nothing, AerobrakingMPCSolution} = nothing
    last_error::Union{Nothing, String} = nothing
end
