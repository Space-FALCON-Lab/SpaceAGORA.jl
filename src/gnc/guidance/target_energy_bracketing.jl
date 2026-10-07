using ..EDGServices

const _EDG_GUIDANCE_MODES = Set((:max_energy_depletion, :targeting))
const _EDG_MAX_ENERGY_SUBMODES = Set((:heat_rate, :structural_load, :heat_load))
const _EDG_HEAT_LOAD_SWITCH_SOLVERS = Set((:closed_form, :tpbvp_integration))

@inline function _edg_symbol_tuple(values)
    values isa Symbol && return (values,)
    return tuple((Symbol(v) for v in values)...)
end

@inline function _edg_validate_symbol_set(values::Tuple, allowed::Set{Symbol}, label::String)
    isempty(values) && throw(ArgumentError("$(label) must not be empty."))
    for value in values
        value in allowed || throw(ArgumentError("Unsupported $(label) value $(value). Allowed values are $(sort!(collect(allowed)))."))
    end
    return values
end

"""
    AerobrakingEnergyDepletionConfig

Configuration for the energy-depletion aerobraking guidance strategy: the
available guidance modes and max-energy submodes, the heat-load switch
solver, the solar-panel links under angle-of-attack control, the target
apoapsis, the commanded angle-of-attack bounds, the heat-rate / heat-load /
structural-load limits, and the planning-horizon, switch-recompute, and
targeting-certification settings.
"""
struct AerobrakingEnergyDepletionConfig
    guidance_modes::Tuple{Vararg{Symbol}}
    max_energy_submodes::Tuple{Vararg{Symbol}}
    heat_load_switch_solver::Symbol
    controlled_panel_links::Tuple{Vararg{Int}}
    target_apoapsis_radius_m::Float64
    max_alpha_rad::Float64
    min_alpha_rad::Float64
    heat_rate_limit_w_cm2::Float64
    heat_load_limit_j_cm2::Float64
    structural_load_limit_pa::Float64
    planning_horizon_s::Float64
    switch_recompute_interval_s::Float64
    targeting_certification_samples::Int
    targeting_energy_order_tolerance_jkg::Float64
    targeting_heat_load_tolerance_j_cm2::Float64
end

function AerobrakingEnergyDepletionConfig(;
    guidance_modes=(:max_energy_depletion,),
    max_energy_submodes=(:heat_rate, :structural_load, :heat_load),
    heat_load_switch_solver::Symbol=:closed_form,
    controlled_panel_links=(2, 3),
    target_apoapsis_radius_m::Real=NaN,
    max_alpha_rad::Real=pi / 2,
    min_alpha_rad::Real=1e-4,
    heat_rate_limit_w_cm2::Real=Inf,
    heat_load_limit_j_cm2::Real=Inf,
    structural_load_limit_pa::Real=Inf,
    planning_horizon_s::Real=5_000.0,
    switch_recompute_interval_s::Real=30.0,
    targeting_certification_samples::Integer=9,
    targeting_energy_order_tolerance_jkg::Real=1e-3,
    targeting_heat_load_tolerance_j_cm2::Real=1e-6,
)
    guidance_modes_t = _edg_validate_symbol_set(_edg_symbol_tuple(guidance_modes), _EDG_GUIDANCE_MODES, "guidance_modes")
    max_energy_submodes_t = _edg_validate_symbol_set(_edg_symbol_tuple(max_energy_submodes), _EDG_MAX_ENERGY_SUBMODES, "max_energy_submodes")
    heat_load_switch_solver in _EDG_HEAT_LOAD_SWITCH_SOLVERS ||
        throw(ArgumentError("Unsupported heat_load_switch_solver $(heat_load_switch_solver). Use :closed_form or :tpbvp_integration."))
    panel_links = tuple((Int(idx) for idx in controlled_panel_links)...)
    isempty(panel_links) && throw(ArgumentError("controlled_panel_links must contain at least one link index."))
    any(<=(0), panel_links) && throw(ArgumentError("controlled_panel_links must be positive 1-based link indices."))

    max_alpha = Float64(max_alpha_rad)
    min_alpha = Float64(min_alpha_rad)
    isfinite(max_alpha) && isfinite(min_alpha) && 0.0 <= min_alpha <= max_alpha ||
        throw(ArgumentError("Expected finite alpha bounds with 0 <= min_alpha_rad <= max_alpha_rad."))
    horizon = Float64(planning_horizon_s)
    recompute = Float64(switch_recompute_interval_s)
    isfinite(horizon) && horizon > 0.0 || throw(ArgumentError("planning_horizon_s must be finite and > 0.0."))
    isfinite(recompute) && recompute > 0.0 || throw(ArgumentError("switch_recompute_interval_s must be finite and > 0.0."))
    certification_samples = Int(targeting_certification_samples)
    certification_samples >= 2 || throw(ArgumentError("targeting_certification_samples must be >= 2."))
    energy_order_tolerance = Float64(targeting_energy_order_tolerance_jkg)
    heat_load_tolerance = Float64(targeting_heat_load_tolerance_j_cm2)
    isfinite(energy_order_tolerance) && energy_order_tolerance >= 0.0 ||
        throw(ArgumentError("targeting_energy_order_tolerance_jkg must be finite and >= 0.0."))
    isfinite(heat_load_tolerance) && heat_load_tolerance >= 0.0 ||
        throw(ArgumentError("targeting_heat_load_tolerance_j_cm2 must be finite and >= 0.0."))

    return AerobrakingEnergyDepletionConfig(
        guidance_modes_t,
        max_energy_submodes_t,
        heat_load_switch_solver,
        panel_links,
        Float64(target_apoapsis_radius_m),
        max_alpha,
        min_alpha,
        Float64(heat_rate_limit_w_cm2),
        Float64(heat_load_limit_j_cm2),
        Float64(structural_load_limit_pa),
        horizon,
        recompute,
        certification_samples,
        energy_order_tolerance,
        heat_load_tolerance,
    )
end

"""
    AerobrakingEnergyDepletionState

Per-spacecraft mutable state for energy-depletion guidance: the selected
mode, targeting/bracketing activity flags and counters, and the target and
bracket specific-energy values (J/kg) maintained across guidance calls.
"""
mutable struct AerobrakingEnergyDepletionState
    selected_mode::Vector{Symbol}
    targeting_active::Vector{Bool}
    safe_low_drag::Vector{Bool}
    energy_bracketing_evaluated::Vector{Bool}
    energy_bracketing_count::Vector{Int}
    target_energy_jkg::Vector{Float64}
    bracket_min_energy_jkg::Vector{Float64}
    bracket_max_energy_jkg::Vector{Float64}
    heat_load_switches_s::Vector{NTuple{2, Float64}}
    heat_load_switch_solved::Vector{Bool}
    heat_load_drag_passage_active::Vector{Bool}
    targeting_switch_s::Vector{Float64}
    last_switch_solve_t::Vector{Float64}
    last_alpha_rad::Vector{Float64}
    last_alpha_heat_rate_rad::Vector{Float64}
    last_alpha_structural_rad::Vector{Float64}
    last_heat_rate_w_cm2::Vector{Float64}
    last_heat_load_j_cm2::Vector{Float64}
    last_dynamic_pressure_pa::Vector{Float64}
end

function AerobrakingEnergyDepletionState(; num_sats::Integer)
    n = Int(num_sats)
    n > 0 || throw(ArgumentError("num_sats must be positive."))
    return AerobrakingEnergyDepletionState(
        fill(:inactive, n),
        falses(n),
        falses(n),
        falses(n),
        zeros(Int, n),
        fill(NaN, n),
        fill(NaN, n),
        fill(NaN, n),
        fill((Inf, Inf), n),
        falses(n),
        falses(n),
        fill(Inf, n),
        fill(-Inf, n),
        fill(NaN, n),
        fill(NaN, n),
        fill(NaN, n),
        fill(NaN, n),
        fill(NaN, n),
        fill(NaN, n),
    )
end

"""
    AerobrakingEnergyDepletionGuidanceModel

Guidance model implementing the energy-depletion aerobraking strategy:
brackets the spacecraft's target specific energy each pass and selects the
guidance mode subject to the configured heat and structural limits. Holds an
[`AerobrakingEnergyDepletionConfig`](@ref) and its mutable
[`AerobrakingEnergyDepletionState`](@ref).
"""
struct AerobrakingEnergyDepletionGuidanceModel <: AbstractGuidanceModel
    config::AerobrakingEnergyDepletionConfig
    state::AerobrakingEnergyDepletionState
end

@inline function _edg_state_index_ok(state::AerobrakingEnergyDepletionState, i::Int)::Bool
    return 1 <= i <= length(state.selected_mode)
end

const _edg_sat_state = EDGServices._edg_control_sat_state

const _edg_pos_vel_mass = EDGServices._edg_control_pos_vel_mass

function calcGuidanceEffect!(
    model::AerobrakingEnergyDepletionGuidanceModel,
    u,
    p::ODEParams,
    t::Float64,
    i::Int64,
)
    return _edg_algorithms().guidance_decision!(model.config, model.state, u, p, t, i)
end

@inline _edg_algorithms() = getfield(_PARENT, :EDGAlgorithms)

function _edg_interpolate_bracket_value(exit_energy::Float64, energy_min::Float64, energy_max::Float64, value_at_min::Float64, value_at_max::Float64)
    return _edg_algorithms()._edg_interpolate_bracket_value(exit_energy, energy_min, energy_max, value_at_min, value_at_max)
end

function _edg_target_energy_from_reachable_bracket(
    planet,
    target_apoapsis_radius_m::Float64,
    energy_min::Float64,
    energy_max::Float64,
    periapsis_at_min::Float64,
    periapsis_at_max::Float64,
)
    return _edg_algorithms()._edg_target_energy_from_reachable_bracket(planet, target_apoapsis_radius_m, energy_min, energy_max, periapsis_at_min, periapsis_at_max)
end

function _edg_set_targeting_fallback!(
    config::AerobrakingEnergyDepletionConfig,
    state::AerobrakingEnergyDepletionState,
    i::Int,
)
    return _edg_algorithms()._edg_set_targeting_fallback!(config, state, i)
end

function _edg_run_target_energy_bracketing!(
    model::AerobrakingEnergyDepletionGuidanceModel,
    u,
    p::ODEParams,
    t::Float64,
    i::Int,
)
    return _edg_algorithms()._edg_run_target_energy_bracketing!(model.config, model.state, u, p, t, i)
end
