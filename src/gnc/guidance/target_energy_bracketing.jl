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
    second_switch_reevaluation::Bool
    heat_load_security_mode::Bool
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
    switch_recompute_interval_s::Real=2.0,
    second_switch_reevaluation::Bool=true,
    heat_load_security_mode::Bool=true,
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
        second_switch_reevaluation,
        heat_load_security_mode,
    )
end

mutable struct AerobrakingEnergyDepletionState
    selected_mode::Vector{Symbol}
    targeting_active::Vector{Bool}
    safe_low_drag::Vector{Bool}
    energy_bracketing_evaluated::Vector{Bool}
    energy_bracketing_count::Vector{Int}
    target_energy_jkg::Vector{Float64}
    target_perturbation_energy_change_jkg::Vector{Float64}
    bracket_min_energy_jkg::Vector{Float64}
    bracket_max_energy_jkg::Vector{Float64}
    heat_load_switches_s::Vector{NTuple{2, Float64}}
    heat_load_switch_solved::Vector{Bool}
    heat_load_drag_passage_active::Vector{Bool}
    heat_load_entry_time_s::Vector{Float64}
    heat_load_last_reevaluation_s::Vector{Float64}
    heat_load_security_active::Vector{Bool}
    heat_load_previous_j_cm2::Vector{Float64}
    targeting_switch_s::Vector{Float64}
    last_switch_solve_t::Vector{Float64}
    last_alpha_rad::Vector{Float64}
    last_alpha_heat_rate_rad::Vector{Float64}
    last_alpha_structural_rad::Vector{Float64}
    last_heat_rate_w_cm2::Vector{Float64}
    last_heat_load_j_cm2::Vector{Float64}
    last_dynamic_pressure_pa::Vector{Float64}
    last_structural_load_pa::Vector{Float64}
    last_minimum_heat_rate_w_cm2::Vector{Float64}
    last_minimum_structural_load_pa::Vector{Float64}
    last_heat_rate_status::Vector{Symbol}
    last_structural_load_status::Vector{Symbol}
    last_heat_budget_status::Vector{Symbol}
    heat_load_entry_j_cm2::Vector{Vector{Float64}}
    heat_load_exit_j_cm2::Vector{Vector{Float64}}
    last_pass_heat_load_j_cm2::Vector{Float64}
    control_clock_origin_s::Float64
    targeting_diagnostics::Vector{NamedTuple}
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
        fill(NaN, n),
        fill((Inf, Inf), n),
        falses(n),
        falses(n),
        fill(NaN, n),
        fill(-Inf, n),
        falses(n),
        zeros(n),
        fill(Inf, n),
        fill(-Inf, n),
        fill(NaN, n),
        fill(NaN, n),
        fill(NaN, n),
        fill(NaN, n),
        fill(NaN, n),
        fill(NaN, n),
        fill(NaN, n),
        fill(NaN, n),
        fill(NaN, n),
        fill(:unobserved, n),
        fill(:unobserved, n),
        fill(:unobserved, n),
        [Float64[] for _ in 1:n],
        [Float64[] for _ in 1:n],
        fill(NaN, n),
        0.0,
        NamedTuple[(status=:unplanned,) for _ in 1:n],
    )
end

struct AerobrakingEnergyDepletionGuidanceModel <: AbstractGuidanceModel
    config::AerobrakingEnergyDepletionConfig
    state::AerobrakingEnergyDepletionState
end

@inline function _edg_state_index_ok(state::AerobrakingEnergyDepletionState, i::Int)::Bool
    return 1 <= i <= length(state.selected_mode)
end

@inline function _edg_sat_state(u, i::Int)
    return hasproperty(u, :sc) ? u.sc[i] : u
end

@inline function _edg_pos_vel_mass(sc)
    pos = hasproperty(sc, :pos) ? SVector{3, Float64}(sc.pos) : SVector{3, Float64}(sc[1], sc[2], sc[3])
    vel = hasproperty(sc, :vel) ? SVector{3, Float64}(sc.vel) : SVector{3, Float64}(sc[4], sc[5], sc[6])
    mass = hasproperty(sc, :mass) ? Float64(sc.mass) : (length(sc) >= 7 ? Float64(sc[7]) : NaN)
    return pos, vel, mass
end

function calcGuidanceEffect!(
    model::AerobrakingEnergyDepletionGuidanceModel,
    u,
    p::ODEParams,
    t::Float64,
    i::Int64,
)
    state = model.state
    _edg_state_index_ok(state, i) || return nothing
    config = model.config

    if :targeting in config.guidance_modes
        _edg_run_target_energy_bracketing!(model, u, p, Float64(t), i)
    elseif :max_energy_depletion in config.guidance_modes
        state.selected_mode[i] = :max_energy_depletion
        state.targeting_active[i] = false
        state.safe_low_drag[i] = false
        state.energy_bracketing_evaluated[i] = false
    else
        state.selected_mode[i] = :safe_low_drag
        state.targeting_active[i] = false
        state.safe_low_drag[i] = true
    end
    return nothing
end

function _edg_set_targeting_fallback!(
    config::AerobrakingEnergyDepletionConfig,
    state::AerobrakingEnergyDepletionState,
    i::Int,
)
    state.targeting_active[i] = false
    state.safe_low_drag[i] = !(:max_energy_depletion in config.guidance_modes)
    state.selected_mode[i] = (:max_energy_depletion in config.guidance_modes) ? :max_energy_depletion : :safe_low_drag
    return nothing
end

function _edg_run_target_energy_bracketing!(
    model::AerobrakingEnergyDepletionGuidanceModel,
    u,
    p::ODEParams,
    t::Float64,
    i::Int,
)
    config = model.config
    state = model.state
    ctrl = _control_module()
    if state.energy_bracketing_evaluated[i]
        # Journal VI.D: replan only before the applied switch, at 0.5 Hz by default.
        state.targeting_active[i] || return nothing
        ctrl._edg_switch_due(t, state.targeting_switch_s[i]) && return nothing
        ctrl._edg_replan_due(t, state.last_switch_solve_t[i], config.switch_recompute_interval_s) || return nothing
    end
    env = ctrl._edg_environment_state(u, p, t, i)
    ctrl._edg_in_drag_passage(p, env) || return nothing
    sc = _edg_sat_state(u, i)
    pos, vel, mass = _edg_pos_vel_mass(sc)
    heat_load = ctrl._edg_pass_heat_load_for_links(sc, config.controlled_panel_links, state, i)
    plan = ctrl._edg_journal_targeting_plan(config, p, p.args.dynamics_model.spacecraft[i],
        pos, vel, mass, t; heat_load_j_cm2=heat_load)
    state.energy_bracketing_evaluated[i] = true
    state.energy_bracketing_count[i] += 1
    state.target_energy_jkg[i] = plan.target_energy_jkg
    state.target_perturbation_energy_change_jkg[i] = plan.coast_energy_change_jkg
    state.bracket_min_energy_jkg[i] = plan.lower_energy_jkg
    state.bracket_max_energy_jkg[i] = plan.upper_energy_jkg
    state.targeting_active[i] = plan.mode == :targeting
    state.safe_low_drag[i] = plan.mode == :safe_low_drag
    state.selected_mode[i] = plan.mode
    state.targeting_switch_s[i] = plan.switch_s
    state.last_switch_solve_t[i] = t
    state.targeting_diagnostics[i] = plan
    return nothing
end
