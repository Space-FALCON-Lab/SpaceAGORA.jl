using ..EDGAlgorithms: _edg_pass_heat_load_for_links
# Passage heat is a difference of cumulative integrals, one baseline per panel.
# Do not reset the physical integrator state or change any guidance policy here.
_edg_heat_state(::Any) = nothing
_edg_heat_state(model::AerobrakingEnergyDepletionControlModel) = model.state
_edg_heat_state(model::AerobrakingEnergyDepletionGuidanceModel) = model.state

function _edg_heat_states(args)
    states = AerobrakingEnergyDepletionState[]
    seen = IdDict{AerobrakingEnergyDepletionState, Nothing}()
    for model in (args.control_model.control_effectors..., args.guidance_model.guidance_effectors...)
        state = _edg_heat_state(model)
        state === nothing && continue
        haskey(seen, state) && continue
        seen[state] = nothing
        push!(states, state)
    end
    return states
end

function _edg_initialize_heat_accounting!(args, u, t_start::Float64)
    states = _edg_heat_states(args)
    # This branch's checkpoint contains u only, not the guidance/pass history.
    # Never silently call a resumed cumulative integral a fresh passage budget.
    if !isempty(states) && t_start != 0.0
        throw(ArgumentError("EDG checkpoint resume requires passage heat baselines and guidance history; this checkpoint format stores neither."))
    end
    for state in states, i in eachindex(state.heat_load_entry_j_cm2)
        sc = _edg_control_sat_state(u, i)
        state.heat_load_entry_j_cm2[i] = hasproperty(sc, :heat_loads) ? zeros(length(sc.heat_loads)) : Float64[]
        state.heat_load_exit_j_cm2[i] = Float64[]
        state.last_pass_heat_load_j_cm2[i] = NaN
        state.heat_boundary_inside[i] = nothing
        state.heat_boundary_entry_time_s[i] = NaN
        state.heat_boundary_exit_time_s[i] = NaN
        state.heat_boundary_observed_at_s[i] = NaN
    end
    return nothing
end

# Use the same time, frame transform and geodetic altitude as
# _edg_environment_state, without sampling density or invoking guidance.
function _edg_heat_boundary_distance(u, p::ODEParams, t::Float64, i::Int)
    sc = _edg_control_sat_state(u, i)
    pos, vel, _ = _edg_control_pos_vel_mass(sc)
    env = p.args.environment_model
    pos_pp, _ = r_intor_p!(pos, vel, env.planet, _edg_ephemeris_time(p, t), env.ephemerides_model)
    return Float64(rtolatlong(pos_pp, env.planet)[1]) - 1e3 * Float64(env.EI)
end

# Crossing direction follows the same ellipsoid and time-dependent frame as
# the heat boundary. Radial inertial velocity can have the opposite sign at
# mid-latitudes and must not decide an omitted simultaneous geodetic event.
function _edg_heat_boundary_vertical_speed(u, p::ODEParams, t::Float64, i::Int)
    sc = _edg_control_sat_state(u, i)
    pos, vel, _ = _edg_control_pos_vel_mass(sc)
    env = p.args.environment_model
    pos_pp, vel_pp = r_intor_p!(pos, vel, env.planet, _edg_ephemeris_time(p, t), env.ephemerides_model)
    down, _, _ = latlongtoNED(rtolatlong(pos_pp, env.planet))
    vertical_speed = -dot(down, vel_pp)
    # Rotation and transport cancellation leave machine-level residue for a
    # stationary body-fixed state. Preserve the tangent-inside convention;
    # this tolerance scales only with floating-point transform roundoff.
    roundoff = 64eps(Float64) * max(1.0, norm(vel), norm(vel_pp))
    return abs(vertical_speed) <= roundoff ? 0.0 : vertical_speed
end

function _edg_capture_entry_heat!(args, u, i::Int; deduplicate::Bool=false, crossing_time_s::Float64=NaN, observed_at_s::Float64=NaN)
    sc = _edg_control_sat_state(u, i)
    for state in _edg_heat_states(args)
        _edg_control_state_index_ok(state, i) || continue
        deduplicate && state.heat_boundary_inside[i] === true && continue
        state.heat_boundary_inside[i] = true
        state.heat_boundary_entry_time_s[i] = crossing_time_s
        state.heat_boundary_exit_time_s[i] = NaN
        state.heat_boundary_observed_at_s[i] = observed_at_s
        GuidanceHooks._edg_reset_pass!(state, i)
        state.heat_load_entry_j_cm2[i] = hasproperty(sc, :heat_loads) ? Float64.(sc.heat_loads) : Float64[]
        state.heat_load_exit_j_cm2[i] = Float64[]

    end
    return nothing
end

function _edg_capture_exit_heat!(args, u, i::Int; deduplicate::Bool=false, crossing_time_s::Float64=NaN, observed_at_s::Float64=NaN)
    sc = _edg_control_sat_state(u, i)
    for state in _edg_heat_states(args)
        _edg_control_state_index_ok(state, i) || continue
        deduplicate && state.heat_boundary_inside[i] === false && continue
        state.heat_boundary_inside[i] = false
        state.heat_boundary_exit_time_s[i] = crossing_time_s
        state.heat_boundary_observed_at_s[i] = observed_at_s
        GuidanceHooks._edg_reset_pass!(state, i)
        state.heat_load_exit_j_cm2[i] = hasproperty(sc, :heat_loads) ? Float64.(sc.heat_loads) : Float64[]
    end
    return nothing
end


# Reconcile undelivered simultaneous members geometrically. An initially unknown
# member records its side without inventing an entry or resetting a t0 plan.
function _edg_reconcile_heat_boundary!(args, u, i::Int, inside::Bool; event::Bool, crossing_time_s::Float64=NaN, observed_at_s::Float64=NaN)
    states = _edg_heat_states(args)
    for state in states
        if !event && state.heat_boundary_inside[i] === nothing
            state.heat_boundary_inside[i] = inside
        end
    end
    if any(state -> state.heat_boundary_inside[i] !== inside, states)
        inside ? _edg_capture_entry_heat!(args, u, i; deduplicate=true, crossing_time_s, observed_at_s) :
            _edg_capture_exit_heat!(args, u, i; deduplicate=true, crossing_time_s, observed_at_s)
    end
    return nothing
end
