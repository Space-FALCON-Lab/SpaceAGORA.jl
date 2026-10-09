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

function _edg_capture_entry_heat!(args, u, i::Int)
    sc = _edg_control_sat_state(u, i)
    for state in _edg_heat_states(args)
        _edg_control_state_index_ok(state, i) || continue
        state.heat_load_entry_j_cm2[i] = hasproperty(sc, :heat_loads) ? Float64.(sc.heat_loads) : Float64[]
        state.heat_load_exit_j_cm2[i] = Float64[]
    end
    return nothing
end

function _edg_capture_exit_heat!(args, u, i::Int)
    sc = _edg_control_sat_state(u, i)
    for state in _edg_heat_states(args)
        _edg_control_state_index_ok(state, i) || continue
        state.heat_load_exit_j_cm2[i] = hasproperty(sc, :heat_loads) ? Float64.(sc.heat_loads) : Float64[]
    end
    return nothing
end

function _edg_pass_heat_load_for_links(sc, links::Tuple{Vararg{Int}}, state, i::Int)::Float64
    hasproperty(sc, :heat_loads) || return 0.0
    baseline = state.heat_load_entry_j_cm2[i]
    loads = isempty(state.heat_load_exit_j_cm2[i]) ? sc.heat_loads : state.heat_load_exit_j_cm2[i]
    value = 0.0
    for idx in links
        if 1 <= idx <= length(loads)
            entry = idx <= length(baseline) ? baseline[idx] : 0.0
            load = Float64(loads[idx]) - entry
            isfinite(load) && (value = max(value, max(0.0, load)))
        end
    end
    return value
end
