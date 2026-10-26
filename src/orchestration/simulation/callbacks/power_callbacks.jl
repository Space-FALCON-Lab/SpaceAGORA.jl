# The power budget during a run (see PowerSystems): the snapshot that the
# battery derivative, the battery events and the result columns share; the
# events that mark the battery full or empty and the step-end update of its
# shedding; and the power result columns.

const _NO_SUN = SVector{3, Float64}(NaN, NaN, NaN)

# Generation, load and net power of spacecraft `sat_idx` in state `sc_view` at
# time `t`: the arrays lit by the Sun of the run's ephemerides (none when they
# cannot place it), the power beams it receives (from the run state `state`,
# which has every spacecraft), the constant loads, and the mounted actuators at
# the levels of every control effector in `:mounted` actuation that holds them,
# as the faults leave them.
function _power_snapshot(sc_view, spacecraft, p, sat_idx::Int, t::Float64; state=nothing)
    args = p.args
    planet = args.environment_model.planet
    generation = 0.0
    if any(link -> !isempty(link.solar_arrays), spacecraft.links)
        et = p.shared_buffers.et_start[] + t
        sun_dir = ephemerides_sun_ii(args.environment_model.ephemerides_model, planet, et)
        if sun_dir !== nothing
            pos = SVector{3, Float64}(sc_view.pos[1], sc_view.pos[2], sc_view.pos[3])
            sun_ii = SVector{3, Float64}(sun_dir) * mean_solar_distance_m(planet)
            lit = Float64(eclipse_area_calc(pos, sun_ii, Float64(planet.Rp_e)))
            to_sun = sun_ii - pos
            distance = norm(to_sun)
            q = hasproperty(sc_view, :q) ?
                SVector{4, Float64}(sc_view.q[1], sc_view.q[2], sc_view.q[3], sc_view.q[4]) :
                SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
            generation = solar_generation_w(spacecraft, rot(q) * (to_sun / distance), solar_flux_w_m2(distance), lit)
        end
    end
    if state !== nothing && link_receives(args.interlink_model, sat_idx)
        generation += received_link_power_w(args.interlink_model, sat_idx, args.dynamics_model.spacecraft, state;
            context=LinkContext(time_s=t, planet=planet), planet_frame=_station_planet_frame(p, t))
    end
    wheel_momentum = hasproperty(sc_view, :h_wheels) ? sc_view.h_wheels : nothing
    effectors = args.control_model.control_effectors
    balance = _mounted_power_balance(spacecraft, effectors, sat_idx, generation, wheel_momentum;
        interlinks=args.interlink_model, members=_scope_members(p, :control, effectors, sat_idx))
    return (generation_w = generation, load_w = balance.load_w, unshed_load_w = balance.unshed_load_w,
            by_tier = balance.by_tier, net_w = balance.net_w)
end

# The power balance (as `PowerSystems.power_balance` gives it for one set of
# levels) with the mounted actuators of every control effector in `:mounted`
# actuation that holds levels for spacecraft `sat_idx`, among `members` (the
# effectors acting on it, `_scope_members`), at the levels its faults leave them:
# the dynamics apply each one's actuators, so each one's draw is a load. Terminal
# hardware an interlink drives draws what its active link sets instead
# (`interlinks`, the run's InterlinkModel).
function _mounted_power_balance(spacecraft, effectors, sat_idx::Int, generation_w::Float64, wheel_momentum;
                                interlinks=nothing, members=eachindex(effectors))
    n = mounted_actuator_count(spacecraft)
    bus = bus_load_w(spacecraft)
    unshed = (0.0, 0.0, 0.0, 0.0)
    shed = 0.0
    if n > 0
        if link_driven(interlinks, sat_idx)
            by_link, link_load = link_terminal_power(interlinks, sat_idx, spacecraft)
            unshed = unshed .+ by_link
            shed += link_load
        end
        for j in members
            effector = effectors[j]
            actuation_mode(effector) === :mounted || continue
            levels = actuator_levels(effector, sat_idx)
            (levels === nothing || length(levels) != n) && continue
            levels = faulted_levels(spacecraft, without_link_driven(levels, interlinks, sat_idx))
            by_tier = actuator_power_by_tier(spacecraft, levels; wheel_momentum=wheel_momentum)
            applied = shed_levels(spacecraft, levels)
            unshed = unshed .+ by_tier
            shed += sum(applied === levels ? by_tier : actuator_power_by_tier(spacecraft, applied; wheel_momentum=wheel_momentum))
        end
    end
    load = bus + shed
    return (load_w = load, unshed_load_w = bus + sum(unshed), by_tier = unshed, net_w = generation_w - load)
end

# Spacecraft with a battery, whose stored energy is a state.
_battery_spacecraft(args) = [i for i in eachindex(args.dynamics_model.spacecraft) if has_battery(args.dynamics_model.spacecraft[i])]

# A discrete or event callback counts as a derivative discontinuity unless it
# clears the flag: declare one when the power state changed (the solver then
# evaluates the derivative again) and clear it otherwise, so the step-end
# update costs no re-evaluation when nothing switched.
@inline function _mark_discontinuity!(integrator, changed::Bool)
    if applicable(DiffEqBase.derivative_discontinuity!, integrator, changed)
        DiffEqBase.derivative_discontinuity!(integrator, changed)
    end
    return nothing
end

# Step-end update of every battery's flags and shedding tier.
function _update_power_states!(integrator, sats::Vector{Int})
    spacecraft = integrator.p.args.dynamics_model.spacecraft
    changed = false
    for i in sats
        sc_view = integrator.u.sc[i]
        snap = _power_snapshot(sc_view, spacecraft[i], integrator.p, i, integrator.t; state=integrator.u)
        changed |= update_power_state!(spacecraft[i], sc_view.battery_wh, snap.generation_w, snap.by_tier)
    end
    _mark_discontinuity!(integrator, changed)
    return changed
end

"""
    get_power_callbacks(num_sats, args) -> Tuple

For runs with a battery: an event on each battery reaching zero or its
capacity, which marks it empty or full, and a step-end update of its flags and
shedding tier (see `PowerSystems.update_power_state!`). Both mark the
derivative discontinuous when they change anything, so no integration step
sees the battery's state change. Empty for runs without a battery.

The flags carry over from one checkpoint segment to the next, and a run
resumed from a checkpoint file takes them from the file (which holds the
stored energy in the state and the flags with the callbacks' run state).
"""
function get_power_callbacks(num_sats::Int, args::SimulationConfiguration)
    sats = _battery_spacecraft(args)
    isempty(sats) && return ()
    spacecraft = args.dynamics_model.spacecraft
    function condition!(out, u, t, integrator)
        @inbounds for (k, i) in enumerate(sats)
            energy = u.sc[i].battery_wh
            out[2k - 1] = energy
            out[2k] = battery(spacecraft[i]).capacity_wh - energy
        end
        return nothing
    end
    # `events` holds one crossing direction per condition (0: no crossing).
    function boundary!(integrator, events)
        hit = false
        for e in eachindex(events)
            events[e] == 0 && continue
            i = sats[(e + 1) ÷ 2]
            hit |= battery_boundary!(integrator.p.args.dynamics_model.spacecraft[i], integrator.u.sc[i].battery_wh)
        end
        hit || (_mark_discontinuity!(integrator, false); return nothing)
        _update_power_states!(integrator, sats)
        _mark_discontinuity!(integrator, true)
        return nothing
    end
    boundary = VectorContinuousCallback(condition!, boundary!, 2 * length(sats); save_positions=(false, false))   # both crossing directions
    # A fresh callback is built for every run. Its initialization also runs at
    # checkpoint-segment boundaries, which keep the flags: they depend on the
    # run so far (an empty battery recharging below `recovery_soc` is still
    # empty), not only on the stored energy; a resumed run's come from its
    # checkpoint.
    initialized = Ref(false)
    function initialize!(cb, u, t, integrator)
        if !initialized[] && !integrator.p.shared_buffers.callback_resume[]
            for i in sats
                sc = integrator.p.args.dynamics_model.spacecraft[i]
                b = battery(sc)
                energy = u.sc[i].battery_wh
                b.full = within_full_tolerance(b, energy)
                b.empty = energy <= 1.0e-9 * b.capacity_wh
                b.shed_tier = 0
            end
            initialized[] = true
        end
        _update_power_states!(integrator, sats)
        _mark_discontinuity!(integrator, true)
        return nothing
    end
    step_end = DiscreteCallback((u, t, integrator) -> true, integrator -> _update_power_states!(integrator, sats);
        initialize=initialize!, save_positions=(false, false))
    return (boundary, step_end)
end

function _save_power(num_sats::Int, u, t, integrator, which::Symbol)
    spacecraft = integrator.p.args.dynamics_model.spacecraft
    out = Vector{Float64}(undef, num_sats)
    @inbounds for i in 1:num_sats
        sc = spacecraft[i]
        if !has_power_parts(sc)
            out[i] = NaN
        elseif which === :power_received_w
            out[i] = received_link_power_w(integrator.p.args.interlink_model, i, spacecraft, u;
                context=LinkContext(time_s=Float64(t), planet=integrator.p.args.environment_model.planet),
                planet_frame=_station_planet_frame(integrator.p, t))
        elseif which === :power_beamed_w
            out[i] = beamed_link_power_w(integrator.p.args.interlink_model, i)
        elseif which === :battery_wh || which === :battery_soc
            b = battery(sc)
            out[i] = b === nothing ? NaN :
                (which === :battery_wh ? Float64(u.sc[i].battery_wh) : Float64(u.sc[i].battery_wh) / b.capacity_wh)
        else
            snap = _power_snapshot(u.sc[i], sc, integrator.p, i, Float64(t); state=u)
            out[i] = which === :power_generated_w ? snap.generation_w :
                     which === :power_load_w ? snap.load_w :
                     which === :power_net_w ? snap.net_w : snap.unshed_load_w - snap.load_w
        end
    end
    return out
end

"""
    power_save_fields(args) -> Vector{SaveField}

The power columns, when a spacecraft carries power-system parts: per
spacecraft `battery_wh` and `battery_soc` (NaN without a battery),
`power_generated_w`, `power_load_w` (with the shedding applied),
`power_shed_w` and `power_net_w` (generation minus load; negative while the
battery cannot cover the load). NaN for a spacecraft without power parts. A
run with power beams adds `power_received_w`, the electrical power its incoming
beams deliver (part of `power_generated_w`), and `power_beamed_w`, the optical
power its outgoing beams carry.
"""
function power_save_fields(args::SimulationConfiguration)
    any(has_power_parts, args.dynamics_model.spacecraft) || return SaveField[]
    num_sats = length(args.dynamics_model.spacecraft)
    names = (:battery_wh, :battery_soc, :power_generated_w, :power_load_w, :power_shed_w, :power_net_w)
    model = hasproperty(args, :interlink_model) ? args.interlink_model : nothing
    if model !== nothing && any(link -> link isa PowerBeamLink, model.graph.edge_properties)
        names = (names..., :power_received_w, :power_beamed_w)
    end
    return SaveField[SaveField(name, (u, t, integrator) -> _save_power(num_sats, u, t, integrator, name); per_satellite=true)
                     for name in names]
end
