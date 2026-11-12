# Externally propagated spacecraft ("shadow entries"): setup-time refusals, per-run runtimes and the shadow
# right-hand side. The hook surface is `SimulationModel.ExternalPropagation`; the sync callback that overwrites
# the shadow entries after every accepted step lives in `callbacks/external_propagation_callbacks.jl`.
#
# A shadow entry stays in `u.sc` so default save fields, guidance, navigation and `bind_spacecraft` see it like
# any other spacecraft. Between syncs it follows the owner's chief acceleration and plain attitude kinematics;
# at each sync it is overwritten from the owner's absolute state.

const _EXTERNAL_ALLOWED_SOLVER_MODES = (:tsit5, :auto_stiff, :rodas5p, :dp8)
# The spacecraft-index fields of a control effector.
const _EXTERNAL_INDEX_FIELDS = (:chaser_index, :spacecraft_index)

@inline _external_propagators(args) = args.external_propagators

function _external_capabilities(ep)
    c = SimulationModel.ExternalPropagation.external_capabilities(ep)
    c isa NamedTuple && keys(c) == (:configured_forces, :checkpoint) &&
        c.configured_forces isa Bool && c.checkpoint isa Bool ||
        throw(ArgumentError("external_capabilities must return (configured_forces::Bool, checkpoint::Bool)."))
    return c
end

_external_stage_safe(model) = model isa Union{SimulationModel.InverseSquaredGravityModel,
    SimulationModel.InverseSquaredJ2GravityModel} ? !model.gravity_gradient : SimulationModel.ExternalPropagation.external_stage_safe(model)

_external_stage_key(model) = model isa Union{SimulationModel.InverseSquaredGravityModel,
    SimulationModel.InverseSquaredJ2GravityModel} ? (kind=string(typeof(model)), fields=Tuple(getfield(model,n) for n in fieldnames(typeof(model)))) :
    SimulationModel.ExternalPropagation.external_stage_key(model)

function _external_force_configuration(args)
    planet = args.environment_model.planet
    fields = (:name, :μ, :Rp_e, :Rp_p, :Rp_m, :ω, :J2)
    ephemeris = args.environment_model.ephemerides_model
    return (planet=Tuple((name, deepcopy(getproperty(planet,name))) for name in fields if hasproperty(planet,name)),
        ephemeris=(kind=string(typeof(ephemeris)),fields=Tuple(getfield(ephemeris,n) for n in fieldnames(typeof(ephemeris)))),
        epoch=Tuple(getfield(args.initial_time,n) for n in fieldnames(typeof(args.initial_time))),
        forces=Tuple(_external_stage_key(m) for m in args.dynamics_model.dynamic_effectors),
        controls=Tuple(_external_stage_key(m) for m in args.control_model.control_effectors),
        spacecraft_mass=Tuple((dry_mass=sc.dry_mass,prop_mass=sc.prop_mass) for sc in args.dynamics_model.spacecraft),
        control_rates=Tuple(args.control_model.control_rates),
        guidance_rates=Tuple(args.guidance_model.guidance_rates), navigation_rates=Tuple(args.navigation_model.navigation_rates))
end

function _validate_external_force_service!(args, ep)
    _external_capabilities(ep).configured_forces || return nothing
    length(args.dynamics_model.spacecraft) == 1 &&
        SimulationModel.ExternalPropagation.external_spacecraft(ep) == [1] ||
        throw(ArgumentError("Experimental configured external forces support exactly one spacecraft."))
    !args.mission_configuration.orientation_sim ||
        throw(ArgumentError("Experimental configured external forces are translation-only; orientation_sim must be false."))
    args.environment_model.ephemerides_model isa SimulationModel.SimpleEphemeridesModel ||
        throw(ArgumentError("Experimental configured external forces currently require SimpleEphemeridesModel; other frame backends are unqualified."))
    spacecraft = only(args.dynamics_model.spacecraft)
    SimulationModel.mounted_actuator_count(spacecraft) == 0 &&
        isempty(spacecraft.root.rw_assembly.mounted_rates) && isempty(spacecraft.root.rw_assembly.cmg_rates) ||
        throw(ArgumentError("Experimental configured external forces do not support mounted actuator hardware, including wheel/CMG rate buffers."))
    isempty(args.guidance_model.guidance_effectors) && isempty(args.navigation_model.navigation_effectors) ||
        throw(ArgumentError("Experimental configured external forces currently require direct-state control without guidance/navigation effectors."))
    args.interlink_model === nothing || throw(ArgumentError("External configured forces do not support interlinks."))
    for model in args.dynamics_model.dynamic_effectors
        _external_stage_safe(model) || throw(ArgumentError("External configured forces do not support dynamic effector $(typeof(model)); stage safety is not declared."))
    end
    for model in args.control_model.control_effectors
        SimulationModel.actuation_mode(model) === :ideal &&
            SimulationModel.SimulationCallbacks.tick_schedule(model) === :periodic &&
            _external_stage_safe(model) || throw(ArgumentError("External configured forces require stage-safe, ideal, periodic control effectors; got $(typeof(model))."))
    end
    _check_external_lifecycle_cadences!(args, SimulationModel.ExternalPropagation.external_step(ep))
    _external_force_configuration(args)
    return nothing
end

# Lifecycle times are exact. With a configured-force owner the owner step, the control
# periods, the output interval, the mission end and (with checkpoints) the checkpoint interval
# must be binary-exact: whole multiples of 2^-20 s below 2^33 s. Sums and whole multiples of such
# values are exact in Float64, so checkpoint segment ends (built by adding the interval), the
# output grid (built by index), the control ticks and the owner's step grid land on the same
# Float64 values, and `floor(t / step)` is exact at every tick, checkpoint and output time. A
# decimal cadence such as 0.1 s drifts by a few ulps between those constructions, which turns an
# output time into a stop one ulp after a segment start; it is refused until it is supported.
const _EXTERNAL_TIME_QUANTUM_EXPONENT = 20

function _binary_exact_time(x::Real)
    v = Float64(x)
    q = ldexp(v, _EXTERNAL_TIME_QUANTUM_EXPONENT)
    return isfinite(v) && v > 0.0 && isinteger(q) && q < 2.0^53
end

function _exact_step_multiple(x::Real, dt::Float64)
    _binary_exact_time(x) || return false
    k = Float64(x) / dt
    return isinteger(k) && k >= 1.0 && k * dt == Float64(x)
end

function _check_external_lifecycle_cadences!(args, dt::Float64)
    _binary_exact_time(dt) || throw(ArgumentError(
        "Experimental configured external forces require a binary-exact owner step (a whole multiple of 2^-20 s, " *
        "such as 1/128 s or 10 s); got $dt s. Decimal cadences are not supported yet."))
    for rate in _periodic_control_rates(args)
        _exact_step_multiple(rate, dt) || throw(ArgumentError(
            "Experimental configured external forces require each control period to be an exact whole multiple of the " *
            "owner step $dt s; got $rate s."))
    end
    mission = args.mission_configuration
    _binary_exact_time(mission.data_rate) || throw(ArgumentError(
        "Experimental configured external forces require a binary-exact output interval data_rate (a whole multiple of " *
        "2^-20 s); got $(mission.data_rate) s."))
    _binary_exact_time(mission.mission_time) || throw(ArgumentError(
        "Experimental configured external forces require a binary-exact mission end (a whole multiple of 2^-20 s); " *
        "got $(mission.mission_time) s."))
    if _typed_checkpoint_enabled(args)
        settings = args.simulation_settings
        _exact_step_multiple(settings.checkpoint_interval_s, dt) && _exact_step_multiple(mission.mission_time, dt) ||
            throw(ArgumentError("Experimental external checkpoints require the checkpoint interval and the mission end to be " *
                "exact whole multiples of the owner step $dt s; got $(settings.checkpoint_interval_s) s and $(mission.mission_time) s."))
    end
    return nothing
end

"""Owned run indices over every external propagator in `args` (empty when there are none)."""
function _external_owned_indices(args)::Vector{Int}
    owned = Int[]
    for ep in _external_propagators(args)
        append!(owned, SimulationModel.ExternalPropagation.external_spacecraft(ep))
    end
    return owned
end

# Does the control effector act on one of the owned spacecraft (or on every spacecraft)?
function _control_effector_targets_owned(effector, owned::Vector{Int})::Bool
    for name in _EXTERNAL_INDEX_FIELDS
        hasproperty(effector, name) || continue
        v = getproperty(effector, name)
        v isa Integer && return Int(v) in owned
    end
    if effector isa SimulationModel.ThrusterModels.BaseThrusterModel
        n = length(effector.thrust)
        return any(i -> i <= n && (effector.thrust[i] != 0.0 || effector.Δv[i] != 0.0), owned)
    end
    return true     # no spacecraft index: the effector acts on every spacecraft
end

@inline function _is_integer_multiple(rate::Real, dt::Float64)::Bool
    r = Float64(rate) / dt
    return isfinite(r) && r >= 1.0 - 1.0e-9 && abs(r - round(r)) <= 1.0e-9 * max(1.0, r)
end

"""
    _validate_external_propagators!(args, solver_mode)

Setup-time guards for externally propagated spacecraft (each throws an `ArgumentError`). A no-op when
`args.external_propagators` is empty. Legacy owners refuse checkpointing, resume and actuation.
Experimental lifecycle owners may explicitly opt in under the narrow force/grid/snapshot guards below.
All owners refuse the solver modes other than
`:tsit5`, `:auto_stiff`, `:rodas5p` and `:dp8`; articulated joints, compliant attachments or a robot arm on an
owned spacecraft; unsupported actuation; a guidance, navigation or control rate that is not an integer multiple of an owner's step;
and the continuous events that would evaluate shadow entries (touchdown, orbit-count and entry-end
termination, atmosphere-interface crossings).
"""
function _validate_external_propagators!(args, solver_mode::Symbol)
    eps = _external_propagators(args)
    isempty(eps) && return nothing
    EP = SimulationModel.ExternalPropagation
    n = length(args.dynamics_model.spacecraft)
    all(ep -> ep isa EP.AbstractExternalPropagator, eps) || throw(ArgumentError(
        "external_propagators must hold SimulationModel.ExternalPropagation.AbstractExternalPropagator values; got $(map(typeof, eps))."))
    owned = Int[]
    steps = Float64[]
    for ep in eps
        idx = EP.external_spacecraft(ep)
        isempty(idx) && throw(ArgumentError("External propagator $(nameof(typeof(ep))) owns no spacecraft."))
        all(i -> 1 <= i <= n, idx) || throw(ArgumentError(
            "External propagator $(nameof(typeof(ep))) owns spacecraft $(idx), but the run has $n spacecraft."))
        append!(owned, idx)
        dt = EP.external_step(ep)
        (isfinite(dt) && dt > 0.0) || throw(ArgumentError("External propagator $(nameof(typeof(ep))) step must be positive and finite; got $dt."))
        push!(steps, dt)
    end
    allunique(owned) || throw(ArgumentError("A spacecraft cannot be owned by more than one external propagator or listed twice; got owned indices $(sort(owned))."))
    if _typed_checkpoint_enabled(args) && any(ep -> !_external_capabilities(ep).checkpoint, eps)
        throw(ArgumentError(
            "Externally propagated spacecraft do not support checkpointing or resume: the external state is not part of the checkpoint. " *
            "Disable simulation_settings.checkpoint_enabled and resume_from_checkpoint."))
    end
    for ep in eps
        caps = _external_capabilities(ep)
        # Lifecycle extension 0.1.0 qualifies checkpoints only with configured forces: a checkpoint-only
        # owner is synchronized only when it advances, so its saved rows between steps would be the
        # shadow integration's.
        caps.checkpoint && !caps.configured_forces && throw(ArgumentError(
            "External propagator $(nameof(typeof(ep))) declares checkpoint without configured_forces; " *
            "the experimental lifecycle extension supports checkpoints only together with configured forces."))
        _validate_external_force_service!(args, ep)
        if _typed_checkpoint_enabled(args)
            dt = EP.external_step(ep)
            settings = args.simulation_settings
            _is_integer_multiple(settings.checkpoint_interval_s, dt) &&
                _is_integer_multiple(args.mission_configuration.mission_time, dt) ||
                throw(ArgumentError("Experimental external checkpoints require checkpoint intervals and mission end on the owner step grid."))
            EP.external_checkpoint_key(ep) # Missing snapshot identity fails before runtime construction.
        end
    end
    if !(solver_mode in _EXTERNAL_ALLOWED_SOLVER_MODES)
        throw(ArgumentError(
            "Externally propagated spacecraft support the first-order single-RHS solver modes :tsit5, :auto_stiff, :rodas5p and :dp8; " *
            "solver mode $(repr(solver_mode)) (split, multirate, symplectic and gravity-backbone routes) is not supported."))
    end
    for i in owned
        sc = args.dynamics_model.spacecraft[i]
        _spacecraft_articulated(sc) && throw(ArgumentError(
            "Spacecraft $i is owned by an external propagator and has articulated joints; articulation belongs to the external scene, not to a SpaceAGORA joint tree."))
        _spacecraft_has_attachments(sc) && throw(ArgumentError(
            "Spacecraft $i is owned by an external propagator and has compliant attachments; the two cannot be combined."))
        _robot_arm_coupling(args, i, 0.0) !== nothing && throw(ArgumentError(
            "Spacecraft $i is owned by an external propagator and has a compliant robot-arm effector; the two cannot be combined."))
    end
    for effector in args.control_model.control_effectors
        _control_effector_targets_owned(effector, owned) &&
            any(ep -> !_external_capabilities(ep).configured_forces, eps) && throw(ArgumentError(
            "Control effector $(nameof(typeof(effector))) acts on a spacecraft owned by an external propagator (or on every spacecraft). " *
            "Actuating externally propagated spacecraft is not supported yet; read their state from guidance or navigation only, " *
            "and bind control effectors to ordinary spacecraft through SpacecraftModel.control."))
    end
    for (kind, rates) in (("guidance", args.guidance_model.guidance_rates), ("navigation", args.navigation_model.navigation_rates),
            ("control", _periodic_control_rates(args)))
        for rate in rates, dt in steps
            _is_integer_multiple(rate, dt) || throw(ArgumentError(
                "A $kind rate of $rate s is not an integer multiple of the external propagator step $dt s; " *
                "periodic callbacks must land on step boundaries."))
        end
    end
    cb = SimulationModel.SimulationCallbacks
    touchdown = cb._touchdown_specs(args, n)
    any(i -> touchdown[i] !== nothing, owned) && throw(ArgumentError(
        "A touchdown event is configured for a spacecraft owned by an external propagator; continuous events on shadow entries are not supported."))
    effectors = args.dynamics_model.dynamic_effectors
    cb._requires_orbit_end_callback(args) && throw(ArgumentError(
        "Orbit-count termination (MissionOrbits or a guidance maneuver_orbit_number) is a continuous event on every spacecraft, shadow entries included; " *
        "use a time-bounded mission with externally propagated spacecraft."))
    (cb._requires_entry_end_callback(effectors, args) || cb._requires_drag_state_callback(effectors, args)) && throw(ArgumentError(
        "Atmosphere-interface and entry events are continuous events on every spacecraft, shadow entries included; " *
        "they are not supported with externally propagated spacecraft."))
    for ep in eps
        EP.external_preflight(ep, args)
    end
    return nothing
end

# A trial view belongs only to this service. The runtime/configuration is this run's, so
# force evaluation sees the current held command without changing committed spacecraft state.
# The narrow preflight excludes history-dependent, mounted and coupled effectors.
struct ExternalPhysicalForceService{P,U}
    parameters::P
    trial::U
    satellite::Int
end

function (service::ExternalPhysicalForceService)(pos, vel, t::Real)
    length(pos) == length(vel) == 3 || throw(ArgumentError("External trial position and velocity must have three components."))
    all(isfinite, pos) && all(isfinite, vel) && isfinite(t) || throw(ArgumentError("External trial state/time must be finite."))
    p = service.parameters; i = service.satellite; x = service.trial
    x.pos .= pos; x.vel .= vel
    forces = MVector{3,Float64}(0.0,0.0,0.0)
    torques = MVector{3,Float64}(0.0,0.0,0.0)
    sample = build_state_sample(x, p.args.dynamics_model.spacecraft[i], false, i)
    # The serial effector sum of the run's right-hand side, including its effector
    # scopes: an effector whose scope leaves this spacecraft out adds nothing.
    effectors = p.args.dynamics_model.dynamic_effectors
    scopes = _effector_scopes(p)
    if scopes === nothing
        _accumulate_effector_chain!(forces, torques, effectors, x, sample, p, i, Float64(t))
    else
        _check_scoped_tuple(scopes, effectors)
        _accumulate_scoped_effector_chain!(forces, torques, effectors, scopes, 1, x, sample, p, i, Float64(t))
    end
    # The ordinary assembly clamps nonfinite propellant rates. The external contract is stricter:
    # a declared constant-mass model must report a finite zero rate, not rely on that fallback.
    for model in p.args.control_model.control_effectors
        propellant = SimulationModel.control_propellant_rate(model, sample, p, Float64(t))
        isfinite(propellant) && iszero(propellant) || throw(ArgumentError("External force service requires a finite zero propellant rate."))
    end
    mass_rate, wheel_torque = _control_effectors_with_wheels!(forces, torques, x, p, i, Float64(t), false)
    mass_rate == 0.0 && all(iszero, torques) && all(iszero, wheel_torque) ||
        throw(ArgumentError("External configured-force service refuses mass loss or torque; the opted-in model violated its stage-safe contract."))
    mass = Float64(x.mass)
    isfinite(mass) && mass > 0 || throw(ArgumentError("External force service requires positive finite constant mass."))
    acceleration = SimulationModel.DynamicsTranslational.acceleration_from_force(forces, mass)
    all(isfinite, acceleration) || throw(ArgumentError("External force service produced nonfinite acceleration."))
    return acceleration
end

function _bind_external_force_services!(p, u)
    EP = SimulationModel.ExternalPropagation
    for (k, ep) in enumerate(_external_propagators(p.args))
        _external_capabilities(ep).configured_forces || continue
        i = only(EP.external_spacecraft(ep))
        service = ExternalPhysicalForceService(p, deepcopy(u.sc[i]), i)
        EP.external_bind_force_service!(p.shared_buffers.external_runtimes[k], service)
    end
    return nothing
end

"""Use the synchronized owner as the authority for a checkpoint segment's physical endpoint.
The solver's saved endpoint can differ from its callback-corrected live state by rounding.
Validate every owner and collect every update before changing the saved endpoint. Never advance an
owner here: a control tick may already have installed the next interval's command.
The returned diagnostic retains the pre-handoff discrepancy; only owned kinematic fields change.
"""
function _capture_external_checkpoint_endpoint!(p, u, t::Float64)
    isempty(_external_propagators(p.args)) && return nothing
    EP = SimulationModel.ExternalPropagation
    updates = NamedTuple[]
    for (k, ep) in enumerate(_external_propagators(p.args))
        rt = p.shared_buffers.external_runtimes[k]
        owner_t = EP.external_time(rt)
        isfinite(owner_t) && abs(owner_t-t) <= 16eps(max(abs(t),1.0)) ||
            throw(ArgumentError("External owner is not synchronized at checkpoint handoff time $t."))
        for (slot, i) in enumerate(EP.external_spacecraft(ep))
            state = EP.external_state(rt,t,slot)
            pos = collect(Float64, state.pos); vel = collect(Float64, state.vel)
            length(pos)==length(vel)==3 && all(isfinite,pos) && all(isfinite,vel) ||
                throw(ArgumentError("Invalid external kinematics at checkpoint handoff."))
            attitude = hasproperty(u.sc[i], :q)
            q = attitude ? collect(Float64,state.q) : nothing
            omega = attitude ? collect(Float64,state.ω) : nothing
            (!attitude || (length(q)==4 && length(omega)==3 && all(isfinite,q) && all(isfinite,omega))) ||
                throw(ArgumentError("Invalid external attitude at checkpoint handoff."))
            push!(updates,(satellite=i,owner_time_s=owner_t,before_pos=collect(u.sc[i].pos),before_vel=collect(u.sc[i].vel),
                pos=pos,vel=vel,before_q=attitude ? collect(u.sc[i].q) : nothing,
                before_omega=attitude ? collect(u.sc[i].ω) : nothing,q=q,omega=omega))
        end
    end
    for update in updates
        sc = u.sc[update.satellite]
        sc.pos .= update.pos; sc.vel .= update.vel
        if update.q !== nothing
            sc.q .= update.q; sc.ω .= update.omega
        end
    end
    return isempty(updates) ? nothing : (version=1,time_s=t,updates=updates)
end

function _external_checkpoint_state(p, t::Float64; handoff=nothing)
    isempty(_external_propagators(p.args)) && return nothing
    EP = SimulationModel.ExternalPropagation
    owners = Any[]
    for (k, ep) in enumerate(_external_propagators(p.args))
        _external_capabilities(ep).checkpoint || throw(ArgumentError("External owner did not opt in to checkpointing."))
        rt = p.shared_buffers.external_runtimes[k]
        abs(EP.external_time(rt)-t) <= 16eps(max(abs(t),1.0)) ||
            throw(ArgumentError("External owner is not synchronized at checkpoint time $t."))
        push!(owners, (key=deepcopy(EP.external_checkpoint_key(ep)), owned=copy(EP.external_spacecraft(ep)),
            step=EP.external_step(ep), force_configuration=_external_capabilities(ep).configured_forces ? deepcopy(_external_force_configuration(p.args)) : nothing,
            state=deepcopy(EP.external_checkpoint_state(rt,t))))
    end
    return isempty(owners) ? nothing : (version=1, time_s=t, owners=owners, handoff=handoff)
end

function _check_external_checkpoint_state(args, checkpoint)
    EP = SimulationModel.ExternalPropagation
    eps = _external_propagators(args)
    runtime = get(checkpoint, :runtime_state, nothing)
    entry = runtime === nothing ? nothing : get(runtime, :external_owners, nothing)
    isempty(eps) && return entry === nothing ? nothing : throw(ArgumentError("Checkpoint contains external owners absent from this configuration."))
    entry !== nothing && entry.version == 1 && entry.time_s == checkpoint.t && length(entry.owners) == length(eps) ||
        throw(ArgumentError("Checkpoint is missing compatible external-owner state."))
    for (ep, saved) in zip(eps, entry.owners)
        _external_capabilities(ep).checkpoint && isequal(saved.key, EP.external_checkpoint_key(ep)) &&
            saved.owned == EP.external_spacecraft(ep) && saved.step == EP.external_step(ep) ||
            throw(ArgumentError("External checkpoint owner identity, ownership or step mismatch."))
        expected = _external_capabilities(ep).configured_forces ? _external_force_configuration(args) : nothing
        isequal(saved.force_configuration, expected) || throw(ArgumentError("External checkpoint configured physics or control parameters differ."))
    end
    return entry
end

function _restore_external_checkpoint_state!(p, entry, u, t::Float64)
    entry === nothing && return nothing
    EP = SimulationModel.ExternalPropagation
    for (k, saved) in enumerate(entry.owners)
        rt = p.shared_buffers.external_runtimes[k]
        EP.external_restore_state!(rt, deepcopy(saved.state), t)
        abs(EP.external_time(rt)-t) <= 16eps(max(abs(t),1.0)) || throw(ArgumentError("Restored external clock differs from checkpoint."))
        for (slot, i) in enumerate(saved.owned)
            state = EP.external_state(rt,t,slot)
            state.pos == u.sc[i].pos && state.vel == u.sc[i].vel || throw(ArgumentError(
                "Restored external physical state differs from checkpoint shadow state at t=$t: " *
                "owner_pos=$(collect(state.pos)), shadow_pos=$(collect(u.sc[i].pos)), " *
                "owner_vel=$(collect(state.vel)), shadow_vel=$(collect(u.sc[i].vel))."))
            if hasproperty(u.sc[i], :q)
                state.q == u.sc[i].q && state.ω == u.sc[i].ω || throw(ArgumentError("Restored external attitude differs from checkpoint shadow state."))
            end
        end
    end
    return nothing
end

function _periodic_control_rates(args)
    rates = Float64[]
    for (effector, rate) in zip(args.control_model.control_effectors, args.control_model.control_rates)
        SimulationModel.SimulationCallbacks.tick_schedule(effector) === :periodic && push!(rates, rate)
    end
    return rates
end

"""Build the run's external runtimes from the initial state `u0` and the per-satellite owner tables."""
function _initialize_external_runtimes!(p, u0)
    args = p.args
    buffers = p.shared_buffers
    n_sats = length(args.dynamics_model.spacecraft)
    empty!(buffers.external_runtimes)
    resize!(buffers.external_owner, n_sats); fill!(buffers.external_owner, 0)
    resize!(buffers.external_local, n_sats); fill!(buffers.external_local, 0)
    buffers.external_present[] = false
    isempty(_external_propagators(args)) && return nothing
    EP = SimulationModel.ExternalPropagation
    for (r, ep) in enumerate(_external_propagators(args))
        push!(buffers.external_runtimes, EP.external_prepare(ep, args, u0))
        for (k, i) in enumerate(EP.external_spacecraft(ep))
            buffers.external_owner[i] = r
            buffers.external_local[i] = k
        end
    end
    buffers.external_present[] = true
    if _rhs_env_config(p).execution_mode == :flat_constellation_effector_queue
        throw(ArgumentError(
            "SPACEAGORA_RHS_EXECUTION_MODE=flat is not supported with externally propagated spacecraft; use auto, serial or satellite."))
    end
    return nothing
end

# Shadow right-hand side: the owner's chief acceleration (at time t) for translation, plain kinematics for attitude with
# constant body rate (the owner's state replaces it at the next sync); mass and heat load are constant.
@inline function _assign_external_shadow_rhs!(du_view, sc_view, p, sat_idx::Int, t::Float64)
    buffers = p.shared_buffers
    rt = buffers.external_runtimes[buffers.external_owner[sat_idx]]
    a = SVector{3, Float64}(SimulationModel.ExternalPropagation.external_acceleration(rt, t, buffers.external_local[sat_idx]))
    du_view .= 0.0
    @inbounds for k in 1:3
        du_view.pos[k] = sc_view.vel[k]
        du_view.vel[k] = a[k]
    end
    if hasproperty(sc_view, :q)
        du_view.q .= SimulationModel.DynamicsRotational.quaternion_derivative(
            SimulationModel.DynamicsRotational.body_angular_velocity(sc_view.ω), sc_view.q)
    end
    return nothing
end
