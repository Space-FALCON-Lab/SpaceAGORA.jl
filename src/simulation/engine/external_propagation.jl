# Externally propagated spacecraft ("shadow entries"): setup-time refusals, per-run runtimes and the shadow
# right-hand side. The hook surface is `SimulationModel.ExternalPropagation`; the sync callback that overwrites
# the shadow entries after every accepted step lives in `callbacks/external_propagation_callbacks.jl`.
#
# A shadow entry stays in `u.sc` so default save fields, guidance, navigation and `bind_spacecraft` see it like
# any other spacecraft. Between syncs it follows the owner's chief acceleration and plain attitude kinematics;
# at each sync it is overwritten from the owner's absolute state.

const _EXTERNAL_ALLOWED_SOLVER_MODES = (:tsit5, :auto_stiff, :rodas5p, :dp8)
const _EXTERNAL_INDEX_FIELDS = (:chaser_idx, :spacecraft_idx, :sat_idx)

@inline _external_propagators(args) = args.external_propagators

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
`args.external_propagators` is empty. Refused: checkpointing or resume; the solver modes other than
`:tsit5`, `:auto_stiff`, `:rodas5p` and `:dp8`; articulated joints, compliant attachments or a robot arm on an
owned spacecraft; a control effector that acts on an owned spacecraft (actuation arrives with the scene
actuation stage); a guidance, navigation or control rate that is not an integer multiple of an owner's step;
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
    if _typed_checkpoint_enabled(args)
        throw(ArgumentError(
            "Externally propagated spacecraft do not support checkpointing or resume: the external state is not part of the checkpoint. " *
            "Disable simulation_settings.checkpoint_enabled and resume_from_checkpoint."))
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
            "Spacecraft $i is owned by an external propagator and has a cloth robot-arm effector; the two cannot be combined."))
    end
    for effector in args.control_model.control_effectors
        _control_effector_targets_owned(effector, owned) && throw(ArgumentError(
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

function _periodic_control_rates(args)
    rates = Float64[]
    for (effector, rate) in zip(args.control_model.control_effectors, args.control_model.control_rates)
        SimulationModel.SimulationCallbacks.control_requires_periodic_callback(effector) && push!(rates, rate)
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
    a = SVector{3, Float64}(SimulationModel.ExternalPropagation.external_acceleration(rt, buffers.external_local[sat_idx], t))
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
