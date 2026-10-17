"""
Closed-loop magnetic momentum management as a discrete-rate control effector.

Reaction wheels absorb the commanded attitude-control torque, so their stored
momentum integrates every secular disturbance and grows without bound unless
it is routed elsewhere. Flight software periodically commands magnetic torque
rods with the cross-product unloading law

    m = mu * (H_w x B) / |B|^2,      tau_rod = m x B = -mu * H_w_perp,

which bleeds the component of wheel momentum perpendicular to the local field
(the parallel component is unreachable by magnetics and must wait for the
field direction to rotate along the orbit).

[`MagneticMomentumManagerModel`](@ref) reproduces that loop as a stateful
zero-order-hold control effector on the discrete `ControlModel` callback path:
`update_control!` fires at the configured control rate (a
`PeriodicCallback`, i.e. only at accepted solver steps), advances the internal
wheel-momentum accumulator from the re-evaluated attitude-control command,
and refreshes the held rod torque; `control_wrench` returns that held
torque on every RHS evaluation. Keeping the accumulator on the discrete path
is what makes it correct under adaptive stepping — a continuous effector
cannot hold integration state because trial evaluations of rejected steps
would corrupt it (see the note in `perturbations.jl` on accumulator state).
The zero-order hold also mirrors how rod commands are actually flown.

The model is deliberately decoupled from any specific attitude controller and
field model through two user-supplied callables:

- `commanded_torque(t, r_ii, v_ii, q_ib, w_body) -> SVector{3}`: the
  body-frame attitude-control torque the wheels are absorbing. Wrap the
  simulation's actual controller (e.g. `_lvlh_cascade_torque` for
  [`LVLHCascadeAttitudeControlModel`](@ref)) so the accumulator sees the same
  command the dynamics apply.
- `b_field_ii(t, r_ii) -> SVector{3}`: the inertial-frame magnetic field in
  tesla (a dipole/IGRF wrapper or an interpolated table).

Wheel momentum uses the convention `dH_w/dt = -tau_cmd` (wheels react to the
body torque they produce). When the spacecraft carries reaction wheels whose
momentum the engine integrates (its root wheel assembly, for example wheels
mounted with `add_actuator!`), the manager reads that momentum at each tick
instead of accumulating its own. The dipole command is capped at `m_max_am2`.

`actuation = :ideal` (the default) holds the rod torque `m x B` computed with
the model's own field callable. `actuation = :mounted` instead drives the
torque rods mounted on the spacecraft (`Magnetorquer` actuators) at the
dipole command, so the torque comes from the run's field (the
`MagneticTorqueRodModel` field source) while the model's callable stays its
estimate of that field.

To read `h_wheels`/`held_dipole_am2` back after a run, call
`run_simulation(args; isolate_state=false)`: the default `isolate_state=true`
deep-copies the configuration — including this model — so the simulation
advances a copy and the caller's handle is left untouched. Each run starts the
accumulator, the held command and the tick count from their values before the
model's first run (`initialize_control!` restores them), so two runs of one
configuration with `isolate_state=false` start equal.

`spacecraft_index` is the run position of the managed spacecraft; `sat_idx`, its
old name, still works with a deprecation warning until SpaceAGORA 0.5.0.
"""
mutable struct MagneticMomentumManagerModel <: AbstractControlEffectorModel
    spacecraft_index::Int
    mu_gain::Float64                          # unloading gain [1/s]
    m_max_am2::Float64                        # rod dipole capability cap [A m^2]
    h_wheels_0::SVector{3, Float64}
    commanded_torque::Any                     # (t, r_ii, v_ii, q_ib, w_body) -> tau_body [N m]
    b_field_ii::Any                           # (t, r_ii) -> B_ii [T]
    # internal discrete state (updated only inside update_control!)
    h_wheels::SVector{3, Float64}
    held_torque_body::SVector{3, Float64}
    held_dipole_am2::SVector{3, Float64}
    last_update_s::Float64
    initialized::Bool
    ticks::Int                                # diagnostic: discrete updates applied
    actuation::Symbol                         # :ideal holds the rod torque; :mounted drives mounted torque rods
    mounted_levels::Vector{Float64}           # :mounted: one level per mounted actuator, refreshed each tick
    run_start::Union{Nothing, NamedTuple}     # the internal state before the first run, restored at each run's start
end
# The deprecated name forwards, as a keyword and as a property.
function MagneticMomentumManagerModel(; spacecraft_index=_NOT_PASSED, sat_idx=_NOT_PASSED, mu_gain=3.0e-4, m_max_am2=Inf,
        h_wheels_0=SVector{3, Float64}(0.0, 0.0, 0.0), commanded_torque=nothing, b_field_ii=nothing,
        h_wheels=SVector{3, Float64}(0.0, 0.0, 0.0), held_torque_body=SVector{3, Float64}(0.0, 0.0, 0.0),
        held_dipole_am2=SVector{3, Float64}(0.0, 0.0, 0.0), last_update_s=NaN, initialized=false, ticks=0,
        actuation=:ideal, mounted_levels=Float64[])
    T = MagneticMomentumManagerModel
    return T(_renamed_keyword(T, :sat_idx, sat_idx, :spacecraft_index, spacecraft_index, 1),
        mu_gain, m_max_am2, h_wheels_0, commanded_torque, b_field_ii, h_wheels, held_torque_body, held_dipole_am2,
        last_update_s, initialized, ticks, actuation, mounted_levels, nothing)
end
@inline _momentum_manager_field(name::Symbol) = name === :sat_idx ?
    _renamed_field_warning(MagneticMomentumManagerModel, :sat_idx, :spacecraft_index) : name
@inline Base.getproperty(m::MagneticMomentumManagerModel, name::Symbol) = getfield(m, _momentum_manager_field(name))
@inline function Base.setproperty!(m::MagneticMomentumManagerModel, name::Symbol, x)
    field = _momentum_manager_field(name)
    return setfield!(m, field, convert(fieldtype(MagneticMomentumManagerModel, field), x))
end

# The internal state a run changes. The first run records it; every later run
# starts from the recorded values (G-8). The mounted levels start empty, as in a
# bound copy.
_manager_run_state(m::MagneticMomentumManagerModel) = (h_wheels=m.h_wheels, held_torque_body=m.held_torque_body,
    held_dipole_am2=m.held_dipole_am2, last_update_s=m.last_update_s, initialized=m.initialized, ticks=m.ticks)

function SimulationLifecycle.initialize_control!(m::MagneticMomentumManagerModel, u, p, t)
    if m.run_start === nothing
        m.run_start = _manager_run_state(m)
    else
        s = m.run_start
        m.h_wheels, m.held_torque_body, m.held_dipole_am2 = s.h_wheels, s.held_torque_body, s.held_dipole_am2
        m.last_update_s, m.initialized, m.ticks = s.last_update_s, s.initialized, s.ticks
    end
    empty!(m.mounted_levels)
    return nothing
end

# The body-frame momentum of the wheels the engine integrates for spacecraft
# `sat_idx` (its root wheel assembly), or `nothing` when it carries none.
function _integrated_wheel_momentum(u, p, sat_idx::Int)
    hasproperty(p, :args) || return nothing
    assembly = p.args.dynamics_model.spacecraft[sat_idx].root.rw_assembly
    (assembly.n_wheels > 0 && hasproperty(u.sc[sat_idx], :h_wheels)) || return nothing
    return SVector{3, Float64}(assembly.J_rw * u.sc[sat_idx].h_wheels)
end

"""Advance the wheel-momentum accumulator (or read the integrated wheels) and refresh the held rod command."""
function update_control!(model::MagneticMomentumManagerModel, u, p, t::Float64, sat_idx::Int)
    sat_idx == model.spacecraft_index || return nothing
    (model.commanded_torque === nothing || model.b_field_ii === nothing) && return nothing

    r_ii = SVector{3, Float64}(u.sc[sat_idx].pos)
    v_ii = SVector{3, Float64}(u.sc[sat_idx].vel)
    q_ib = SVector{4, Float64}(u.sc[sat_idx].q)
    w_body = SVector{3, Float64}(u.sc[sat_idx].ω)

    model.ticks += 1
    tau_cmd = SVector{3, Float64}(model.commanded_torque(t, r_ii, v_ii, q_ib, w_body))
    integrated = _integrated_wheel_momentum(u, p, sat_idx)
    if integrated !== nothing
        # The vehicle carries wheels whose momentum the engine integrates: use it.
        model.h_wheels = integrated
        model.initialized = true
    elseif !model.initialized
        model.h_wheels = model.h_wheels_0
        model.initialized = true
    else
        dt = t - model.last_update_s
        if isfinite(dt) && dt > 0.0
            model.h_wheels = model.h_wheels - tau_cmd * dt
        end
    end
    model.last_update_s = t

    b_ii = SVector{3, Float64}(model.b_field_ii(t, r_ii))
    b_body = rot(q_ib) * b_ii
    b2 = dot(b_body, b_body)
    if !(b2 > 0.0)
        model.held_dipole_am2 = SVector{3, Float64}(0.0, 0.0, 0.0)
        model.held_torque_body = SVector{3, Float64}(0.0, 0.0, 0.0)
        model.actuation === :mounted && _manager_fill_mounted_levels!(model, p.args.dynamics_model.spacecraft[sat_idx], model.held_dipole_am2)
        return nothing
    end
    m = model.mu_gain * cross(model.h_wheels, b_body) / b2
    nm = norm(m)
    if isfinite(model.m_max_am2) && nm > model.m_max_am2
        m = m * (model.m_max_am2 / nm)
    end
    model.held_dipole_am2 = m
    model.held_torque_body = cross(m, b_body)
    model.actuation === :mounted && _manager_fill_mounted_levels!(model, p.args.dynamics_model.spacecraft[sat_idx], m)
    return nothing
end

# `:mounted`: the dipole command as levels of the mounted torque rods (their
# least-norm split, so three orthogonal rods carry its components); every other
# mounted actuator idle.
function _manager_fill_mounted_levels!(model::MagneticMomentumManagerModel, spacecraft, dipole_body::SVector{3, Float64})
    n_rods = mounted_magnetorquer_count(spacecraft)
    n_rods > 0 || throw(ArgumentError("Mounted momentum management needs torque rods mounted on the spacecraft: add Magnetorquer actuators."))
    axes = zeros(Float64, 3, n_rods)
    k = 0
    for link in spacecraft.links, rod in link.magnetorquers
        k += 1
        axes[:, k] = actuator_context(link).link_to_body * rod.axis
    end
    split = pinv(axes) * dipole_body
    levels = model.mounted_levels
    resize!(levels, mounted_actuator_count(spacecraft))
    fill!(levels, 0.0)
    offset = mounted_thruster_count(spacecraft) + mounted_wheel_count(spacecraft)
    for j in 1:n_rods
        levels[offset + j] = split[j]
    end
    return levels
end

function actuation_mode(model::MagneticMomentumManagerModel)
    model.actuation === :ideal || model.actuation === :mounted ||
        throw(ArgumentError("MagneticMomentumManagerModel actuation must be :ideal or :mounted, got $(repr(model.actuation))."))
    return model.actuation
end

function actuator_levels(model::MagneticMomentumManagerModel, i::Int)
    (i == model.spacecraft_index && !isempty(model.mounted_levels)) || return nothing
    return model.mounted_levels
end

"""Return the held rod torque (body frame); no force."""
function control_wrench(model::MagneticMomentumManagerModel, x::StateSample, p, t::Float64)
    zero3 = SVector{3, Float64}(0.0, 0.0, 0.0)
    x.spacecraft_index == model.spacecraft_index || return zero3, zero3
    return zero3, model.held_torque_body
end

# No control_propellant_rate method: torque rods consume no propellant, and
# the dynamics then take the rate as zero.
