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
`calcControlEffect!` fires at the configured control rate (a
`PeriodicCallback`, i.e. only at accepted solver steps), advances the internal
wheel-momentum accumulator from the re-evaluated attitude-control command,
and refreshes the held rod torque; `calcControlForceTorque` returns that held
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
body torque they produce). The dipole command is capped at `m_max_am2`.

To read `h_wheels`/`held_dipole_am2` back after a run, call
`run_simulation(args; isolate_state=false)`: the default `isolate_state=true`
deep-copies the configuration — including this model — so the simulation
advances a copy and the caller's handle is left untouched.
"""
Base.@kwdef mutable struct MagneticMomentumManagerModel <: AbstractControlEffectorModel
    sat_idx::Int = 1
    mu_gain::Float64 = 3.0e-4                 # unloading gain [1/s]
    m_max_am2::Float64 = Inf                  # rod dipole capability cap [A m^2]
    h_wheels_0::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 0.0)
    commanded_torque::Any = nothing           # (t, r_ii, v_ii, q_ib, w_body) -> tau_body [N m]
    b_field_ii::Any = nothing                 # (t, r_ii) -> B_ii [T]
    # internal discrete state (updated only inside calcControlEffect!)
    h_wheels::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 0.0)
    held_torque_body::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 0.0)
    held_dipole_am2::SVector{3, Float64} = SVector{3, Float64}(0.0, 0.0, 0.0)
    last_update_s::Float64 = NaN
    initialized::Bool = false
    ticks::Int = 0                            # diagnostic: discrete updates applied
end

"""Advance the wheel-momentum accumulator and refresh the held rod command."""
function calcControlEffect!(model::MagneticMomentumManagerModel, u, p, t::Float64, sat_idx::Int)
    sat_idx == model.sat_idx || return nothing
    (model.commanded_torque === nothing || model.b_field_ii === nothing) && return nothing

    r_ii = SVector{3, Float64}(u.sc[sat_idx].pos)
    v_ii = SVector{3, Float64}(u.sc[sat_idx].vel)
    q_ib = SVector{4, Float64}(u.sc[sat_idx].q)
    w_body = SVector{3, Float64}(u.sc[sat_idx].ω)

    model.ticks += 1
    tau_cmd = SVector{3, Float64}(model.commanded_torque(t, r_ii, v_ii, q_ib, w_body))
    if !model.initialized
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
        return nothing
    end
    m = model.mu_gain * cross(model.h_wheels, b_body) / b2
    nm = norm(m)
    if isfinite(model.m_max_am2) && nm > model.m_max_am2
        m = m * (model.m_max_am2 / nm)
    end
    model.held_dipole_am2 = m
    model.held_torque_body = cross(m, b_body)
    return nothing
end

"""Return the held rod torque (body frame); no force."""
function calcControlForceTorque(model::MagneticMomentumManagerModel, u, p, i::Int64, t::Float64)
    zero3 = SVector{3, Float64}(0.0, 0.0, 0.0)
    i == model.sat_idx || return zero3, zero3
    return zero3, model.held_torque_body
end

# No calcControlMassFlowRate method: torque rods consume no propellant, and
# the AbstractControlEffectorModel fallback (propulsive_maneuvers.jl) already
# returns 0.0 — defining one here with looser argument types is ambiguous
# against that fallback.
