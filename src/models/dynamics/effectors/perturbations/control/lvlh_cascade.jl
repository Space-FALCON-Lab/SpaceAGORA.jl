"""
    LVLHCascadeAttitudeControlModel <: AbstractForceTorqueModel

Closed-loop LVLH-pointing attitude controller as a dynamic effector: an outer
proportional attitude loop commands a rate-limited body rate toward a fixed
LVLH-relative attitude setpoint, and an inner proportional rate loop converts
the rate error into a body torque with a per-axis saturation. This cascade
(rate-command) architecture is the common small-satellite ADCS structure
(e.g. Blue Canyon XACT-class flight software).

    omega_cmd = clamp.(-k_out .* theta_err, -w_max, w_max)
    tau       = clamp.(k_rate .* (omega_cmd - omega_rel), -tau_max, tau_max)

`theta_err` is the small-angle rotation vector (body coordinates) from the
commanded attitude to the current attitude, extracted from
`R_err = R_bl * R_cmd'`; `omega_rel` is the body rate relative to the rotating
LVLH frame. The LVLH triad is +Z nadir, +Y = nadir x velocity (negative orbit
normal), +X completing (near along-track), and `q_cmd_lb` is the commanded
LVLH-to-body quaternion (scalar-last, identity = LVLH-aligned).

The torque is applied directly to the body (zero force): the actuator stack
(wheel allocation, momentum management) is abstracted into `tau_max`. The
control law is evaluated continuously at the integrator rate — a
flight-software discrete update (a few Hz) is well inside the closed-loop
bandwidth this law can express, so the continuous approximation is benign for
tracking studies. Signs follow the engine's own conventions
([`quaternion_derivative`](@ref) gives `theta_dot = +omega_body` near
identity), which the probe suite pins with an end-to-end convergence test.
Returns zero wrench when attitude state is unavailable.
"""
struct LVLHCascadeAttitudeControlModel <: AbstractForceTorqueModel
    q_cmd_lb::SVector{4, Float64}
    k_out::SVector{3, Float64}
    w_max::Float64
    k_rate::SVector{3, Float64}
    tau_max::Float64
    tau_ff::SVector{3, Float64}

    function LVLHCascadeAttitudeControlModel(;
        q_cmd_lb::AbstractVector{<:Real}=SVector{4, Float64}(0.0, 0.0, 0.0, 1.0),
        k_out::AbstractVector{<:Real},
        w_max::Real,
        k_rate::AbstractVector{<:Real},
        tau_max::Real,
        tau_ff::AbstractVector{<:Real}=SVector{3, Float64}(0.0, 0.0, 0.0),
    )
        q = SVector{4, Float64}(q_cmd_lb...)
        qn = norm(q)
        (isfinite(qn) && qn > 0) ||
            throw(ArgumentError("q_cmd_lb must be a finite, nonzero quaternion"))
        ko = SVector{3, Float64}(k_out...)
        kr = SVector{3, Float64}(k_rate...)
        all(isfinite, ko) && all(>=(0.0), ko) ||
            throw(ArgumentError("k_out must be finite and nonnegative [1/s]"))
        all(isfinite, kr) && all(>=(0.0), kr) ||
            throw(ArgumentError("k_rate must be finite and nonnegative [N m s]"))
        (isfinite(w_max) && w_max > 0) ||
            throw(ArgumentError("w_max must be finite and positive [rad/s]"))
        (isfinite(tau_max) && tau_max > 0) ||
            throw(ArgumentError("tau_max must be finite and positive [N m]"))
        tf = SVector{3, Float64}(tau_ff...)
        all(isfinite, tf) ||
            throw(ArgumentError("tau_ff must be finite [N m]"))
        return new(q / qn, ko, Float64(w_max), kr, Float64(tau_max), tf)
    end
end

"""
    _lvlh_cascade_torque(model, pos_ii, vel_ii, q_ib, w_body)

Body-frame control torque of [`LVLHCascadeAttitudeControlModel`](@ref) for the
given inertial position/velocity and attitude state. Returns zero for
degenerate orbits (near-zero radius or angular momentum).
"""
@inline function _lvlh_cascade_torque(
    model::LVLHCascadeAttitudeControlModel,
    pos_ii::SVector{3, Float64},
    vel_ii::SVector{3, Float64},
    q_ib::SVector{4, Float64},
    w_body::SVector{3, Float64},
)::SVector{3, Float64}
    r2 = dot(pos_ii, pos_ii)
    h_ii = cross(pos_ii, vel_ii)
    if !(r2 > 0.0) || !(dot(h_ii, h_ii) > 0.0)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    z_l = -pos_ii / sqrt(r2)
    y_l = cross(z_l, vel_ii)
    y_l = y_l / norm(y_l)
    x_l = cross(y_l, z_l)
    R_li = SMatrix{3, 3, Float64, 9}(
        x_l[1], y_l[1], z_l[1],
        x_l[2], y_l[2], z_l[2],
        x_l[3], y_l[3], z_l[3],
    )                                    # inertial -> LVLH (rows = triad)
    R_bi = rot(q_ib)                     # inertial -> body
    R_bl = R_bi * R_li'                  # LVLH -> body
    R_e = R_bl * rot(model.q_cmd_lb)'    # commanded -> current, body coords
    # Rotation-log error extraction. The plain vee-map (0.5*(R_e[2,3]-R_e[3,2]),
    # ...) returns sin(theta)*axis, which vanishes at a half-turn — an
    # undesired zero-command equilibrium at 180 deg and weak authority near it
    # (Codex review). The quaternion log 2*atan2(|qv|, qs)*qv/|qv| equals the
    # vee-map to second order at small angles (2*qs*qv == sin(theta)*axis)
    # but grows monotonically to pi*axis at the antipode. Index order/signs
    # follow the engine conventions pinned empirically in the probe suite:
    # rotating the body by +theta yields rot(q) ≈ I - S(theta); the opposite
    # sign turns the attitude loop into positive feedback.
    tr_e = R_e[1, 1] + R_e[2, 2] + R_e[3, 3]
    q_e = if tr_e > 0.0
        s = sqrt(tr_e + 1.0) * 2
        SVector{4, Float64}((R_e[2, 3] - R_e[3, 2]) / s, (R_e[3, 1] - R_e[1, 3]) / s,
                            (R_e[1, 2] - R_e[2, 1]) / s, s / 4)
    elseif R_e[1, 1] > R_e[2, 2] && R_e[1, 1] > R_e[3, 3]
        s = sqrt(1.0 + R_e[1, 1] - R_e[2, 2] - R_e[3, 3]) * 2
        SVector{4, Float64}(s / 4, (R_e[1, 2] + R_e[2, 1]) / s,
                            (R_e[1, 3] + R_e[3, 1]) / s, (R_e[2, 3] - R_e[3, 2]) / s)
    elseif R_e[2, 2] > R_e[3, 3]
        s = sqrt(1.0 + R_e[2, 2] - R_e[1, 1] - R_e[3, 3]) * 2
        SVector{4, Float64}((R_e[1, 2] + R_e[2, 1]) / s, s / 4,
                            (R_e[2, 3] + R_e[3, 2]) / s, (R_e[3, 1] - R_e[1, 3]) / s)
    else
        s = sqrt(1.0 + R_e[3, 3] - R_e[1, 1] - R_e[2, 2]) * 2
        SVector{4, Float64}((R_e[1, 3] + R_e[3, 1]) / s, (R_e[2, 3] + R_e[3, 2]) / s,
                            s / 4, (R_e[1, 2] - R_e[2, 1]) / s)
    end
    q_e = q_e[4] < 0.0 ? -q_e : q_e      # shortest rotation (qs >= 0)
    qv = SVector{3, Float64}(q_e[1], q_e[2], q_e[3])
    sv = norm(qv)
    theta_err = sv > 1e-12 ? (2.0 * atan(sv, q_e[4]) / sv) * qv : 2.0 * qv
    # Body-frame LVLH feed-forward: with the body-rate quaternion kinematics
    # (see quaternion_derivative), the state that tracks the LVLH frame is
    # its rotation rate expressed in BODY coordinates — pinned by the
    # torque-free tracking and closed-loop convergence probes.
    w_lvlh_body = R_bi * (h_ii / r2)
    w_rel = w_body - w_lvlh_body
    w_cmd = clamp.(-model.k_out .* theta_err, -model.w_max, model.w_max)
    # tau_ff: constant body-frame feedforward torque, the stateless equivalent
    # of converged integral action against a quasi-steady disturbance (a pure
    # PD cascade holds a steady offset error of |tau_dist| / (k_rate*k_out)
    # under a constant torque; the feedforward nulls it). True discrete-rate
    # integral action belongs in a stateful control effector - a continuous
    # wrench cannot hold accumulator state correctly under adaptive stepping.
    return clamp.(model.k_rate .* (w_cmd - w_rel) + model.tau_ff, -model.tau_max, model.tau_max)
end

@inline environment_requirements(::LVLHCascadeAttitudeControlModel) = EffectorEnvironmentRequirements()

@inline function wrench(
    model::LVLHCascadeAttitudeControlModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    zero_wrench = (SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))
    if x.q_ib === nothing || x.ω_body === nothing
        return zero_wrench
    end
    torque = _lvlh_cascade_torque(model, x.pos_ii, x.vel_ii, x.q_ib, x.ω_body)
    return SVector{3, Float64}(0.0, 0.0, 0.0), torque
end
