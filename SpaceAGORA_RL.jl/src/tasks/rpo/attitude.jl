function _quat_normalize(q::AbstractVector{<:Real})
    out = Vector{Float64}(q)
    n = norm(out)
    n > eps(Float64) || return [0.0, 0.0, 0.0, 1.0]
    return out ./ n
end

_quat_conjugate(q) = [-q[1], -q[2], -q[3], q[4]]

function _quat_multiply(a, b)
    av = view(a, 1:3)
    bv = view(b, 1:3)
    return vcat(a[4] .* bv .+ b[4] .* av .+ cross(av, bv),
                a[4] * b[4] - dot(av, bv))
end

function _axis_angle_quaternion(axis::Integer, angle::Real)
    q = zeros(4)
    q[Int(axis)] = sin(Float64(angle) / 2.0)
    q[4] = cos(Float64(angle) / 2.0)
    return q
end

function _quat_slerp(q0, q1, fraction::Real)
    a = _quat_normalize(q0)
    b = _quat_normalize(q1)
    d = dot(a, b)
    if d < 0.0
        b .*= -1.0
        d = -d
    end
    f = clamp(Float64(fraction), 0.0, 1.0)
    if d > 0.9995
        return _quat_normalize((1.0 - f) .* a .+ f .* b)
    end
    theta = acos(clamp(d, -1.0, 1.0))
    return (sin((1.0 - f) * theta) .* a .+ sin(f * theta) .* b) ./ sin(theta)
end

function _quat_rotation_matrix(q)
    x, y, z, w = _quat_normalize(q)
    return [
        1 - 2(y*y + z*z)  2(x*y - z*w)      2(x*z + y*w)
        2(x*y + z*w)      1 - 2(x*x + z*z)  2(y*z - x*w)
        2(x*z - y*w)      2(y*z + x*w)      1 - 2(x*x + y*y)
    ]
end

"""Return the minimum-rotation quaternion that maps one direction onto another."""
function _quat_between_directions(from, to)
    source = Vector{Float64}(from)
    target = Vector{Float64}(to)
    source_norm = norm(source)
    target_norm = norm(target)
    source_norm > eps(Float64) || return [0.0, 0.0, 0.0, 1.0]
    target_norm > eps(Float64) || return [0.0, 0.0, 0.0, 1.0]
    source ./= source_norm
    target ./= target_norm
    alignment = clamp(dot(source, target), -1.0, 1.0)
    alignment >= 1.0 - 1.0e-12 && return [0.0, 0.0, 0.0, 1.0]
    if alignment <= -1.0 + 1.0e-12
        helper = abs(source[1]) <= abs(source[2]) &&
                 abs(source[1]) <= abs(source[3]) ?
            [1.0, 0.0, 0.0] :
            (abs(source[2]) <= abs(source[3]) ?
                [0.0, 1.0, 0.0] : [0.0, 0.0, 1.0])
        axis = cross(source, helper)
        axis ./= norm(axis)
        return vcat(axis, 0.0)
    end
    return _quat_normalize(vcat(cross(source, target), 1.0 + alignment))
end

"""
Point the nearest available body thruster along the requested RTN acceleration.

For the configured equal-Isp, orthogonal thruster pairs this minimizes
instantaneous translation impulse while requiring the smallest attitude slew.
The smooth endpoint weight keeps the prescribed initial and final attitudes.
"""
function _rpo_fuel_optimal_pointing_attitude(
    attitude_rtn_to_body,
    baseline_attitude_rtn_to_body,
    command_rtn,
    progress::Real,
    config::RPOHyPRRLConfig,
)
    baseline = _quat_normalize(baseline_attitude_rtn_to_body)
    fraction = clamp(Float64(progress), 0.0, 1.0)
    (fraction <= 0.0 || fraction >= 1.0) && return baseline

    command = Vector{Float64}(command_rtn)
    command_norm = norm(command)
    command_norm > 1.0e-12 || return baseline
    current = _quat_normalize(attitude_rtn_to_body)
    force_body = _quat_rotation_matrix(current) * (command ./ command_norm)

    best_index = 0
    best_alignment = -Inf
    @inbounds for index in axes(config.thruster_directions_body, 2)
        config.thruster_max_thrust_n[index] > 0.0 || continue
        direction = view(config.thruster_directions_body, :, index)
        direction_norm = norm(direction)
        direction_norm > 1.0e-12 || continue
        alignment = dot(force_body, direction) / direction_norm
        if alignment > best_alignment
            best_alignment = alignment
            best_index = index
        end
    end
    best_index > 0 || return baseline

    target_direction = view(config.thruster_directions_body, :, best_index)
    delta = _quat_between_directions(force_body, target_direction)
    aligned = _quat_normalize(_quat_multiply(delta, current))
    pointing_weight = sinpi(fraction)^2
    return _quat_slerp(baseline, aligned, pointing_weight)
end

function _attitude_reference(progress, quaternions, sample_progress)
    n = length(sample_progress)
    result = zeros(4, n)
    for (column, value) in pairs(sample_progress)
        upper = searchsortedfirst(progress, value)
        if upper <= 1
            result[:, column] .= quaternions[:, 1]
        elseif upper > length(progress)
            result[:, column] .= quaternions[:, end]
        else
            lower = upper - 1
            span = max(progress[upper] - progress[lower], eps(Float64))
            result[:, column] .= _quat_slerp(
                view(quaternions, :, lower), view(quaternions, :, upper),
                (value - progress[lower]) / span,
            )
        end
    end
    return result
end

function _integrate_attitude(q, omega, desired_q, dt, config::RPOHyPRRLConfig)
    q_error = _quat_multiply(desired_q, _quat_conjugate(q))
    sign_scalar = q_error[4] < 0.0 ? -1.0 : 1.0
    torque = config.reaction_wheel_kp .* (2.0 * sign_scalar .* q_error[1:3]) .-
             config.reaction_wheel_kd .* omega
    torque = clamp.(torque, -config.reaction_wheel_max_torque_nm,
                    config.reaction_wheel_max_torque_nm)
    inertia = config.reaction_wheel_inertia_kgm2
    omega_dot = (torque .- cross(omega, inertia .* omega)) ./ inertia
    next_omega = omega .+ Float64(dt) .* omega_dot
    q_dot = 0.5 .* _quat_multiply(vcat(next_omega, 0.0), q)
    next_q = _quat_normalize(q .+ Float64(dt) .* q_dot)
    return next_q, next_omega, torque
end
