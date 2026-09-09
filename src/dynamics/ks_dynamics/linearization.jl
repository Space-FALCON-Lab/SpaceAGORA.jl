"""Return Cartesian-position and velocity Jacobians with respect to `(u, u′)`."""
function ks_kinematics_jacobians(u, u_prime)
    lambda_u = _ks_lambda(u)
    radius = dot(u, u)
    radius > eps(Float64) || throw(ArgumentError("KS kinematics are undefined at the origin."))
    velocity = (2.0 / radius) * (lambda_u * u_prime)
    position_u = 2.0 * lambda_u
    velocity_u = (2.0 / radius) * _ks_lambda(u_prime) -
        (4.0 / radius^2) * ((lambda_u * u_prime) * transpose(u))
    velocity_u_prime = (2.0 / radius) * lambda_u
    return position_u, velocity_u, velocity_u_prime, radius, velocity, lambda_u
end

"""Return the analytic Cartesian Jacobian of the J2 perturbing acceleration."""
function ks_j2_acceleration_jacobian_si(rvec, params)
    x, y, z = rvec
    radius2 = dot(rvec, rvec)
    radius = sqrt(radius2)
    radius <= eps(Float64) && return @SMatrix zeros(3, 3)
    radius4 = radius2^2
    radius5 = radius^5
    radius7 = radius^7

    coefficient = 1.5 * params.J2 * params.μ * params.Re^2
    z2_r2 = z^2 / radius2
    xy_factor = 5.0 * z2_r2 - 1.0
    z_factor = 5.0 * z2_r2 - 3.0
    shape = @SVector [x * xy_factor, y * xy_factor, z * z_factor]

    d_factor_dx = -10.0 * z^2 * x / radius4
    d_factor_dy = -10.0 * z^2 * y / radius4
    d_factor_dz = 10.0 * z * (x^2 + y^2) / radius4
    column_x = @SVector [xy_factor + x * d_factor_dx, y * d_factor_dx, z * d_factor_dx]
    column_y = @SVector [x * d_factor_dy, xy_factor + y * d_factor_dy, z * d_factor_dy]
    column_z = @SVector [x * d_factor_dz, y * d_factor_dz, z_factor + z * d_factor_dz]
    shape_jacobian = SMatrix{3, 3, Float64}(hcat(column_x, column_y, column_z))

    radial_gradient = (-5.0 * coefficient / radius7) * rvec
    return (coefficient / radius5) * shape_jacobian + shape * transpose(radial_gradient)
end

function _ks_central_jacobian(f, state::AbstractVector; relative_step::Real=eps(Float64)^(1 / 5))
    x = Float64.(collect(state))
    output_length = length(f(x))
    jacobian = Matrix{Float64}(undef, output_length, length(x))
    for column in eachindex(x)
        step = max(abs(x[column]), 1.0) * Float64(relative_step)
        x_plus = copy(x)
        x_minus = copy(x)
        x_plus[column] += step
        x_minus[column] -= step
        coarse = (f(x_plus) .- f(x_minus)) ./ (2.0 * step)
        half_step = 0.5 * step
        x_plus[column] = x[column] + half_step
        x_minus[column] = x[column] - half_step
        fine = (f(x_plus) .- f(x_minus)) ./ (2.0 * half_step)
        jacobian[:, column] .= (4.0 .* fine .- coarse) ./ 3.0
    end
    return jacobian
end

"""
    ks_rhs_jacobian(state, params, area_m2=0; kwargs...)

Numerically linearize the complete reusable KS right-hand side, including the
energy and physical-time states. Keyword arguments are forwarded to `ks_rhs`.
"""
function ks_rhs_jacobian(state::AbstractVector, params, area_m2::Real=0.0;
    relative_step::Real=eps(Float64)^(1 / 5), kwargs...)
    return _ks_central_jacobian(
        x -> ks_rhs(x, params, area_m2; kwargs...),
        state;
        relative_step=relative_step,
    )
end

"""
    ks_step_jacobian(state, params, area_m2, delta_s; kwargs...)

Numerically linearize one complete KS RK4 step. This discrete state-transition
Jacobian is independent of MPC output definitions and is reusable by estimation,
sensitivity analysis, and other guidance or control algorithms.
"""
function ks_step_jacobian(state::AbstractVector, params, area_m2::Real, delta_s::Real;
    relative_step::Real=eps(Float64)^(1 / 5), kwargs...)
    return _ks_central_jacobian(
        x -> ks_rk4_step(x, params, area_m2, delta_s; kwargs...),
        state;
        relative_step=relative_step,
    )
end
