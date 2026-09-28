"""Return Cartesian position and velocity Jacobians with respect to `(u, u′)`."""
function ks_kinematics_jacobians(u, u_prime)
    lambda_u = _ks_lambda(u)
    radius = dot(u, u)
    radius > eps(Float64) || throw(ArgumentError(
        "KS kinematics are undefined at the origin."))
    velocity = (2.0 / radius) * (lambda_u * u_prime)
    position_u = 2.0 * lambda_u
    velocity_u = (2.0 / radius) * _ks_lambda(u_prime) -
        (4.0 / radius^2) * ((lambda_u * u_prime) * transpose(u))
    velocity_u_prime = (2.0 / radius) * lambda_u
    return position_u, velocity_u, velocity_u_prime, radius, velocity, lambda_u
end

"""Return the analytic Cartesian Jacobian of the J2 acceleration."""
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
    column_x = @SVector [xy_factor + x * d_factor_dx,
        y * d_factor_dx, z * d_factor_dx]
    column_y = @SVector [x * d_factor_dy,
        xy_factor + y * d_factor_dy, z * d_factor_dy]
    column_z = @SVector [x * d_factor_dz,
        y * d_factor_dz, z_factor + z * d_factor_dz]
    shape_jacobian = SMatrix{3,3,Float64}(
        hcat(column_x, column_y, column_z))
    radial_gradient = (-5.0 * coefficient / radius7) * rvec
    return (coefficient / radius5) * shape_jacobian +
        shape * transpose(radial_gradient)
end

@inline function _ks_perturbation_matrix(u)
    u1, u2, u3, u4 = u
    return @SMatrix [
         u1  u2  u3
        -u2  u1  u4
        -u3 -u4  u1
         u4 -u3  u2
    ]
end

@inline function _ks_perturbation_state_jacobian(acceleration)
    ax, ay, az = acceleration
    return @SMatrix [
        ax  ay  az  0.0
        ay -ax  0.0 az
        az  0.0 -ax -ay
        0.0 az -ay ax
    ]
end

@inline function _ks_density_value(source, position, elapsed_time_s, params)
    source isa Real && return max(0.0, Float64(source))
    altitude_m = norm(position) - params.Re
    value = applicable(source, altitude_m, elapsed_time_s, position) ?
        source(altitude_m, elapsed_time_s, position) :
        source(altitude_m, elapsed_time_s)
    return max(0.0, Float64(value))
end

function _ks_density_gradient(source, position, elapsed_time_s, params;
    step_m::Real=1.0)
    source isa Real && return @SVector [0.0, 0.0, 0.0]
    if applicable(source, Val(:gradient), position, elapsed_time_s)
        return SVector{3,Float64}(source(
            Val(:gradient), position, elapsed_time_s))
    end
    step = max(abs(Float64(step_m)), 0.1)
    return SVector{3,Float64}(ntuple(axis -> begin
        offset = SVector{3,Float64}(ntuple(
            index -> index == axis ? step : 0.0, 3))
        (_ks_density_value(source, position + offset, elapsed_time_s, params) -
         _ks_density_value(source, position - offset, elapsed_time_s, params)) /
            (2.0 * step)
    end, 3))
end

"""Return atmospheric density and its Cartesian position gradient."""
function ks_density_value_gradient(source, position, elapsed_time_s, params;
    gradient_step_m::Real=1.0)
    position_vector = SVector{3,Float64}(position)
    density = _ks_density_value(
        source, position_vector, Float64(elapsed_time_s), params)
    gradient = _ks_density_gradient(
        source, position_vector, Float64(elapsed_time_s), params;
        step_m=gradient_step_m)
    return density, gradient
end

function _ks_acceleration_partials(state, params, area_m2;
    density_kg_m3=0.0,
    drag_coefficient::Real=0.0, mass_kg::Real=1.0,
    use_drag::Bool=false, density_gradient_step_m::Real=1.0)
    u = SVector{4,Float64}(state[1:4])
    u_prime = SVector{4,Float64}(state[5:8])
    position_u, velocity_u, velocity_u_prime, radius, velocity, _ =
        ks_kinematics_jacobians(u, u_prime)
    position = ks_position(u)
    acceleration = ks_j2_acceleration_si(position, params)
    acceleration_position = Matrix(
        ks_j2_acceleration_jacobian_si(position, params))
    acceleration_velocity = zeros(3, 3)
    acceleration_area = @SVector [0.0, 0.0, 0.0]

    if use_drag
        coefficient = Float64(drag_coefficient)
        mass = Float64(mass_kg)
        mass > 0.0 || throw(ArgumentError("mass_kg must be positive."))
        elapsed_time_s = Float64(state[10])
        density = _ks_density_value(
            density_kg_m3, position, elapsed_time_s, params)
        density_gradient = _ks_density_gradient(
            density_kg_m3, position, elapsed_time_s, params;
            step_m=density_gradient_step_m)
        rotation_cross = ks_rotation_cross_matrix(params)
        relative_velocity = velocity - rotation_cross * position
        relative_speed = norm(relative_velocity)
        relative_speed > eps(Float64) || throw(ArgumentError(
            "Drag Jacobian is undefined at zero relative speed."))
        velocity_product_jacobian = relative_speed *
            Matrix{Float64}(I, 3, 3) +
            relative_velocity * transpose(relative_velocity) / relative_speed
        speed_weighted_velocity = relative_speed * relative_velocity
        drag_scale = 0.5 * coefficient * Float64(area_m2) / mass
        acceleration -= drag_scale * density * speed_weighted_velocity
        acceleration_position -= drag_scale * (
            speed_weighted_velocity * transpose(density_gradient) -
            density * velocity_product_jacobian * rotation_cross)
        acceleration_velocity -= drag_scale * density *
            velocity_product_jacobian
        acceleration_area = -0.5 * coefficient / mass * density *
            speed_weighted_velocity
    end

    acceleration_u = acceleration_position * position_u +
        acceleration_velocity * velocity_u
    acceleration_u_prime = acceleration_velocity * velocity_u_prime
    return (; acceleration, acceleration_u, acceleration_u_prime,
        acceleration_area, radius, velocity, velocity_u, velocity_u_prime)
end

"""Return continuous augmented-KS state and exposed-area Jacobians."""
function ks_rhs_jacobians(state::AbstractVector, params, area_m2::Real=0.0;
    density_kg_m3=0.0,
    drag_coefficient::Real=0.0, mass_kg::Real=1.0,
    use_drag::Bool=false, density_gradient_step_m::Real=1.0)
    length(state) == 10 || throw(ArgumentError(
        "KS state must have 10 components."))
    u = SVector{4,Float64}(state[1:4])
    partials = _ks_acceleration_partials(
        state, params, area_m2;
        density_kg_m3=density_kg_m3,
        drag_coefficient=drag_coefficient, mass_kg=mass_kg,
        use_drag=use_drag,
        density_gradient_step_m=density_gradient_step_m)
    G = _ks_perturbation_matrix(u)
    lifted_acceleration = G * partials.acceleration
    F = zeros(10, 10)
    F[1:4, 5:8] .= Matrix{Float64}(I, 4, 4)
    F[5:8, 1:4] .= -0.5 * Float64(state[9]) *
        Matrix{Float64}(I, 4, 4) + 0.5 * (
            lifted_acceleration * transpose(2.0 * u) +
            partials.radius *
                _ks_perturbation_state_jacobian(partials.acceleration) +
            partials.radius * G * partials.acceleration_u)
    F[5:8, 5:8] .= 0.5 * partials.radius * G *
        partials.acceleration_u_prime
    F[5:8, 9] .= -0.5 * u
    velocity_dot_acceleration = dot(
        partials.velocity, partials.acceleration)
    F[9, 1:4] .= vec(
        -2.0 * velocity_dot_acceleration * transpose(u) -
        partials.radius * (
            transpose(partials.acceleration) * partials.velocity_u +
            transpose(partials.velocity) * partials.acceleration_u))
    F[9, 5:8] .= vec(-partials.radius * (
        transpose(partials.acceleration) * partials.velocity_u_prime +
        transpose(partials.velocity) * partials.acceleration_u_prime))
    F[10, 1:4] .= 2.0 * u
    Gamma = zeros(10, 1)
    if use_drag
        lifted_area_derivative = G * partials.acceleration_area
        Gamma[5:8, 1] .= 0.5 * partials.radius * lifted_area_derivative
        Gamma[9, 1] = -partials.radius *
            dot(partials.velocity, partials.acceleration_area)
    end
    return (; state=F, input=Gamma)
end

"""Return the analytical Jacobian of the continuous augmented-KS dynamics."""
ks_rhs_jacobian(state::AbstractVector, params, area_m2::Real=0.0; kwargs...) =
    ks_rhs_jacobians(state, params, area_m2; kwargs...).state

function _ks_implicit_midpoint_solution(state, params, area_m2, delta_s;
    maximum_iterations::Integer=12,
    nonlinear_tolerance::Real=2.0e-13, kwargs...)
    x = Float64.(collect(state))
    step = Float64(delta_s)
    step > 0.0 || throw(ArgumentError("delta_s must be positive."))
    next_state = x + step * ks_rhs(x, params, area_m2; kwargs...)
    identity10 = Matrix{Float64}(I, 10, 10)
    for iteration in 1:Int(maximum_iterations)
        midpoint = 0.5 * (x + next_state)
        residual = next_state - x -
            step * ks_rhs(midpoint, params, area_m2; kwargs...)
        jacobians = ks_rhs_jacobians(
            midpoint, params, area_m2; kwargs...)
        factorization = lu!(identity10 - 0.5 * step * jacobians.state)
        correction = factorization \ (-residual)
        next_state += correction
        scale = max.(max.(abs.(x), abs.(next_state)), 1.0)
        if maximum(abs.(correction) ./ scale) <= nonlinear_tolerance
            final_midpoint = 0.5 * (x + next_state)
            final_jacobians = ks_rhs_jacobians(
                final_midpoint, params, area_m2; kwargs...)
            final_factorization = lu!(identity10 -
                0.5 * step * final_jacobians.state)
            return (; state=next_state, jacobians=final_jacobians,
                factorization=final_factorization, iterations=iteration)
        end
    end
    throw(ErrorException(
        "KS implicit-midpoint solve did not converge in $(maximum_iterations) iterations."))
end

"""Advance the nonlinear augmented-KS state with implicit midpoint."""
function ks_implicit_midpoint_step(state::AbstractVector, params,
    area_m2::Real, delta_s::Real; kwargs...)
    return _ks_implicit_midpoint_solution(
        state, params, area_m2, delta_s; kwargs...).state
end

"""Return one implicit-midpoint step and its state and input tangent maps."""
function ks_implicit_midpoint_linearization(state::AbstractVector, params,
    area_m2::Real, delta_s::Real; kwargs...)
    result = _ks_implicit_midpoint_solution(
        state, params, area_m2, delta_s; kwargs...)
    step = Float64(delta_s)
    identity10 = Matrix{Float64}(I, 10, 10)
    solution = result.factorization \ hcat(
        identity10 + 0.5 * step * result.jacobians.state,
        step .* result.jacobians.input)
    return (; state=result.state,
        transition=solution[:, 1:10],
        input_transition=solution[:, 11:11],
        iterations=result.iterations)
end

"""Return the discrete implicit-midpoint state-transition Jacobian."""
function ks_step_jacobian(state::AbstractVector, params, area_m2::Real,
    delta_s::Real; kwargs...)
    return ks_implicit_midpoint_linearization(
        state, params, area_m2, delta_s; kwargs...).transition
end
