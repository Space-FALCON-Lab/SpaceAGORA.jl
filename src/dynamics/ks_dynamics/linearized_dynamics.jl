const KS_IDENTITY3 = SMatrix{3,3,Float64}(I)
const KS_IDENTITY4 = SMatrix{4,4,Float64}(I)
const KS_IDENTITY10 = SMatrix{10,10,Float64}(I)
const KS_ZERO_MATRIX3 = zero(SMatrix{3,3,Float64,9})
const KS_ZERO_VECTOR3 = zero(SVector{3,Float64})

## Cartesian and perturbation derivatives

"""Return Cartesian position and velocity Jacobians with respect to `(u, u′)`."""
function ks_kinematics_jacobians(u, u_prime)
    lambda_u = ks_lambda_matrix(u)
    radius = dot(u, u)
    radius > eps(Float64) || throw(ArgumentError(
        "KS kinematics are undefined at the origin."))
    velocity = (2.0 / radius) * (lambda_u * u_prime)
    position_u = 2.0 * lambda_u
    velocity_u = (2.0 / radius) * ks_lambda_matrix(u_prime) -
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

@inline function ks_perturbation_matrix(u)
    u1, u2, u3, u4 = u
    return @SMatrix [
         u1  u2  u3
        -u2  u1  u4
        -u3 -u4  u1
         u4 -u3  u2
    ]
end

@inline function ks_perturbation_state_jacobian(acceleration)
    ax, ay, az = acceleration
    return @SMatrix [
        ax  ay  az  0.0
        ay -ax  0.0 az
        az  0.0 -ax -ay
        0.0 az -ay ax
    ]
end

@inline function ks_density_value(source, position, elapsed_time_s, params)
    source isa Real && return max(0.0, Float64(source))
    altitude_m = norm(position) - params.Re
    value = applicable(source, altitude_m, elapsed_time_s, position) ?
        source(altitude_m, elapsed_time_s, position) :
        source(altitude_m, elapsed_time_s)
    return max(0.0, Float64(value))
end

function ks_density_gradient(source, position, elapsed_time_s, params;
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
        (ks_density_value(source, position + offset, elapsed_time_s, params) -
         ks_density_value(source, position - offset, elapsed_time_s, params)) /
            (2.0 * step)
    end, 3))
end

"""Return atmospheric density and its Cartesian position gradient."""
function ks_density_value_gradient(source, position, elapsed_time_s, params;
    gradient_step_m::Real=1.0)
    position_vector = SVector{3,Float64}(position)
    density = ks_density_value(
        source, position_vector, Float64(elapsed_time_s), params)
    gradient = ks_density_gradient(
        source, position_vector, Float64(elapsed_time_s), params;
        step_m=gradient_step_m)
    return density, gradient
end

function ks_acceleration_partials(state, params, area_m2;
    density_kg_m3=0.0,
    drag_coefficient::Real=0.0, mass_kg::Real=1.0,
    use_drag::Bool=false, density_gradient_step_m::Real=1.0,
    compute_area_partial::Val{ComputeArea}=Val(true)) where {ComputeArea}
    u = SVector{4,Float64}(state[1], state[2], state[3], state[4])
    u_prime = SVector{4,Float64}(state[5], state[6], state[7], state[8])
    position_u, velocity_u, velocity_u_prime, radius, velocity, _ =
        ks_kinematics_jacobians(u, u_prime)
    position = ks_position(u)
    acceleration = ks_j2_acceleration_si(position, params)
    acceleration_position = ks_j2_acceleration_jacobian_si(position, params)
    acceleration_velocity = KS_ZERO_MATRIX3
    acceleration_area = KS_ZERO_VECTOR3

    if use_drag
        coefficient = Float64(drag_coefficient)
        mass = Float64(mass_kg)
        mass > 0.0 || throw(ArgumentError("mass_kg must be positive."))
        elapsed_time_s = Float64(state[10])
        density = ks_density_value(
            density_kg_m3, position, elapsed_time_s, params)
        density_gradient = ks_density_gradient(
            density_kg_m3, position, elapsed_time_s, params;
            step_m=density_gradient_step_m)
        rotation_cross = ks_rotation_cross_matrix(params)
        relative_velocity = velocity - rotation_cross * position
        relative_speed = norm(relative_velocity)
        relative_speed > eps(Float64) || throw(ArgumentError(
            "Drag Jacobian is undefined at zero relative speed."))
        velocity_product_jacobian = relative_speed * KS_IDENTITY3 +
            relative_velocity * transpose(relative_velocity) / relative_speed
        speed_weighted_velocity = relative_speed * relative_velocity
        drag_scale = 0.5 * coefficient * Float64(area_m2) / mass
        acceleration -= drag_scale * density * speed_weighted_velocity
        acceleration_position -= drag_scale * (
            speed_weighted_velocity * transpose(density_gradient) -
            density * velocity_product_jacobian * rotation_cross)
        acceleration_velocity -= drag_scale * density *
            velocity_product_jacobian
        if ComputeArea
            acceleration_area = -0.5 * coefficient / mass * density *
                speed_weighted_velocity
        end
    end

    acceleration_u = acceleration_position * position_u +
        acceleration_velocity * velocity_u
    acceleration_u_prime = acceleration_velocity * velocity_u_prime
    return (; acceleration, acceleration_u, acceleration_u_prime,
        acceleration_area, radius, velocity, velocity_u, velocity_u_prime)
end

## Continuous KS linearization

"""Evaluate the KS derivative and its analytical state and area Jacobians."""
function evaluate_ks_dynamics_and_jacobians(state::AbstractVector, params,
    area_m2::Real=0.0;
    density_kg_m3=0.0,
    drag_coefficient::Real=0.0, mass_kg::Real=1.0,
    use_drag::Bool=false, density_gradient_step_m::Real=1.0,
    compute_input_jacobian::Val{ComputeInput}=Val(true)) where {ComputeInput}
    length(state) == 10 || throw(ArgumentError(
        "KS state must have 10 components."))
    u = SVector{4,Float64}(state[1], state[2], state[3], state[4])
    u_prime = SVector{4,Float64}(state[5], state[6], state[7], state[8])
    partials = ks_acceleration_partials(
        state, params, area_m2;
        density_kg_m3=density_kg_m3,
        drag_coefficient=drag_coefficient, mass_kg=mass_kg,
        use_drag=use_drag,
        density_gradient_step_m=density_gradient_step_m,
        compute_area_partial=compute_input_jacobian)
    G = ks_perturbation_matrix(u)
    lifted_acceleration = G * partials.acceleration
    Fqp = -0.5 * Float64(state[9]) * KS_IDENTITY4 + 0.5 * (
        lifted_acceleration * transpose(2.0 * u) +
        partials.radius *
            ks_perturbation_state_jacobian(partials.acceleration) +
        partials.radius * G * partials.acceleration_u)
    Fqq = 0.5 * partials.radius * G * partials.acceleration_u_prime
    velocity_dot_acceleration = dot(
        partials.velocity, partials.acceleration)
    Fhp = -2.0 * velocity_dot_acceleration * u - partials.radius * (
        transpose(partials.velocity_u) * partials.acceleration +
        transpose(partials.acceleration_u) * partials.velocity)
    Fhq = -partials.radius * (
        transpose(partials.velocity_u_prime) * partials.acceleration +
        transpose(partials.acceleration_u_prime) * partials.velocity)
    F = MMatrix{10,10,Float64}(undef)
    fill!(F, 0.0)
    @inbounds for row in 1:4
        F[row, row + 4] = 1.0
        F[row + 4, 9] = -0.5 * u[row]
        F[9, row] = Fhp[row]
        F[9, row + 4] = Fhq[row]
        F[10, row] = 2.0 * u[row]
        for column in 1:4
            F[row + 4, column] = Fqp[row, column]
            F[row + 4, column + 4] = Fqq[row, column]
        end
    end
    derivative = MVector{10,Float64}(undef)
    @inbounds for row in 1:4
        derivative[row] = u_prime[row]
        derivative[row + 4] = -0.5 * Float64(state[9]) * u[row] +
            0.5 * partials.radius * lifted_acceleration[row]
    end
    derivative[9] = -partials.radius * velocity_dot_acceleration
    derivative[10] = partials.radius
    if !ComputeInput
        return (; derivative=SVector{10,Float64}(derivative),
            state=SMatrix{10,10,Float64}(F))
    end
    Gamma = MVector{10,Float64}(undef)
    fill!(Gamma, 0.0)
    if use_drag
        lifted_area_derivative = G * partials.acceleration_area
        @inbounds for row in 1:4
            Gamma[row + 4] = 0.5 * partials.radius *
                lifted_area_derivative[row]
        end
        Gamma[9] = -partials.radius *
            dot(partials.velocity, partials.acceleration_area)
    end
    return (; derivative=SVector{10,Float64}(derivative),
        state=SMatrix{10,10,Float64}(F),
        input=SVector{10,Float64}(Gamma))
end

"""Return continuous augmented-KS state and exposed-area Jacobians."""
function ks_rhs_jacobians(state::AbstractVector, params, area_m2::Real=0.0;
    kwargs...)
    result = evaluate_ks_dynamics_and_jacobians(state, params, area_m2; kwargs...)
    return (; state=result.state, input=result.input)
end

"""Return the analytical Jacobian of the continuous augmented-KS dynamics."""
ks_rhs_jacobian(state::AbstractVector, params, area_m2::Real=0.0; kwargs...) =
    evaluate_ks_dynamics_and_jacobians(state, params, area_m2;
        compute_input_jacobian=Val(false), kwargs...).state

"""
Advance the KS state with one linearly implicit midpoint correction.

The continuous dynamics and state Jacobian are evaluated once at the current
state. This propagation routine neither iterates a nonlinear residual nor
constructs discrete state or input maps.
"""
function ks_linear_implicit_midpoint_step(state::AbstractVector, params,
    area_m2::Real, delta_s::Real; kwargs...)
    x = SVector{10,Float64}(state)
    step = Float64(delta_s)
    step > 0.0 || throw(ArgumentError("delta_s must be positive."))
    evaluation = evaluate_ks_dynamics_and_jacobians(
        x, params, area_m2; compute_input_jacobian=Val(false), kwargs...)
    midpoint_matrix = KS_IDENTITY10 - 0.5 * step * evaluation.state
    increment = midpoint_matrix \ (step * evaluation.derivative)
    return x + increment
end

"""Return the first-order discrete state and area tangent maps for MPC."""
function ks_first_order_tangent_map(state::AbstractVector, params,
    area_m2::Real, delta_s::Real; kwargs...)
    step = Float64(delta_s)
    step > 0.0 || throw(ArgumentError("delta_s must be positive."))
    evaluation = evaluate_ks_dynamics_and_jacobians(
        state, params, area_m2; kwargs...)
    return (;
        transition=KS_IDENTITY10 + step * evaluation.state,
        input_transition=step * evaluation.input,
    )
end
