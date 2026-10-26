const KS_IDENTITY3 = SMatrix{3,3,Float64}(I)
const KS_IDENTITY4 = SMatrix{4,4,Float64}(I)
const KS_IDENTITY10 = SMatrix{10,10,Float64}(I)
const KS_ZERO_MATRIX3 = zero(SMatrix{3,3,Float64,9})
const KS_ZERO_VECTOR3 = zero(SVector{3,Float64})

## Cartesian and perturbation derivatives

"""
Return the Jacobians of the Cartesian position and velocity with respect to
`(u, u′)`: `∂r/∂u`, `∂v/∂u`, `∂v/∂u′`, plus `|u|²`, `v` and `Λ(u)`.
"""
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

# ∂(Λ(u)ᵀ a)/∂u for a fixed acceleration `a`: Λ(u)ᵀ is linear in u.
@inline function ks_lifted_acceleration_jacobian(acceleration)
    ax, ay, az = acceleration
    return @SMatrix [
        ax  ay  az  0.0
        ay -ax  0.0 az
        az  0.0 -ax -ay
        0.0 az -ay ax
    ]
end

function ks_density_gradient(source, position, elapsed_time_s, params;
    step_m::Real=1.0)
    source isa Real && return @SVector [0.0, 0.0, 0.0]
    if applicable(source, Val(:gradient), position, elapsed_time_s)
        return SVector{3,Float64}(source(
            Val(:gradient), position, elapsed_time_s))
    end
    step = max(abs(Float64(step_m)), 0.1)
    center = ks_density_value(source, position, elapsed_time_s, params)
    return SVector{3,Float64}(ntuple(axis -> begin
        offset = SVector{3,Float64}(ntuple(
            index -> index == axis ? step : 0.0, 3))
        ks_density_difference(
            ks_density_value(source, position + offset, elapsed_time_s, params),
            center,
            ks_density_value(source, position - offset, elapsed_time_s, params),
            step)
    end, 3))
end

# Central difference, or a one-sided one when a neighbour lies past a density
# cutoff (zero density next to a positive one, as at a model's top altitude),
# where the central difference would report a spurious gradient.
@inline function ks_density_difference(plus, center, minus, step)
    if center > 0.0 && plus <= 0.0 && minus > 0.0
        return (center - minus) / step
    elseif center > 0.0 && minus <= 0.0 && plus > 0.0
        return (plus - center) / step
    end
    return (plus - minus) / (2.0 * step)
end

"""
    ks_density_value_gradient(source, position, elapsed_time_s, params; gradient_step_m=1)

Density (kg/m³) of a density source at a position and its Cartesian gradient:
the source's own `(Val(:gradient), position, elapsed_s)` method when it has one,
otherwise a finite difference with step `gradient_step_m` (m).
"""
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

# The perturbing acceleration (J2 and, with `use_drag`, drag) at a KS state, its
# Jacobians with respect to u, u′ and the exposed area, and the kinematic and
# atmospheric quantities they are built from.
function ks_acceleration_partials(state, params, area_m2;
    density_kg_m3=0.0,
    drag_coefficient::Real=0.0, mass_kg::Real=1.0,
    use_drag::Bool=false, density_gradient_step_m::Real=1.0,
    compute_area_partial::Val{ComputeArea}=Val(true)) where {ComputeArea}
    u = SVector{4,Float64}(state[1], state[2], state[3], state[4])
    u_prime = SVector{4,Float64}(state[5], state[6], state[7], state[8])
    position_u, velocity_u, velocity_u_prime, radius, velocity, lambda_u =
        ks_kinematics_jacobians(u, u_prime)
    position = lambda_u * u
    acceleration = ks_j2_acceleration_si(position, params)
    acceleration_position = ks_j2_acceleration_jacobian_si(position, params)
    acceleration_velocity = KS_ZERO_MATRIX3
    acceleration_area = KS_ZERO_VECTOR3
    density = 0.0
    density_gradient = KS_ZERO_VECTOR3
    relative_velocity = ks_relative_velocity(position, velocity, params)

    if use_drag
        coefficient = Float64(drag_coefficient)
        mass = Float64(mass_kg)
        mass > 0.0 || throw(ArgumentError("mass_kg must be positive."))
        density, density_gradient = ks_density_value_gradient(
            density_kg_m3, position, Float64(state[10]), params;
            gradient_step_m=density_gradient_step_m)
        rotation_cross = ks_rotation_cross_matrix(params)
        relative_speed = norm(relative_velocity)
        relative_speed > eps(Float64) || throw(ArgumentError(
            "Drag Jacobian is undefined at zero relative speed."))
        # ∂(|V| V)/∂V, with V = v - Ω × r.
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
        acceleration_area, radius, position, velocity, relative_velocity,
        position_u, velocity_u, velocity_u_prime, lambda_u,
        density, density_gradient)
end

## Continuous KS linearization

"""
    evaluate_ks_dynamics_and_jacobians(state, params, area_m2=0; density_kg_m3=0,
        drag_coefficient=0, mass_kg=1, use_drag=false, density_gradient_step_m=1,
        compute_input_jacobian=Val(true)) -> NamedTuple

The KS derivative ([`nonlinear_ks_rhs`](@ref)) at `state` and its analytic
Jacobians with respect to the state (`state`, 10×10) and the exposed area
(`input`, 10). `partials` holds the acceleration, kinematic and density
quantities they were built from.
"""
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
    h_ks = Float64(state[9])
    partials = ks_acceleration_partials(
        state, params, area_m2;
        density_kg_m3=density_kg_m3,
        drag_coefficient=drag_coefficient, mass_kg=mass_kg,
        use_drag=use_drag,
        density_gradient_step_m=density_gradient_step_m,
        compute_area_partial=compute_input_jacobian)
    radius = partials.radius
    lift = transpose(partials.lambda_u)   # Λ(u)ᵀ
    lifted_acceleration = lift * partials.acceleration
    # Rows of u″ = -(h/2) u + (|u|²/2) Λ(u)ᵀ a with respect to u and u′.
    Fqp = -0.5 * h_ks * KS_IDENTITY4 + 0.5 * (
        lifted_acceleration * transpose(2.0 * u) +
        radius * ks_lifted_acceleration_jacobian(partials.acceleration) +
        radius * lift * partials.acceleration_u)
    Fqq = 0.5 * radius * lift * partials.acceleration_u_prime
    # Row of h′ = -|u|² v·a.
    velocity_dot_acceleration = dot(partials.velocity, partials.acceleration)
    Fhp = -2.0 * velocity_dot_acceleration * u - radius * (
        transpose(partials.velocity_u) * partials.acceleration +
        transpose(partials.acceleration_u) * partials.velocity)
    Fhq = -radius * (
        transpose(partials.velocity_u_prime) * partials.acceleration +
        transpose(partials.acceleration_u_prime) * partials.velocity)
    F = MMatrix{10,10,Float64}(undef)
    fill!(F, 0.0)
    @inbounds for row in 1:4
        F[row, row + 4] = 1.0
        F[row + 4, 9] = -0.5 * u[row]
        F[9, row] = Fhp[row]
        F[9, row + 4] = Fhq[row]
        F[10, row] = 2.0 * u[row]     # t′ = |u|²
        for column in 1:4
            F[row + 4, column] = Fqp[row, column]
            F[row + 4, column + 4] = Fqq[row, column]
        end
    end
    derivative = ks_rhs_from_acceleration(
        u, u_prime, h_ks, radius, partials.velocity, partials.acceleration)
    ComputeInput || return (; derivative, state=SMatrix{10,10,Float64}(F), partials)
    Gamma = MVector{10,Float64}(undef)
    fill!(Gamma, 0.0)
    if use_drag
        lifted_area_derivative = lift * partials.acceleration_area
        @inbounds for row in 1:4
            Gamma[row + 4] = 0.5 * radius * lifted_area_derivative[row]
        end
        Gamma[9] = -radius * dot(partials.velocity, partials.acceleration_area)
    end
    return (; derivative, state=SMatrix{10,10,Float64}(F),
        input=SVector{10,Float64}(Gamma), partials)
end

"""
    ks_linear_implicit_midpoint_step(state, params, area_m2, delta_s; kwargs...) -> SVector{10}

Advance the KS state by one linearly implicit midpoint step,
`x + (I - Δs/2 F)⁻¹ Δs f(x)`, with the derivative `f` and its Jacobian `F`
evaluated once at `x` (no nonlinear iteration). `kwargs` are those of
[`evaluate_ks_dynamics_and_jacobians`](@ref).
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

"""
    ks_first_order_tangent_map(state, params, area_m2, delta_s; kwargs...) -> NamedTuple

First-order (forward Euler) discrete maps of a KS step of fictitious time
`delta_s` about `state`: `transition = I + Δs F` and `input_transition = Δs Γ`
for the exposed area, with the `evaluation` they come from (see
[`evaluate_ks_dynamics_and_jacobians`](@ref)).
"""
function ks_first_order_tangent_map(state::AbstractVector, params,
    area_m2::Real, delta_s::Real; kwargs...)
    step = Float64(delta_s)
    step > 0.0 || throw(ArgumentError("delta_s must be positive."))
    evaluation = evaluate_ks_dynamics_and_jacobians(
        state, params, area_m2; kwargs...)
    return (;
        transition=KS_IDENTITY10 + step * evaluation.state,
        input_transition=step * evaluation.input,
        evaluation,
    )
end
