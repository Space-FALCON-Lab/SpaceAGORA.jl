## Propagation parameters and energy convention

"""Planet constants required by the reusable KS orbital propagator."""
Base.@kwdef struct KSPropagationParams
    Re::Float64
    μ::Float64
    J2::Float64 = 0.0
    Ω::Float64 = 0.0
end

"""Return the KS energy parameter `h_KS = -specific_energy`."""
ks_energy_parameter(specific_energy::Real) = -Float64(specific_energy)

"""Recover specific orbital energy from the KS energy parameter."""
specific_energy_from_ks(h_ks::Real) = -Float64(h_ks)

"""Return the cross-product matrix for the planet rotation vector `[0, 0, Ω]`."""
function ks_rotation_cross_matrix(params)
    omega = Float64(params.Ω)
    return @SMatrix [0.0 -omega 0.0; omega 0.0 0.0; 0.0 0.0 0.0]
end

## Perturbing accelerations

"""Return the J2 perturbing acceleration in SI units."""
function ks_j2_acceleration_si(rvec, params)
    r = norm(rvec)
    r <= eps(Float64) && return @SVector [0.0, 0.0, 0.0]
    return ks_j2_acceleration_with_radius(rvec, r, params)
end

@inline function ks_j2_acceleration_with_radius(rvec, r, params)
    scale = 3.0 * params.J2 * params.μ * params.Re^2 / (2.0 * r^5)
    z2_r2 = rvec[3]^2 / r^2
    return scale * @SVector [
        rvec[1] * (5.0 * z2_r2 - 1.0),
        rvec[2] * (5.0 * z2_r2 - 1.0),
        rvec[3] * (5.0 * z2_r2 - 3.0),
    ]
end

"""Return atmospheric drag acceleration in SI units."""
function ks_drag_acceleration_si(
    rvec,
    vvec,
    params,
    area_m2::Real;
    density_kg_m3::Real,
    drag_coefficient::Real,
    mass_kg::Real,
)
    density = max(0.0, Float64(density_kg_m3))
    mass = Float64(mass_kg)
    mass > 0.0 || throw(ArgumentError("mass_kg must be positive."))
    # Evaluate Ω × r directly to avoid constructing a skew matrix in the RHS.
    omega = Float64(params.Ω)
    relative_velocity = @SVector [
        vvec[1] + omega * rvec[2],
        vvec[2] - omega * rvec[1],
        vvec[3],
    ]
    speed = norm(relative_velocity)
    speed <= eps(Float64) && return @SVector [0.0, 0.0, 0.0]
    return -0.5 * density * Float64(drag_coefficient) * Float64(area_m2) /
        mass * speed * relative_velocity
end

## State conversion

"""
    cartesian_to_ks_state(position_ii_m, velocity_ii_m, params; elapsed_time_s=0)

Create the 10-component KS state `[u, u′, h, t]` from an inertial Cartesian
state. Distances, velocity, gravitational parameter, and time use SI units.
"""
function cartesian_to_ks_state(position_ii_m, velocity_ii_m, params; elapsed_time_s::Real=0.0)
    r = SVector{3, Float64}(position_ii_m)
    v = SVector{3, Float64}(velocity_ii_m)
    u = SVector{4, Float64}(cartesian_position_to_ks_coordinate(position_ii_m)...)
    u_prime = SVector{4, Float64}(cartesian_velocity_to_ks_derivative(velocity_ii_m, u)...)
    energy = 0.5 * dot(v, v) - params.μ / norm(r)
    h_ks = ks_energy_parameter(energy)
    return MVector{10,Float64}(vcat(u, u_prime, h_ks, Float64(elapsed_time_s)))
end

"""Convert a 10-component KS state back to inertial Cartesian state."""
function ks_state_to_cartesian(state)
    u = SVector{4, Float64}(state[1], state[2], state[3], state[4])
    u_prime = SVector{4, Float64}(state[5], state[6], state[7], state[8])
    return (
        position_ii_m=ks_position(u),
        velocity_ii_m=ks_velocity(u, u_prime),
        h_ks=Float64(state[9]),
        energy_parameter=Float64(state[9]),
        specific_energy_j_kg=specific_energy_from_ks(state[9]),
        elapsed_time_s=Float64(state[10]),
    )
end

## Nonlinear KS dynamics

"""Evaluate the KS fictitious-time right-hand side."""
function nonlinear_ks_rhs(state::AbstractVector, params,
    area_m2::Real=0.0;
    density_kg_m3=0.0,
    drag_coefficient::Real=0.0,
    mass_kg::Real=1.0,
    use_drag::Bool=false,
)
    length(state) == 10 || throw(ArgumentError("KS state must have 10 components."))
    u = SVector{4, Float64}(state[1], state[2], state[3], state[4])
    u_prime = SVector{4, Float64}(state[5], state[6], state[7], state[8])
    h_ks = Float64(state[9])
    rvec = ks_position(u)
    r = dot(u, u)
    vvec = ks_velocity_with_radius(u, u_prime, r)
    acceleration = ks_j2_acceleration_with_radius(rvec, r, params)
    if use_drag
        density_value = density_kg_m3 isa Real ? density_kg_m3 :
            applicable(density_kg_m3, norm(rvec) - params.Re,
                Float64(state[10]), rvec) ?
                density_kg_m3(norm(rvec) - params.Re, Float64(state[10]), rvec) :
                density_kg_m3(norm(rvec) - params.Re, Float64(state[10]))
        acceleration += ks_drag_acceleration_si(
            rvec,
            vvec,
            params,
            area_m2;
            density_kg_m3=density_value,
            drag_coefficient=drag_coefficient,
            mass_kg=mass_kg,
        )
    end
    # With h_KS = -ε, the unperturbed oscillator term is -(h_KS/2)u.
    ax, ay, az = acceleration
    u1, u2, u3, u4 = u
    perturbation_scale = 0.5 * r
    return @SVector [
        u_prime[1],
        u_prime[2],
        u_prime[3],
        u_prime[4],
        -0.5 * h_ks * u1 + perturbation_scale * (u1 * ax + u2 * ay + u3 * az),
        -0.5 * h_ks * u2 + perturbation_scale * (-u2 * ax + u1 * ay + u4 * az),
        -0.5 * h_ks * u3 + perturbation_scale * (-u3 * ax - u4 * ay + u1 * az),
        -0.5 * h_ks * u4 + perturbation_scale * (u4 * ax - u3 * ay + u2 * az),
        -r * dot(vvec, acceleration),
        r,
    ]
end

"""Evaluate the KS fictitious-time right-hand side in place."""
function nonlinear_ks_rhs!(derivative::AbstractVector,
    state::AbstractVector, params, area_m2::Real=0.0; kwargs...)
    length(derivative) == 10 || throw(ArgumentError(
        "KS derivative must have 10 components."))
    derivative .= nonlinear_ks_rhs(state, params, area_m2; kwargs...)
    return derivative
end

const ks_rhs! = nonlinear_ks_rhs!
const ks_rhs = nonlinear_ks_rhs

"""Advance the augmented KS state by one classical fourth-order RK step."""
function ks_rk4_step(state::AbstractVector, params, area_m2::Real,
    delta_s::Real; kwargs...)
    x = SVector{10,Float64}(state)
    step = Float64(delta_s)
    step > 0.0 || throw(ArgumentError("delta_s must be positive."))
    half_step = 0.5 * step
    k1 = nonlinear_ks_rhs(x, params, area_m2; kwargs...)
    k2 = nonlinear_ks_rhs(x + half_step * k1, params, area_m2; kwargs...)
    k3 = nonlinear_ks_rhs(x + half_step * k2, params, area_m2; kwargs...)
    k4 = nonlinear_ks_rhs(x + step * k3, params, area_m2; kwargs...)
    return x + (step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
end
