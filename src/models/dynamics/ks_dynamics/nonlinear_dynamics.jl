## Propagation parameters and energy convention

"""
Planet constants of the KS propagator, in SI units: equatorial radius `Re` (m),
gravitational parameter `μ` (m³/s²), the `J2` zonal coefficient and the rotation
rate `Ω` (rad/s). The KS frame is non-rotating with its z axis along the planet's
pole: J2 acts about that axis and the atmosphere co-rotates about it at `Ω`.
"""
Base.@kwdef struct KSPropagationParams
    Re::Float64
    μ::Float64
    J2::Float64 = 0.0
    Ω::Float64 = 0.0
end

"""Return the KS energy parameter `h_KS = -ε` of a specific orbital energy `ε` (J/kg)."""
ks_energy_parameter(specific_energy::Real) = -Float64(specific_energy)

"""Return the specific orbital energy `ε = -h_KS` (J/kg) of a KS energy parameter."""
specific_energy_from_ks(h_ks::Real) = -Float64(h_ks)

"""Return the cross-product matrix of the planet rotation vector `[0, 0, Ω]`."""
function ks_rotation_cross_matrix(params)
    omega = Float64(params.Ω)
    return @SMatrix [0.0 -omega 0.0; omega 0.0 0.0; 0.0 0.0 0.0]
end

"""Return the velocity relative to the co-rotating atmosphere, `v - Ω × r`."""
@inline function ks_relative_velocity(rvec, vvec, params)
    omega = Float64(params.Ω)
    return SVector{3,Float64}(vvec[1] + omega * rvec[2], vvec[2] - omega * rvec[1], vvec[3])
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

"""Return the drag acceleration `-½ ρ C_D A/m |V| V` (SI), with `V = v - Ω × r`."""
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
    relative_velocity = ks_relative_velocity(rvec, vvec, params)
    speed = norm(relative_velocity)
    speed <= eps(Float64) && return @SVector [0.0, 0.0, 0.0]
    return -0.5 * density * Float64(drag_coefficient) * Float64(area_m2) /
        mass * speed * relative_velocity
end

"""
    ks_density_value(source, position, elapsed_time_s, params) -> Float64

Nonnegative density (kg/m³) of a density source: a number, or a callable
`(altitude_m, elapsed_s, position)` or `(altitude_m, elapsed_s)`, with the
altitude `|r| - Re` (m) above the equatorial radius.
"""
@inline function ks_density_value(source, position, elapsed_time_s, params)
    source isa Real && return max(0.0, Float64(source))
    altitude_m = norm(position) - params.Re
    value = applicable(source, altitude_m, elapsed_time_s, position) ?
        source(altitude_m, elapsed_time_s, position) :
        source(altitude_m, elapsed_time_s)
    return max(0.0, Float64(value))
end

## State conversion

"""
    cartesian_to_ks_state(position_ii_m, velocity_ii_m, params; elapsed_time_s=0)

Create the 10-component KS state `[u, u′, h_KS, t]` from a Cartesian state in
the KS frame. Distances, velocity, gravitational parameter, and time use SI units.
"""
function cartesian_to_ks_state(position_ii_m, velocity_ii_m, params; elapsed_time_s::Real=0.0)
    r = SVector{3, Float64}(position_ii_m)
    v = SVector{3, Float64}(velocity_ii_m)
    u = cartesian_position_to_ks_coordinate(r)
    u_prime = cartesian_velocity_to_ks_derivative(v, u)
    h_ks = ks_energy_parameter(0.5 * dot(v, v) - params.μ / norm(r))
    return MVector{10,Float64}(u..., u_prime..., h_ks, Float64(elapsed_time_s))
end

"""
Convert a 10-component KS state to its Cartesian position and velocity, KS
energy parameter, specific orbital energy (J/kg) and elapsed time.
"""
function ks_state_to_cartesian(state)
    u = SVector{4, Float64}(state[1], state[2], state[3], state[4])
    u_prime = SVector{4, Float64}(state[5], state[6], state[7], state[8])
    return (
        position_ii_m=ks_position(u),
        velocity_ii_m=ks_velocity(u, u_prime),
        h_ks=Float64(state[9]),
        specific_energy_j_kg=specific_energy_from_ks(state[9]),
        elapsed_time_s=Float64(state[10]),
    )
end

## Nonlinear KS dynamics

# The KS right-hand side for a perturbing acceleration `a` (SI):
#   u″ = -(h_KS/2) u + (|r|/2) Λ(u)ᵀ a,  h_KS′ = -|r| v·a,  t′ = |r|.
# With h_KS = -ε, -(h_KS/2) u is the unperturbed harmonic-oscillator term, and
# the energy parameter changes by the work of the perturbation.
@inline function ks_rhs_from_acceleration(u, u_prime, h_ks, radius, velocity, acceleration)
    u_second = -0.5 * h_ks * u + 0.5 * radius * (transpose(ks_lambda_matrix(u)) * acceleration)
    return SVector{10,Float64}(u_prime..., u_second...,
        -radius * dot(velocity, acceleration), radius)
end

"""
    nonlinear_ks_rhs(state, params, area_m2=0; density_kg_m3=0, drag_coefficient=0,
        mass_kg=1, use_drag=false) -> SVector{10}

Derivative of the KS state `[u, u′, h_KS, t]` with respect to fictitious time,
under J2 and, with `use_drag`, drag on `area_m2`. `density_kg_m3` is a number or
a density callable (see `ks_density_value`).
"""
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
    rvec = ks_position(u)
    r = dot(u, u)
    vvec = ks_velocity_with_radius(u, u_prime, r)
    acceleration = ks_j2_acceleration_with_radius(rvec, r, params)
    if use_drag
        acceleration += ks_drag_acceleration_si(
            rvec,
            vvec,
            params,
            area_m2;
            density_kg_m3=ks_density_value(density_kg_m3, rvec, Float64(state[10]), params),
            drag_coefficient=drag_coefficient,
            mass_kg=mass_kg,
        )
    end
    return ks_rhs_from_acceleration(u, u_prime, Float64(state[9]), r, vvec, acceleration)
end

"""Evaluate [`nonlinear_ks_rhs`](@ref) into `derivative`."""
function nonlinear_ks_rhs!(derivative::AbstractVector,
    state::AbstractVector, params, area_m2::Real=0.0; kwargs...)
    length(derivative) == 10 || throw(ArgumentError(
        "KS derivative must have 10 components."))
    derivative .= nonlinear_ks_rhs(state, params, area_m2; kwargs...)
    return derivative
end

"""
    ks_rk4_step(state, params, area_m2, delta_s; kwargs...) -> SVector{10}

Advance the KS state by one classical fourth-order Runge-Kutta step of
fictitious time `delta_s`; `kwargs` are those of [`nonlinear_ks_rhs`](@ref).
"""
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
