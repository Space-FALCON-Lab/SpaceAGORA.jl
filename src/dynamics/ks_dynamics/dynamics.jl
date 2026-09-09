"""Planet constants required by the reusable KS orbital propagator."""
Base.@kwdef struct KSPropagationParams
    Re::Float64
    μ::Float64
    J2::Float64 = 0.0
    Ω::Float64 = 0.0
end

"""Return the paper-convention KS energy parameter `h_KS = -specific_energy`."""
ks_energy_parameter(specific_energy::Real) = -Float64(specific_energy)

"""Recover specific orbital energy from the paper-convention KS parameter."""
specific_energy_from_ks(h_ks::Real) = -Float64(h_ks)

function _ks_skew_rotation(params)
    omega = Float64(params.Ω)
    return @SMatrix [0.0 -omega 0.0; omega 0.0 0.0; 0.0 0.0 0.0]
end

"""Return the J2 perturbing acceleration in SI units."""
function ks_j2_acceleration_si(rvec, params)
    r = norm(rvec)
    r <= eps(Float64) && return @SVector [0.0, 0.0, 0.0]
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
    relative_velocity = vvec - _ks_skew_rotation(params) * rvec
    speed = norm(relative_velocity)
    speed <= eps(Float64) && return @SVector [0.0, 0.0, 0.0]
    return -0.5 * density * Float64(drag_coefficient) * Float64(area_m2) /
        mass * speed * relative_velocity
end

"""
    cartesian_to_ks_state(position_ii_m, velocity_ii_m, params; elapsed_time_s=0)

Create the 10-component KS state `[u, u′, h, t]` from an inertial Cartesian
state. Distances, velocity, gravitational parameter, and time use SI units.
"""
function cartesian_to_ks_state(position_ii_m, velocity_ii_m, params; elapsed_time_s::Real=0.0)
    u = SVector{4, Float64}(_ks_coordinate_from_position(position_ii_m)...)
    u_prime = SVector{4, Float64}(_ks_derivative_from_velocity(velocity_ii_m, u)...)
    r = _ks_lambda(u) * u
    v = ks_velocity(u, u_prime)
    energy = 0.5 * dot(v, v) - params.μ / norm(r)
    h_ks = ks_energy_parameter(energy)
    return collect(vcat(u, u_prime, h_ks, Float64(elapsed_time_s)))
end

"""Convert a 10-component KS state back to inertial Cartesian state."""
function ks_state_to_cartesian(state)
    u = SVector{4, Float64}(state[1:4])
    u_prime = SVector{4, Float64}(state[5:8])
    return (
        position_ii_m=_ks_lambda(u) * u,
        velocity_ii_m=ks_velocity(u, u_prime),
        h_ks=Float64(state[9]),
        energy_parameter=Float64(state[9]),
        specific_energy_j_kg=specific_energy_from_ks(state[9]),
        elapsed_time_s=Float64(state[10]),
    )
end

"""Evaluate the KS fictitious-time right-hand side."""
function ks_rhs(
    state::AbstractVector,
    params,
    area_m2::Real=0.0;
    config=nothing,
    density_kg_m3=0.0,
    drag_coefficient::Real=0.0,
    mass_kg::Real=1.0,
    use_drag::Bool=false,
)
    u = SVector{4, Float64}(state[1:4])
    u_prime = SVector{4, Float64}(state[5:8])
    h_ks = Float64(state[9])
    rvec = _ks_lambda(u) * u
    r = dot(u, u)
    vvec = ks_velocity(u, u_prime)
    acceleration = ks_j2_acceleration_si(rvec, params)
    if use_drag
        if config !== nothing
            drag_coefficient = config.drag_coefficient
            mass_kg = config.mass_kg
        end
        density_value = density_kg_m3 isa Function ?
            density_kg_m3(norm(rvec) - params.Re, Float64(state[10])) :
            density_kg_m3
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
    acceleration4 = SVector(acceleration[1], acceleration[2], acceleration[3], 0.0)
    du = u_prime
    # Paper convention: h_KS = -ε and ω_KS² = h_KS/2. The older
    # implementation stored -2ε, which required the equivalent -h*u/4 term.
    du_prime = -0.5 * h_ks .* u + 0.5 * r .* (transpose(_ks_L(u)) * acceleration4)
    dh_ks = -r * dot(vvec, acceleration)
    dt = r
    return collect(vcat(du, du_prime, dh_ks, dt))
end

"""Advance a KS state by one fixed classical-RK4 fictitious-time step."""
function ks_rk4_step(
    state::AbstractVector,
    params,
    area_m2::Real,
    delta_s::Real;
    config=nothing,
    density_kg_m3=0.0,
    drag_coefficient::Real=0.0,
    mass_kg::Real=1.0,
    use_drag::Bool=false,
)
    x = Float64.(collect(state))
    step = Float64(delta_s)
    f(z) = ks_rhs(
        z,
        params,
        area_m2;
        config=config,
        density_kg_m3=density_kg_m3,
        drag_coefficient=drag_coefficient,
        mass_kg=mass_kg,
        use_drag=use_drag,
    )
    k1 = f(x)
    k2 = f(x .+ 0.5 * step .* k1)
    k3 = f(x .+ 0.5 * step .* k2)
    k4 = f(x .+ step .* k3)
    return x .+ (step / 6.0) .* (k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4)
end
