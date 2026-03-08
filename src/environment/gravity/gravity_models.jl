
# using PythonCall
using StaticArrays
using LinearAlgebra
using ComponentArrays
# sys = pyimport("sys")

# Gravity models
@kwdef struct ConstantGravityModel <: AbstractForceTorqueModel
    # μ::Float64 = 3.986004418e14  # Standard gravitational parameter for Earth, m^3/s^2
    gravity_gradient::Bool = false
end

@kwdef struct InverseSquaredGravityModel <: AbstractForceTorqueModel
    # μ::Float64 = 3.986004418e14  # Standard gravitational parameter for Earth, m^3/s^2
    gravity_gradient::Bool = false
end

@kwdef struct InverseSquaredJ2GravityModel <: AbstractForceTorqueModel
    # μ::Float64 = 3.986004418e14  # Standard gravitational parameter for Earth, m^3/s^2
    # J2::Float64 = 1.08263e-3     # J2 coefficient for Earth
    # Rp_m::Float64 = 6378137.0    # Equatorial radius of Earth in meters
    gravity_gradient::Bool = false
end

@inline function _gravity_runtime_field(args, name::Symbol, default)
    if args !== nothing && hasproperty(args, name)
        return getproperty(args, name)
    end
    if args !== nothing && applicable(get, args, name, default)
        return get(args, name, default)
    end
    if args !== nothing && applicable(getindex, args, name)
        return getindex(args, name)
    end
    return default
end

@inline function _inverse_squared_gravity_accel(pos_ii::SVector{3, Float64}, planet)::SVector{3, Float64}
    r = norm(pos_ii)
    return -planet.μ / r^2 * normalize(pos_ii)
end

@inline function _inverse_squared_j2_gravity_accel(pos_ii::SVector{3, Float64}, planet)::SVector{3, Float64}
    r = norm(pos_ii)
    μ = planet.μ
    J2 = planet.J2
    Rp_m = planet.Rp_m
    x, y, z = pos_ii
    r_squared = r^2
    gravity_ii_mag_spherical = -μ / r_squared
    pos_ii_hat = normalize(pos_ii)
    j2_term = SVector{3, Float64}(
        x / r * (5 * z^2 / r_squared - 1),
        y / r * (5 * z^2 / r_squared - 1),
        z / r * (5 * z^2 / r_squared - 3),
    )
    return gravity_ii_mag_spherical * pos_ii_hat + 3 / 2 * J2 * μ * Rp_m^2 / r^4 * j2_term
end

@inline function _gravity_gradient_torque_body(
    model,
    pos_ii::SVector{3, Float64},
    x,
    param::ODEParams,
    i::Int64
)::SVector{3, Float64}
    if !model.gravity_gradient || !param.args.mission_configuration.orientation_sim
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if i < 1 || i > length(param.args.dynamics_model.spacecraft)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    if !hasproperty(x, :q)
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    r = norm(pos_ii)
    if !isfinite(r) || r <= 0.0
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    inertia_tensor = param.args.dynamics_model.spacecraft[i].inertia_tensor
    q_ib = getproperty(x, :q)
    q_body = SVector{4, Float64}(Float64(q_ib[1]), Float64(q_ib[2]), Float64(q_ib[3]), Float64(q_ib[4]))
    r_body = rot(q_body) * pos_ii
    return gravity_gradient(inertia_tensor, r_body, param.args.environment_model.planet.μ)
end

@inline function _gravity_gradient_torque_body(
    model,
    x::StateSample,
    planet
)::SVector{3, Float64}
    if !model.gravity_gradient || x.q_ib === nothing || x.spacecraft === nothing
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    r = norm(x.pos_ii)
    if !isfinite(r) || r <= 0.0
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    inertia_tensor = x.spacecraft.inertia_tensor
    r_body = rot(x.q_ib) * x.pos_ii
    return gravity_gradient(inertia_tensor, r_body, planet.μ)
end

"""
    j2_secular_rates(a, e, i, planet) -> (Ωdot, ωdot)

First-order secular drift rates for the J2-only zonal gravity model, using the
standard expressions given in Vallado and Montenbruck & Gill.

Returns the right ascension of ascending node drift `Ωdot` and argument of
periapsis drift `ωdot` in rad/s.
"""
@inline function j2_secular_rates(
    a::Float64,
    e::Float64,
    i::Float64,
    planet
)::Tuple{Float64, Float64}
    if !isfinite(planet.J2) || planet.J2 == 0.0 || !isfinite(a) || a <= 0.0 || !isfinite(e) || e < 0.0 || e >= 1.0
        return 0.0, 0.0
    end
    p = a * (1.0 - e^2)
    if !isfinite(p) || p <= 0.0
        return 0.0, 0.0
    end
    n = sqrt(planet.μ / a^3)
    if !isfinite(n) || n <= 0.0
        return 0.0, 0.0
    end
    scale = planet.J2 * (planet.Rp_e / p)^2
    Ωdot = -1.5 * n * scale * cos(i)
    ωdot = 0.75 * n * scale * (5.0 * cos(i)^2 - 1.0)
    return Ωdot, ωdot
end

function aerobraking_gravity_force_ii(
    gm_code::Integer,
    mass::Float64,
    pos_ii::SVector{3, Float64},
    vel_ii::SVector{3, Float64},
    pos_pp::SVector{3, Float64},
    lat::Float64,
    lon::Float64,
    alt::Float64,
    planet,
    et::Float64,
    args,
    gram_atmosphere,
    gram,
    n_bodies_list,
)::SVector{3, Float64}
    gravity_ii = if gm_code == 2
        _inverse_squared_j2_gravity_accel(pos_ii, planet)
    else
        # The legacy gm=0/1/3 aerobraking paths all reduced to inverse-squared gravity.
        _inverse_squared_gravity_accel(pos_ii, planet)
    end

    force_ii = mass * gravity_ii

    n_bodies = _gravity_runtime_field(args, :n_bodies, ())
    if !isempty(n_bodies)
        for k in eachindex(n_bodies)
            force_ii += mass * PerturbationEffectors.gravity_n_bodies(et, pos_ii, planet, n_bodies_list[k])
        end
    end

    if _gravity_runtime_field(args, :gravity_harmonics, 0) == 1
        L = Int(_gravity_runtime_field(args, :L, 0))
        M = Int(_gravity_runtime_field(args, :M, 0))
        force_ii += mass * planet.L_PI' * PerturbationEffectors.acc_gravity_pines!(pos_pp, planet.Clm, planet.Slm, L, M, planet.μ, planet.Rp_e, planet)
    end

    return force_ii
end

function calcForceTorque(model::ConstantGravityModel, x::ComponentVector, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    pos_ii = SVector{3, Float64}(x.pos)
    mass = x.mass
    gravity_ii = _inverse_squared_gravity_accel(pos_ii, param.args.environment_model.planet)
    force_ii = mass * gravity_ii
    torque_body = _gravity_gradient_torque_body(model, pos_ii, x, param, i)
    return force_ii, torque_body
end

@inline function wrench(
    model::ConstantGravityModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    gravity_ii = _inverse_squared_gravity_accel(x.pos_ii, env.planet)
    force_ii = x.mass_kg * gravity_ii
    torque_body = _gravity_gradient_torque_body(model, x, env.planet)
    return force_ii, torque_body
end

function calcForceTorque(model::InverseSquaredGravityModel, x::ComponentVector, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    pos_ii = SVector{3, Float64}(x.pos)
    mass = x.mass
    gravity_ii = _inverse_squared_gravity_accel(pos_ii, param.args.environment_model.planet)
    force_ii = mass * gravity_ii
    torque_body = _gravity_gradient_torque_body(model, pos_ii, x, param, i)
    return force_ii, torque_body
end

@inline function wrench(
    model::InverseSquaredGravityModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    gravity_ii = _inverse_squared_gravity_accel(x.pos_ii, env.planet)
    force_ii = x.mass_kg * gravity_ii
    torque_body = _gravity_gradient_torque_body(model, x, env.planet)
    return force_ii, torque_body
end

function calcForceTorque(model::InverseSquaredJ2GravityModel, x::ComponentVector, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    pos_ii = SVector{3, Float64}(x.pos)
    mass = x.mass
    gravity_ii = _inverse_squared_j2_gravity_accel(pos_ii, param.args.environment_model.planet)
    force_ii = mass * gravity_ii
    torque_body = _gravity_gradient_torque_body(model, pos_ii, x, param, i)
    return force_ii, torque_body
end

@inline environment_requirements(::InverseSquaredJ2GravityModel) = EffectorEnvironmentRequirements(planet_frame=true)

@inline function wrench(
    model::InverseSquaredJ2GravityModel,
    x::StateSample,
    env::EnvironmentSample,
    t::Float64,
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    planet_frame = env.planet_frame
    planet_frame === nothing && throw(ArgumentError("InverseSquaredJ2GravityModel wrench requires env.planet_frame."))
    gravity_pp = _inverse_squared_j2_gravity_accel(planet_frame.pos_pp, env.planet)
    gravity_ii = planet_frame.l_pi' * gravity_pp
    force_ii = x.mass_kg * gravity_ii
    torque_body = _gravity_gradient_torque_body(model, x, env.planet)
    return force_ii, torque_body
end

"""
    gravity_gradient(J::SMatrix{3,3,Float64}, rVec::SVector{3,Float64}, μ::Float64)

Calculate the body-frame gravity-gradient torque for a rigid spacecraft in a
central field.

# Args
- `J`: The spacecraft inertia tensor expressed in the body frame [kg·m²].
- `rVec`: The spacecraft position vector expressed in the same body frame [m].
- `μ`: The gravitational parameter of the central body [m³/s²].

# Returns
- A 3-element `SVector` representing the gravity-gradient torque in the body
  frame `[τ_x, τ_y, τ_z]` [N·m].
"""
function gravity_gradient(J::SMatrix{3,3,Float64}, rVec::SVector{3,Float64}, μ::Float64)
    r = norm(rVec)
    if !isfinite(r) || r <= 0.0
        return SVector{3, Float64}(0.0, 0.0, 0.0)
    end
    r_hat = rVec / r
    return 3*μ/r^3 * cross(r_hat, J * r_hat)
end
