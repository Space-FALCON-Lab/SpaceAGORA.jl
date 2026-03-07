
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
            force_ii += mass * gravity_n_bodies(et, pos_ii, planet, n_bodies_list[k])
        end
    end

    if _gravity_runtime_field(args, :gravity_harmonics, 0) == 1
        L = Int(_gravity_runtime_field(args, :L, 0))
        M = Int(_gravity_runtime_field(args, :M, 0))
        force_ii += mass * planet.L_PI' * acc_gravity_pines!(pos_pp, planet.Clm, planet.Slm, L, M, planet.μ, planet.Rp_e, planet)
    end

    return force_ii
end

# Calculate force/torque functions
# Model is the gravity model struct and x is the state vector from Complete_passage
function calcForceTorque!(model::ConstantGravityModel, x::AbstractVector{Float64}, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    m = param.m
    cnf = param.cnf

    pos_ii = SVector{3, Float64}(x[1:3]) # Position in inertial frame, change to x.r if using StructArrays in Complete_passage
    mass = x[7]               # Mass of the spacecraft, change to x.m if using StructArrays in Complete_passage
    gravity_ii = -m.planet.μ / norm(pos_ii)^2 * normalize(pos_ii)
    force_ii = mass * gravity_ii
    if model.gravity_gradient && param.orientation_sim
        tau_gg = gravity_gradient(param.inertia_matrix, pos_ii, model.planet.μ)
        return force_ii, tau_gg
    else
        return force_ii, SVector{3, Float64}(0.0, 0.0, 0.0)
    end
end

function calcForceTorque(model::InverseSquaredGravityModel, x::ComponentVector, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    # m = param.m
    # cnf = param.cnf
    # planet = param.args.environment_model.planet

    pos_ii = SVector{3, Float64}(x.pos) # Position in inertial frame, change to x.r if using StructArrays in Complete_passage
    mass = x.mass               # Mass of the spacecraft, change to x.m if using StructArrays in Complete_passage
    gravity_ii = -param.args.environment_model.planet.μ / norm(pos_ii)^2 * normalize(pos_ii)

    # cnf.gravity_cent_ii = mass * gravity_ii # Store gravity in config for other uses

    force_ii = mass * gravity_ii
    torque_ii = SVector{3, Float64}(zeros(3))
    return force_ii, torque_ii
end

function calcForceTorque(model::InverseSquaredJ2GravityModel, x::ComponentVector, param::ODEParams, i::Int64)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
        # m = param.m
        # cnf = param.cnf

    pos_ii = SVector{3, Float64}(x.pos) # Position in inertial frame, change to x.r if using StructArrays in Complete_passage
    mass = x.mass               # Mass of the spacecraft, change to x.m if using StructArrays in Complete_passage
    r = norm(pos_ii)
    μ = param.args.environment_model.planet.μ
    J2 = param.args.environment_model.planet.J2
    Rp_m = param.args.environment_model.planet.Rp_m

    pos_ii_hat = normalize(pos_ii)
    r_squared = r^2
    gravity_ii_mag_spherical = -μ / r_squared

    x,y,z = pos_ii

    gravity_ii = gravity_ii_mag_spherical * pos_ii_hat + 3/2 * J2 * μ * Rp_m^2 / r^4 * [x/r*(5*z^2/r_squared - 1), y/r*(5*z^2/r_squared - 1), z/r*(5*z^2/r_squared - 3)]

    # cnf.gravity_cent_ii = mass * gravity_ii # Store gravity in config for other uses

    force_ii = mass * gravity_ii
    torque_ii = SVector{3, Float64}(zeros(3))
    return force_ii, torque_ii
    
end

"""
    gravity_gradient(J::SMatrix{3,3,Float64}, rVec::SVector{3,Float64}, μ::Float64)

Calculates the gravity gradient torque exerted on a spacecraft due to the size of the spacecraft

# Args
- `J`: The inertia matrix of the spacecraft in the inertial frame [kg·m²].
- `rVec`: The position vector of the spacecraft in the inertial frame [meters].
- `μ`: The gravitational parameter of the central body [m³/s²].

# Returns
- A 3-element `SVector` representing the gravity gradient torque in the body frame `[τ_x, τ_y, τ_z]` [N·m].
"""
function gravity_gradient(J::SMatrix{3,3,Float64}, rVec::SVector{3,Float64}, μ::Float64)
    r = norm(rVec)
    r_hat = rVec / r
    return 3*μ/r^3 * cross(r_hat, J * r_hat)
end
