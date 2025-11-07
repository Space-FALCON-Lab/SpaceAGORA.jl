
# using PythonCall
using StaticArrays
using LinearAlgebra
# sys = pyimport("sys")

# Gravity models
@kwdef struct ConstantGravityModel <: AbstractForceTorqueModel
    μ::Float64 = 3.986004418e14  # Standard gravitational parameter for Earth, m^3/s^2
end

@kwdef struct InverseSquaredGravityModel <: AbstractForceTorqueModel
    μ::Float64 = 3.986004418e14  # Standard gravitational parameter for Earth, m^3/s^2
end

@kwdef struct InverseSquaredJ2GravityModel <: AbstractForceTorqueModel
    μ::Float64 = 3.986004418e14  # Standard gravitational parameter for Earth, m^3/s^2
    J2::Float64 = 1.08263e-3     # J2 coefficient for Earth
    Rp_m::Float64 = 6378137.0    # Equatorial radius of Earth in meters
end

# Calculate force/torque functions
# Model is the gravity model struct and x is the state vector from Complete_passage
function calcForceTorque(model::ConstantGravityModel, x::AbstractVector{Float64}, param::ODEParams)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    pos_ii = SVector{3, Float64}(x[1:3]) # Position in inertial frame, change to x.r if using StructArrays in Complete_passage
    mass = x[7]               # Mass of the spacecraft, change to x.m if using StructArrays in Complete_passage
    gravity_ii = -model.μ / norm(pos_ii)^2 * normalize(pos_ii)
    force_ii = mass * gravity_ii
    torque_ii = SVector{3, Float64}(zeros(3))
    return force_ii, torque_ii
end

function calcForceTorque(model::InverseSquaredGravityModel, x::AbstractVector{Float64}, param::ODEParams)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    pos_ii = SVector{3, Float64}(x[1:3]) # Position in inertial frame, change to x.r if using StructArrays in Complete_passage
    mass = x[7]               # Mass of the spacecraft, change to x.m if using StructArrays in Complete_passage
    gravity_ii = -model.μ / norm(pos_ii)^2 * normalize(pos_ii)
    force_ii = mass * gravity_ii
    torque_ii = SVector{3, Float64}(zeros(3))
    return force_ii, torque_ii
end

function calcForceTorque(model::InverseSquaredJ2GravityModel, x::AbstractVector{Float64}, param::ODEParams)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    pos_ii = SVector{3, Float64}(x[1:3]) # Position in inertial frame, change to x.r if using StructArrays in Complete_passage
    mass = x[7]               # Mass of the spacecraft, change to x.m if using StructArrays in Complete_passage
    r = norm(pos_ii)
    μ = model.μ
    J2 = model.J2
    Rp_m = model.Rp_m

    pos_ii_hat = normalize(pos_ii)
    r_squared = r^2
    gravity_ii_mag_spherical = -μ / r_squared

    x,y,z = pos_ii

    gravity_ii = gravity_ii_mag_spherical * pos_ii_hat + 3/2 * J2 * μ * Rp_m^2 / r^4 * [x/r*(5*z^2/r_squared - 1), y/r*(5*z^2/r_squared - 1), z/r*(5*z^2/r_squared - 3)] 

    force_ii = mass * gravity_ii
    torque_ii = SVector{3, Float64}(zeros(3))
    return force_ii, torque_ii
    
end