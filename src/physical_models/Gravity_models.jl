
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
function calcForceTorque(model::ConstantGravityModel, x::AbstractVector{Float64})::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    pos_ii = SVector{3, Float64}(x[1:3]) # Position in inertial frame, change to x.r if using StructArrays in Complete_passage
    mass = x[7]               # Mass of the spacecraft, change to x.m if using StructArrays in Complete_passage
    gravity_ii = -model.μ / norm(pos_ii)^2 * normalize(pos_ii)
    force_ii = mass * gravity_ii
    torque_ii = SVector{3, Float64}(zeros(3))
    return force_ii, torque_ii
end

function calcForceTorque(model::InverseSquaredGravityModel, x::AbstractVector{Float64})::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    pos_ii = SVector{3, Float64}(x[1:3]) # Position in inertial frame, change to x.r if using StructArrays in Complete_passage
    mass = x[7]               # Mass of the spacecraft, change to x.m if using StructArrays in Complete_passage
    gravity_ii = -model.μ / norm(pos_ii)^2 * normalize(pos_ii)
    force_ii = mass * gravity_ii
    torque_ii = SVector{3, Float64}(zeros(3))
    return force_ii, torque_ii
end

function calcForceTorque(model::InverseSquaredJ2GravityModel, x::AbstractVector{Float64})::Tuple{SVector{3, Float64}, SVector{3, Float64}}
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

# function gravity_const(pos_ii_mag::Float64, pos_ii::SVector{3, Float64}, p)
#     """

#     """

#     μ = p.μ
#     pos_ii_hat = normalize(pos_ii)

#     if config.cnf.drag_state == false
#         gravity_ii_mag = -μ / pos_ii_mag^2
#     else
#         gravity_ii_mag = -p.g_ref
#     end

#     return gravity_ii_mag * pos_ii_hat
# end

# function gravity_invsquared(pos_ii_mag::Float64, pos_ii::SVector{3, Float64}, p)
#     """

#     """
#     return -p.μ / pos_ii_mag^2 * normalize(pos_ii)
# end

# function gravity_invsquared_J2(r::Float64, pos_ii::SVector{3, Float64}, p)
#     """

#     """

#     μ = p.μ
#     J2 = p.J2

#     pos_ii_hat = normalize(pos_ii)
#     r_squared = r^2
#     gravity_ii_mag_spherical = -μ / r_squared

#     x,y,z = pos_ii

#     # gx = (-μ * x / r_cubed) * (1 + 3/2 * J2 * (p.Rp_e/r)^2 * (1 - 5*(z/r)^2))
#     # gy = (-μ * y / r_cubed) * (1 + 3/2 * J2 * (p.Rp_e/r)^2 * (1 - 5*(z/r)^2))
#     # gz = (-μ * z / r_cubed) * (1 + 3/2 * J2 * (p.Rp_e/r)^2 * (3 - 5*(z/r)^2))

#     return gravity_ii_mag_spherical * pos_ii_hat + 3/2 * J2 * μ * p.Rp_m^2 / r^4 * [x/r*(5*z^2/r_squared - 1), y/r*(5*z^2/r_squared - 1), z/r*(5*z^2/r_squared - 3)] 
# end

# function gravity_GRAM(pos_ii, lat, lon, alt, p, mass, vel_ii, el_time, atmosphere, args, gram)
#     """
#     Function to calculate the gravity of the planet using the GRAM model
#     :param pos_ii: position of the spacecraft in the inertial frame
#     :param lat: latitude of the spacecraft
#     :param lon: longitude of the spacecraft
#     :param alt: altitude of the spacecraft
#     :param timereal: current time of the simulation, used to calculate the elapsed time since the GRAM atmosphere start time
#     :param atmosphere: atmosphere of the planet
#     :return: gravity of the planet
#     """
#     if norm(pos_ii) - p.Rp_e > args[:EI] * 1e3
#         return gravity_invsquared_J2(norm(pos_ii), pos_ii, p)
#     end

#     position = gram.Position()
#     position.lat = rad2deg(lat)
#     position.lon = rad2deg(lon)
#     position.height = alt*1e-3

#     position.elapsedTime = el_time
#     atmosphere.setPosition(position)
    
#     atmosphere.update()
#     pos = atmosphere.getPosition()

#     gravity = pos.gravity
#     return -pyconvert(Any, gravity) * pos_ii / norm(pos_ii)
# end