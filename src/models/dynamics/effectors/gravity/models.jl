
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
