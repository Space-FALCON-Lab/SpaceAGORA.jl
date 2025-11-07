using SPICE
using LoopVectorization
using AssociatedLegendrePolynomials
using LinearAlgebra
using SatelliteToolbox
using SatelliteToolboxGeomagneticField
using DateFormats
using CSV
using DataFrames

include("../utils/quaternion_utils.jl")
include("Planet_data.jl")
# import .config
# Define delta function
δ(x,y) = ==(x,y)
δ(x) = δ(x,0)

# Constants for the Tilted Dipole Model (Epoch 2020.0)
# Reference magnetic field strength at the equator on the Earth's surface.
const B0_2020 = 3.12e-5  # Tesla
# WGS84 Earth radius for the model.
const R_EARTH_MODEL = 6371.2e3 # meters
# North Magnetic Pole location (geocentric).
const POLE_LAT_2020 = deg2rad(80.7)  # Radians
const POLE_LON_2020 = deg2rad(-72.7) # Radians

# Pre-calculate the magnetic pole axis vector in ECEF for efficiency.
const M_HAT_ECEF = SVector{3, Float64}(
    cos(POLE_LAT_2020) * cos(POLE_LON_2020),
    cos(POLE_LAT_2020) * sin(POLE_LON_2020),
    sin(POLE_LAT_2020)
)

# N-body gravity perturbation model
@kwdef struct NBodyGravityModel <: AbstractForceTorqueModel
    body_names::Vector{String} = String[]  # Names of celestial bodies to include
    primary_body_name::String = "Earth" # Name of the primary body
    planet::Planet = Planet() # Planet data for primary body
end

struct GravitationalHarmonicsModel <: AbstractForceTorqueModel
    L::Int64 # Degree
    M::Int64 # Order
    C::Matrix{Float64} # Cosine coefficients
    S::Matrix{Float64} # Sine coefficients
    A::Matrix{Float64} # Preallocated ALF array
    R::Vector{Float64} # Preallocated real terms powers
    I::Vector{Float64} # Preallocated imaginary terms vector
    VR01::Matrix{Float64} # Preallocated VR01 array
    VR11::Matrix{Float64} # Preallocated VR11 array
    N1::Matrix{Float64} # Preallocated N1 array
    N2::Matrix{Float64} # Preallocated N2 array
    planet::Planet # Planet data for primary body
end

struct SolarRadiationPressureModel <: AbstractForceTorqueModel
    Cr::Float64 # Reflectivity coefficient
    A::Float64  # Cross-sectional area in m^2
    AU_m::Float64 # Astronomical unit in meters
end
# Constructor to get planet data
function NBodyGravityModel(body_names::Vector{String}, primary_body_name::String="Earth")
    planet = planet_data(primary_body_name)
    return NBodyGravityModel(body_names=body_names, primary_body_name=primary_body_name, planet=planet)
end

function GravitationalHarmonicsModel(L::Int64, M::Int64, coefficients_file::String, planet_name::String)
    harmonics_data = CSV.read(coefficients_file, DataFrame)
    total_data_size = size(harmonics_data, 1)
    degree = maximum(harmonics_data[:, 1]) + 1
    C = zeros(degree, degree)
    S = zeros(degree, degree)
    for i=1:total_data_size
        l = harmonics_data[i, 1] + 1 # Get the degree, l, from the data and convert to an index (subtract 1 because the data starts at 2nd degree coefficient)
        m = harmonics_data[i, 2] + 1 # Get the order, m, from the data and convert to an index (add 1 because the data starts at 0th order coefficient)
        C[l, m] = harmonics_data[i, 3]
        S[l, m] = harmonics_data[i, 4]
    end

    N1 = MMatrix{L+4, L+4, Float64}(zeros(L+4, L+4))
    N2 = MMatrix{L+4, L+4, Float64}(zeros(L+4, L+4))
    VR01 = MMatrix{L+1, L+1, Float64}(zeros(L+1, L+1))
    VR11 = MMatrix{L+1, L+1, Float64}(zeros(L+1, L+1))
    sqrt_2 = sqrt(2.0)
    for m = 0:M+2
        j = m + 1
        for l = m+2:L+2
            i = l + 1
            N1[i, j] = √((2*l+1)*(2*l-1)/(l+m)/(l-m))
            N2[i, j] = √((l+m-1)*(2*l+1)*(l-m-1)/(2*l-3)/(l+m)/(l-m))
        end
    end
    N1 = SMatrix{L+4, L+4, Float64}(N1)
    N2 = SMatrix{L+4, L+4, Float64}(N2)

    for l = 0:L
        i = l + 1
        for m = 0:min(M, l)
            j = m + 1
            divisor = m == 0 ? sqrt_2 : 1
            VR01[i, j] = sqrt((l-m)*(l+m+1)) / divisor
            VR11[i, j] = sqrt((2*l+1)*(l+m+2)*(l+m+1)/(2*l+3)) / divisor
        end
    end
    VR01 = SMatrix{L+1, L+1, Float64}(VR01)
    VR11 = SMatrix{L+1, L+1, Float64}(VR11)

    A = MMatrix{L+4, L+4, Float64}(zeros(L+4, L+4))
    R = MVector{L+4, Float64}(zeros(L+4))
    I = MVector{L+4, Float64}(zeros(L+4))
    A[1, 1] = 1
    # Fill the diagonal elements of A
    for l = 1:L+2
        i = l + 1
        A[i, i] = sqrt((2*l+1)/(2*l))*A[i-1, i-1]
    end

    planet = planet_data(planet_name)

    return GravitationalHarmonicsModel(
        L,
        M,
        C[1:L, 1:L],
        S[1:L, 1:L],
        A,
        R,
        I,
        VR01,
        VR11,
        N1,
        N2,
        planet
    )
end

"""
    calcForceTorque(model::NBodyGravityModel, x::AbstractVector{Float64}, ODEParams)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
"""
function calcForceTorque(model::NBodyGravityModel, x::AbstractVector{Float64}, param::ODEParams)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    pos_ii = SVector{3, Float64}(x[1:3]) # Position in inertial frame, change to x.r if using StructArrays in Complete_passage
    mass = x[7]               # Mass of the spacecraft, change to x.m if using StructArrays in Complete_passage
    primary_body_name = model.primary_body_name
    force_ii = MVector{3, Float64}(0.0, 0.0, 0.0) # Initialize force vector
    for body_name in model.body_names
        barycenter_bodies = ["mars", "jupiter", "saturn", "uranus", "neptune", "earth"]
        if !isnothing(findfirst(==(primary_body_name), barycenter_bodies))
            primary_body_name *= "_barycenter"
        end

        if !isnothing(findfirst(==(body_name), barycenter_bodies))
            body_name *= "_barycenter"
        end

        pos_primary_k = model.planet.J2000_to_pci*SVector{3, Float64}(spkpos(body_name, param.cnf.et, "J2000", "none", primary_body_name)[1]) * 1e3
        pos_spacecraft_k = pos_primary_k - pos_ii
        pos_spacecraft_k_mag = norm(pos_spacecraft_k)

        force_ii += mass * model.planet.μ * ((pos_spacecraft_k / pos_spacecraft_k_mag^3) - (pos_primary_k / norm(pos_primary_k)^3))
    end
    return force_ii, SVector{3, Float64}(0.0, 0.0, 0.0)
end

function calcForceTorque(model::GravitationalHarmonicsModel, x::AbstractVector{Float64}, param::ODEParams)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    rVec_cart = SVector{3, Float64}(x[1:3])
    mass = x[7]               # Mass of the spacecraft, change to x.m if using StructArrays in Complete_passage

    RE = model.planet.Rp_e
    r = norm(rVec_cart)
    s,t,u=normalize(rVec_cart)
    L = model.L
    M = model.M
    sqrt_2 = sqrt(2.0)
    @fastmath begin
        model.A[2, 1] = u*sqrt(3.0)
        # Fill the off diagonal elements of A
        @inbounds @simd for n = 1:L+1
            i = n + 1
            model.A[i+1, i] = u*sqrt(2.0*n+3.0)*model.A[i, i]
        end
        # Fill the rest of A
        @inbounds for m = 0:M+1
            j = m + 1
            @inbounds for l = m+2:L+1
                i = l + 1
                model.A[i, j] = u*model.N1[i, j]*model.A[i-1, j] - model.N2[i, j]*model.A[i-2, j]
            end
            R_term = model.R[j-1]
            I_term = model.I[j-1]
            model.R[j] = ifelse(m == 0, 1, s*R_term - t*I_term)
            model.I[j] = ifelse(m == 0, 0, s*I_term + t*R_term)
        end

        ρ = RE/r
        ρ_np1 = -model.planet.μ/r * ρ
        a1 = a2 = a3 = a4 = 0.0
        @inbounds for l = 1:L
            i = l + 1
            ρ_np1 *= ρ
            sum1 = 0.0
            sum2 = 0.0
            sum3 = 0.0
            sum4 = 0.0
            @inbounds for m = 0:min(l, M)
                j = m + 1
                C = model.C[i, j]
                S = model.S[i, j]
                R_term = model.R[j-1]
                I_term = model.I[j-1]
                D =                   (C*model.R[j] + S*model.I[j])     * sqrt_2
                E = ifelse(m == 0, 0.0, (C*R_term + S*I_term) * sqrt_2)
                F = ifelse(m == 0, 0.0, (S*R_term - C*I_term) * sqrt_2)

                Avv00, Avv01, Avv11 = model.A[i, j], model.VR01[i, j]*model.A[i, j+1], model.VR11[i, j]*model.A[i+1, j+1]

                sum1 += m * Avv00 * E
                sum2 += m * Avv00 * F
                sum3 +=     Avv01 * D
                sum4 +=     Avv11 * D
            end
            rr = ρ_np1/RE
            a1 += rr * sum1
            a2 += rr * sum2
            a3 += rr * sum3
            a4 -= rr * sum4
        end
    end
    return mass * SVector{3, Float64}(-a1 - s*a4, -a2 - t*a4, -a3 - u*a4), SVector{3, Float64}(0.0, 0.0, 0.0)
end
"""
    get_magnetic_field_dipole(r_ecef::AbstractVector)

Calculates the Earth's magnetic field using a fast tilted dipole approximation.

This model is significantly faster than the full WMM and is suitable for use
inside performance-critical code like numerical integrators.

# Args

- `r_ecef`: The position vector of the spacecraft in an Earth-Centered,
            Earth-Fixed (ECEF) frame [meters].

# Returns

- A 3-element `SVector` representing the magnetic field `[Bx, By, Bz]` in the
  ECEF frame [Tesla].
"""
function get_magnetic_field_dipole(r_ecef::AbstractVector, L_PI::MMatrix{3, 3, Float64})
    r_norm = norm(r_ecef)
    r_hat = r_ecef / r_norm

    # Cosine of the magnetic colatitude
    cos_colat = dot(r_hat, M_HAT_ECEF)

    # Dipole field equation. This is a standard formulation.
    B_ecef = -B0_2020 * (R_EARTH_MODEL / r_norm)^3 * (M_HAT_ECEF - 3 * cos_colat * r_hat)

    return L_PI' * B_ecef
end

"""
    get_magnetic_field(date::DateTime, lat_deg::Number, lon_deg::Number, alt_m::Number)

Computes the Earth's magnetic field vector in the local North-East-Down (NED)
frame using the World Magnetic Model (WMM).

The function automatically uses the correct WMM version based on the input `date`.

# Args

- `date`: The `DateTime` of the measurement.
- `lat_deg`: The geodetic latitude of the observer [degrees].
- `lon_deg`: The longitude of the observer [degrees].
- `alt_m`: The altitude above the WGS84 ellipsoid [meters].

# Returns

- A 3-element `SVector` representing the magnetic field `[B_north, B_east, B_down]`
  in nanoTeslas [nT].
"""
function get_magnetic_field(date::DateTime, lat_rad::Number, lon_rad::Number, alt_m::Number, L_PI::MMatrix{3, 3, Float64})
    # println("Calculating magnetic field at lat: $lat_rad, lon: $lon_rad, alt: $alt_m, date: $date")
    # Calculate the magnetic field vector using the World Magnetic Model.
    # The result is in the NED frame and has units of nT.
    B_ned = igrf(yeardecimal(date), alt_m, lat_rad, lon_rad, Val(:geodetic))
    B_pp = ned_to_ecef(B_ned, lat_rad, lon_rad, alt_m)
    B_ii = L_PI' * B_pp
    return B_ii
end

"""
    calculate_magnetic_torque(m::AbstractVector, B::AbstractVector)

Calculates the magnetic torque exerted on a magnetic dipole by an external
magnetic field.

**τ** = **m** × **B**

The magnetic dipole moment `m` and the magnetic field `B` vectors **must** be
expressed in the same reference frame. The resulting torque vector `τ` will be
in that same frame.

# Args

- `m`: The magnetic dipole moment vector of the spacecraft [A·m²].
- `B`: The external magnetic field vector [Tesla].

# Returns

- A 3-element `SVector` representing the magnetic torque `[τ_x, τ_y, τ_z]` [N·m].
"""
function calculate_magnetic_torque(m::AbstractVector, B::AbstractVector)
    # Ensure inputs are StaticArrays for performance
    m_svector = SVector{3, Float64}(m)
    B_svector = SVector{3, Float64}(B)

    # Calculate the torque using the cross product
    # τ = m × B
    τ = cross(m_svector, B_svector)

    return τ
end

function eclipse_area_calc(r_sat::SVector{3, Float64}, r_sun::SVector{3, Float64}, rp::Float64)
    """
    Calculate the exposed area of the satellite. Translated from Python to Julia. 
    See equations and diagrams in Basilisk documentation: https://avslab.github.io/basilisk/Documentation/simulation/environment/eclipse/eclipse.html

    Parameters 
    ----------
    r_sat : Vector{Float64}
        Position vector of the satellite relative to the planet.
    r_sun : Vector{Float64}
        Position vector of the Sun relative to the planet.
    A : Float64
        Area of the satellite.
    rp : Float64
        Radius of the planet.


    Returns
    -------
    eclipse_ratio : Float64
        Eclipse ratio of the satellite.
    
    """
    rs = 695000e3 # Radius of the Sun in meters
    if dot(r_sun, r_sat) >= 0 # check the cos of the angle between the satellite and the Sun. If positive (angle less than 90 degrees), the satellite is not in eclipse
        return 1.0 # If the satellite is not in eclipse, return 1.0
    end
    # Eclipse conditions
    f1 = asin((rs + rp) / norm(r_sun)) # Penumbra angle
    f2 = asin((rs - rp) / norm(r_sun)) # Umbra angle
    s0 = -dot(r_sat, r_sun) / (norm(r_sun)) # Plane-axis intersection and planet center distance
    c1 = s0 + rp * sin(f1) # Distance from fundamental plane to cone vertex V1
    c2 = s0 - rp * sin(f2) # Distance from fundamental plane to cone vertex V2
    l1 = c1*tan(f1) # Radius of penumbra cone in fundamental plane
    l2 = c2*tan(f2) # Radius of umbra cone in fundamental plane
    l = √(norm(r_sat)^2 - s0^2) # Distance from fundamental plane to satellite
    
    # Apparent radii of sun, planet, and apparent separation of sun and planet, respectively
    a = asin(rs / norm(r_sun)) # Apparent radius of the Sun
    b = asin(rp / norm(r_sat)) # Apparent radius of the planet
    c = acos(-dot(r_sun, r_sat) / (norm(r_sun) * norm(r_sat))) # Apparent separation of the Sun and planet
    if c < b - a # Total eclipse condition
        return 0.0 # If the satellite is in total eclipse, return 0.0
    elseif c < a - b # Annular eclipse condition
        A_sun = π * a^2 # Apparent area of the Sun
        A_planet = π * b^2 # Apparent area of the planet
        return b^2 / a^2 # Shadow fraction
    elseif c < a + b # Partial eclipse condition
        x = (c^2 + a^2 - b^2) / (2 * c)
        y = √(a^2 - x^2)
        A = a^2 * acos(x / a) + b^2 * acos((c - x) / b) - c * y
        return 1 - A / (π * a^2) # Shadow fraction
    else # No eclipse condition
        return 1.0 # If the satellite is not in eclipse, return 1.0
    end

    # shadow = "none"
    # rs = 6.9634e8 # Radius of the Sun in meters 
    # RP = norm(r_sun) # Distance from Sun to the planet 
    # alpha_umb = asin((rs - rp) / RP)  # Umbra angle
    # alpha_pen = asin((rs + rp) / RP)  # Penumbra angle

    # if dot(r_sun, r_sat) < 0 # if the angle is greater than 90 degrees, satellite is potentially in an eclipse
    #     # Compute the satellite's horizontal and vertical distances
    #     sigma = acos(dot(-r_sat, r_sun) / (norm(r_sat) * norm(r_sun)))
    #     sat_horiz = norm(r_sat) * cos(sigma)
    #     sat_vert = norm(r_sat) * sin(sigma)
    #     # Determine the eclipse conditions
    #     x = rp / sin(alpha_pen)
    #     pen_vert = tan(alpha_pen) * (x + sat_horiz)

    #     if sat_vert <= pen_vert # if true, the satellite is in partial shadow(penumbra)
    #         shadow = "penumbra"
    #         y = rp / sin(alpha_umb)
    #         umb_vert = tan(alpha_umb) * (y - sat_horiz)
    #         if sat_vert <= umb_vert
    #             shadow = "umbra"
    #         end
    #     end
    # end

    # if shadow == "none"
    #     eclipse_ratio = 1
    # elseif shadow == "penumbra"
    #     eclipse_ratio = (1 - (1 - (sat_vert / pen_vert))^2)
    # else 
    #     eclipse_ratio = 0
    # end
    
    # return eclipse_ratio
end

# function srp(p::config.Planet, facet::config.Facet, b::config.Body)
#     """
#     Calculate SRP force on a single facet, i.e., a single face of a rigid body.

#     Parameters
#     ----------
#     p : Planet struct
#         Contains planetary parameters, including equatorial radius.
#     facet : Facet struct
#         Contains information about the facet.
#     b : Body struct
#         Contains physical information about the entire rigid body.
    
#     Returns
#     -------
#     F_srp : SVector{3, Float64}
#         Force on the facet in the inertial frame
#     """

# end

function srp!(model, root_index::Int64, sun_dir_ii::SVector{3, Float64}, body, P_srp::Float64, eclipse_ratio::Float64, orientation::Bool)
    """
    Calculate force on a body due to solar radiation pressure.

    Parameters
    ----------
    pos_ii : SVector{3, Float64}
        Position of the body in the inertial frame (J2000)
    sun_dir_ii : SVector{3, Float64}
        Unit vector in the direction of the Sun expressed in the inertial frame
    body : Body struct
        Struct containing physical information about the body
    r_sun_norm : Float64
        Magnitude of the spacecraft distance to the Sun
    P_srp : Float64
        Magnitude of the solar radiation pressure force at r_sun_norm meters from the Sun
    
    Returns
    -------
    F_srp : SVector{3, Float64}
        Force on the body in the inertial frame
    """
    rot_inertial = config.rotate_to_inertial(model, body, root_index)
    rot_body_to_inertial = rot(model.links[root_index].q)
    @inbounds for facet in body.SRP_facets
        rot_RF = rot_inertial * rot(facet.attitude)' # Rotation matrix from facet frame to inertial frame
        n = normalize(rot_RF * facet.normal_vector) # Normal vector of the facet in the inertial frame
        cos_α_srp = dot(n, sun_dir_ii) / norm(n) / norm(sun_dir_ii)

        if cos_α_srp > 0 && eclipse_ratio != 0.0 # If the facet is illuminated by the Sun
            F_SRP = -P_srp * facet.area * cos_α_srp * ((1 - facet.δ) * sun_dir_ii + 2 * (facet.ρ / 3 + facet.δ * cos_α_srp) * n) * eclipse_ratio
            body.net_force += F_SRP # Rotate F_SRP from body frame to inertial frame

            if orientation
                R_facet = rot_inertial*facet.cp + rot_body_to_inertial'*body.r # Vector from CoM of spacecraft to facet Cp in inertial frame
                # R_facet_body = config.rotate_to_body(body)*facet.cp + body.r
                body.net_torque += rot_body_to_inertial * cross(R_facet, F_SRP) # Calculate body frame net torque
            end
        end
    end
    # CSV.write("facet_forces.csv", df)
    # return F_SRP_tracker
    # println("Total F_SRP: $F_SRP_tracker")
end

