# include("Ref_system_conf.jl")

using LinearAlgebra
using AstroTime
using StaticArrays
using SatelliteToolboxTransformations
using SatelliteToolbox
using .SimulationModel: Model, Planet
# eop_iau2000a = fetch_iers_eop(Val(:IAU2000A))

function r_intor_p!(r_i::SVector{3, Float64}, v_i::SVector{3, Float64}, planet, et)
    # From PCI (planet centered inertial) to PCPF (planet centered/planet fixed)

    # rot_angle = norm(planet.ω) * (t + t_prev)
    # et = utc2et(to_utc(DateTime(date_initial + t0*seconds)))
    # current_time = value(seconds(date_initial + t0*seconds - from_utc(DateTime(2000, 1, 1, 12, 0, 0.0)))) # current time in seconds since J2000
    # et = utc2et(to_utc(current_time)) # convert to Julian Ephemeris Time (ET)
    # primary_body_name = planet.name
    # planet.L_PI .= SMatrix{3, 3, Float64}(pxform("J2000", "IAU_"*uppercase(primary_body_name), et))*planet.J2000_to_pci' # Construct a rotation matrix from J2000 (Planet-fixed frame 0.0 seconds past the J2000 epoch) to planet-fixed frame
    
    # L_pi = [cos(rot_angle) sin(rot_angle) 0; 
    #         -sin(rot_angle) cos(rot_angle) 0; 
    #         0 0 1]

    r_p = SVector{3, Float64}(planet.L_PI * r_i)

    v_p = SVector{3, Float64}(planet.L_PI * (v_i - cross(planet.ω, r_i)))

    return r_p, v_p
end

function r_pintor_i(r_p::SVector{3, Float64}, v_p::SVector{3, Float64}, planet)
    # From PCPF (planet centered/planet fixed) to PCI (planet centered inertial)

    r_i = planet.L_PI' * r_p
    v_i = planet.L_PI' * SVector{3, Float64}(v_p + cross(planet.ω, r_p))

    return r_i, v_i
end

function orbitalelemtorv(oe::SVector{7, Float64}, planet)
    # From orbital element to ECI (Planet Centered Inertial)

    a, e, i, Ω, ω, vi = oe[1], oe[2], oe[3], oe[4], oe[5], oe[6]

    p = a*(1 - e^2)
    h = sqrt(planet.μ * p)

    r_x = (h^2) / planet.μ * (1 / (1 + e * cos(vi))) * SVector{3, Float64}([cos(vi); sin(vi); 0])
    v_x = planet.μ / h * SVector{3, Float64}([-sin(vi); e + cos(vi); 0])
    
    Q = SMatrix{3, 3, Float64}([-sin(Ω)*cos(i)*sin(ω)+cos(Ω)*cos(ω) cos(Ω)*cos(i)*sin(ω)+sin(Ω)*cos(ω) sin(i)*sin(ω); 
         -sin(Ω)*cos(i)*cos(ω)-cos(Ω)*sin(ω) cos(Ω)*cos(i)*cos(ω)-sin(Ω)*sin(ω) sin(i)*cos(ω);
          sin(Ω)*sin(i) -cos(Ω)*sin(i) cos(i)])

    R = Q' * r_x
    V = Q' * v_x

    return collect(R), collect(V)
end

function rvtoorbitalelement(r::SVector, v::SVector, m::Float64, planet::Planet)
    # From ECI (Planet Centered Inertial) to orbital element
    i_x = SVector{3, Int64}([1, 0, 0])
    i_y = SVector{3, Int64}([0, 1, 0])
    i_z = SVector{3, Int64}([0, 0, 1])

    Energy = dot(v,v)/2 - planet.μ/norm(r)
    a = -planet.μ / (2 * Energy)
    h = cross(r, v)
    index = 0

    r_ver = r/norm(r)
    e = cross(v, h)/planet.μ - r_ver

    i = acos(dot(i_z, h)/norm(h))
    e_vers = e / norm(e)
    if i == 0 || i == pi
        if dot(e, i_y) >= 0
            periapsis_longitude = acos(dot(i_x, e_vers))
        elseif dot(e, i_y) < 0
            periapsis_longitude = 2*pi - acos(dot(i_x, e_vers))
        end
        Ω = periapsis_longitude
        ω = 0
    else
        n = cross(i_z, h)/norm(cross(i_z, h))
        if dot(n, i_y) >= 0
            Ω = acos(dot(i_x, n))
        elseif dot(n, i_y) < 0
            Ω = 2*pi - acos(dot(i_x, n))
        end

        if dot(e, i_z) >= 0
            if dot(n, e_vers) > 1 && dot(n, e_vers) < 1 + 1e-4
                ω = acos(1)
            elseif dot(n, e_vers) < -1 && dot(n, e_vers) > -1 - 1e-4
                ω = acos(-1)
            else
                ω = acos(dot(n, e_vers))
            end
        elseif dot(e, i_z) < 0
            if dot(n, e_vers) > 1 && dot(n, e_vers) < 1 + 1e-4
                ω = 2*pi - acos(1)
            elseif dot(n, e_vers) < -1 && dot(n, e_vers) > -1 - 1e-4
                ω = 2*pi - acos(-1)
            else
                ω = 2*pi - acos(dot(n, e_vers))
            end
        end
    end

    if dot(r, v) > 0
        value = dot(e_vers, r_ver)
        if abs(value) >= 1
            value = round(value)
        end
        vi = acos(value)
    elseif dot(r, v) <= 0
        value = dot(e_vers, r_ver)
        if abs(value) >= 1
            value = round(value)
        end
        vi = 2*pi - acos(value)
    end

    e = norm(e)
    if Ω == pi
        Ω = 0
    end

    return SVector{7, Float64}([a, e, i, Ω, ω, vi, m])
end

function rtoalfadeltar(r)
    # From PCI (Planet Centered Inertial) to Geocentric Celestial Reference Frame (GCRF)
    # Conversion between x,y,z and right ascension (RA), declination (dec), distance from the center of the planet (r)
    x = r[1]
    y = r[2]
    z = r[3]
    r = sqrt(x^2 + y^2 + z^2)
    l, m, n = x/r, y/r, z/r

    dec = asin(n)
    if m > 0
        RA = acos(l/cos(dec))
    else
        RA = 2*pi - acos(l/cos(dec))
    end

    return [r, RA, dec]
end

function alfadeltartor(R_RA_DEC)
    # From Geocentric Celestial Reference Frame (GCRF) to PCI (Planet Centered Inertial)
    R = R_RA_DEC[1]
    RA = R_RA_DEC[2]
    DEC = R_RA_DEC[3]

    x = R * cos(DEC) * cos(RA)
    y = R * cos(DEC) * sin(RA)
    z = R * sin(DEC)

    return [x, y, z]
end

function latlongtor(LATLONGH, planet, α_g0, t, t0)
    # From Geodetic to PCI
    ϕ = LATLONGH[1]
    λ = LATLONGH[2]
    h = LATLONGH[3]

    a = planet.Rp_e
    b = planet.Rp_p
    e = sqrt(1 - b^2/a^2)
    α = λ + α_g0 + planet.ω[3]*(t - t0)
    cnst = a / (1 - e^2 * sin(ϕ)^2) + h

    x = cnst * cos(ϕ) * cos(α)
    y = cnst * cos(ϕ) * sin(α)
    z = cnst * sin(ϕ)

    return [x, y, z]
end

function latlongtoOE(LATLONGH, planet, γ, α, v)
    """
    From Geodetic to Orbital Elements
    
    This function converts geodetic coordinates (latitude, longitude, and altitude) to orbital elements.
    It uses the planet's parameters and the velocity vector to compute the orbital elements.
    The function returns the orbital elements as a vector.
    Parameters:
    - LATLONGH: A vector containing latitude, longitude, and altitude.
    - planet: A structure containing the planet's parameters.
    - γ: The flight path angle.
    - α: The azimuth of the velocity vector.
    - v: The velocity magnitude.
    Returns:
    - A vector containing the orbital elements: semi-major axis, eccentricity, inclination, right ascension of the ascending node, argument of periapsis, and true anomaly.
    """
    # Geodetic to Orbital Elements
    ϕ = LATLONGH[1]
    λ = LATLONGH[2]
    h = LATLONGH[3]
    println("ϕ: ", ϕ)
    println("λ: ", λ)
    println("h: ", h)
    println("γ: ", γ)
    println("α: ", α)
    println("v: ", v)
    
    f = (planet.Rp_e - planet.Rp_p) / planet.Rp_e
    # e = 1 - (1 - f)^2 # ellipticity (NOTE =  considered as square)
    e = sqrt(1 - (planet.Rp_p/planet.Rp_e)^2) # ellipticity (NOTE =  considered as square)
    N = planet.Rp_e / sqrt(1 - e^2 * sin(ϕ)^2)
    x = (N + h) * cos(ϕ) * cos(λ)
    y = (N + h) * cos(ϕ) * sin(λ)
    z = ((1 - e^2) * N + h) * sin(ϕ)
    r = sqrt(x^2 + y^2 + z^2)
    println("x: ", x)
    println("y: ", y)
    println("z: ", z)
    r_i, _ = r_pintor_i(SVector{3, Float64}([x, y, z]), SVector{3, Float64}([0, 0, 0]), planet)
    # println("r_i: ", r_i)
    # Define local unit vectors in ECEF:
    Z_ecef = SVector{3,Float64}([cos(ϕ)*cos(λ), cos(ϕ)*sin(λ), sin(ϕ)])
    E_ecef = SVector{3,Float64}([-sin(λ), cos(λ), 0])
    N_ecef = -SVector{3,Float64}([-sin(ϕ)*cos(λ), -sin(ϕ)*sin(λ), cos(ϕ)])  # since S = [-sinϕ*cosλ, -sinϕ*sinλ, cosϕ]

    # Alternatively, you can define N_ecef directly as:
    N_ecef = SVector{3,Float64}([-cos(λ)*sin(ϕ), -sin(λ)*sin(ϕ), cos(ϕ)])

    # Decompose the local velocity (assumed provided as v magnitude, γ, and α):
    v_N = v * cos(γ) * cos(α)  # North component
    v_E = v * cos(γ) * sin(α)  # East component
    v_Z = v * sin(γ)           # Up/Zenith component

    # Form the velocity in ECEF:
    v_ecef = v_N * N_ecef + v_E * E_ecef + v_Z * Z_ecef

    # Convert to ECI using your transformation matrix:
    v_eci = planet.L_PI' * v_ecef
    # v_vec = SVector{3, Float64}([v*cos(γ)*cos(α), v*cos(γ)*sin(α), v*sin(γ)])

    OE = rvtoorbitalelement(r_i, v_eci, 0, planet)[1:6]
    return OE
end

function rtolatlong(r_p::SVector{3, Float64}, planet, spherical_harmonic_topography::Bool=false)
    # From PCPF to LLA through Bowring's method https://www.mathworks.com/help/aeroblks/ecefpositiontolla.html;jsessionid=2ae36964c7d5f2115d2c21286db0?nocookie=true
    x_p = r_p[1]
    y_p = r_p[2]
    z_p = r_p[3]

    f = (planet.Rp_e - planet.Rp_p) / planet.Rp_e
    e = 1 - (1 - f)^2 # ellipticity (NOTE =  considered as square)
    r = sqrt(x_p^2 + y_p^2)

    # Calculate initial guesses for reduced latitude (latr) and planet-detic latitude (latd)
    latr = atan(z_p / ((1-f)*r)) # reduced latitude
    latd = atan((z_p + (e*(1-f)*planet.Rp_e*sin(latr)^3)/(1 - e)) / (r - e*planet.Rp_e*cos(latr)^3))

    # Recalculate reduced latitude based on planet-detic latitude
    latr2 = atan((1 - f)*sin(latd) / cos(latd))
    diff = latr - latr2

    # Iterate until reduced latitude converges
    while diff > 1e-10
        latr = latr2
        latd = atan((z_p + (e*(1-f)*planet.Rp_e*sin(latr)^3) / (1 - e)) / (r - e*planet.Rp_e*cos(latr)^3))
        latr2 = atan((1 - f)*sin(latd) / cos(latd))
        diff = latr - latr2
    end
    lat = latd

    #Calculate longitude
    lon = atan(y_p, x_p)
    # Calculate Altitude
    if !spherical_harmonic_topography
        N = planet.Rp_e / sqrt(1 - e*sin(lat)^2)
        alt = r*cos(lat) + (z_p + e*N*sin(lat)^2)*sin(lat) - N
    else
        alt = norm(r_p) - planet.topography_function(args, 
                                                    planet.Clm_topo, 
                                                    planet.Slm_topo, 
                                                    lat, 
                                                    lon,
                                                    planet.A_topo)
    end
    
    return SVector{3, Float64}([alt, lat, lon])
end
function rtolatlongrad(r_p, planet)
    # Same as previous function, but returns radius instead of altitude and planetocentric latitude and longitude
    # From PCPF to LLA through Bowring's method https://www.mathworks.com/help/aeroblks/ecefpositiontolla.html;jsessionid=2ae36964c7d5f2115d2c21286db0?nocookie=true
    x_p = r_p[1]
    y_p = r_p[2]
    z_p = r_p[3]

    r = sqrt(x_p^2 + y_p^2)
    lat = asin(z_p / norm(r_p))

    #Calculate longitude
    lon = atan(y_p, x_p)
    
    r = norm(r_p)
    
    return SVector{3, Float64}([r, lat, lon])
end

function latlongtoNED(H_LAN_LON)
    lon = H_LAN_LON[3]
    lat = H_LAN_LON[2]

    # Compute first in xyz coordinates(z: north pole, x - z plane: contains r, y: completes right - handed set)
    uDxyz = SVector{3, Float64}([-cos(lat), 0.0, -sin(lat)])
    uNxyz = SVector{3, Float64}([-sin(lat), 0.0, cos(lat)])
    uExyz = SVector{3, Float64}([0.0, 1.0, 0.0])

    # Rotate by longitude to change to PCPF frame
    L3 = SMatrix{3, 3, Float64}([cos(lon) -sin(lon) 0.0;
          sin(lon) cos(lon) 0.0;
          0.0 0.0 1.0])

    uN = L3 * uNxyz
    uE = L3 * uExyz
    uD = L3 * uDxyz

    return SVector{3, SVector{3, Float64}}([uD, uN, uE])
end

# function inertialToLVLHFromOE(a::Float64, e::Float64, i::Float64, ω::Float64, Ω::Float64, ν::Float64)::SVector{4, Float64}
#     """
#     Determine the quaternion from the inertial frame to the LVLH frame based on the orbital elements

#     Args:
#         a::Float64 : Semimajor axis
#         e::Float64 : Eccentricity
#         i::Float64 : Inclination, rad
#         ω::Float64 : Argument of periapsis, rad
#         Ω::Float64 : RAAN, rad 
#         ν::Float64 : True anomaly, rad

#     Returns:
#         q::SVector{4, Float64} : Quaternion expressing rotation from inertial to body frame
    # """


"""
    orbital_elements_to_lvlh_quaternion(raan, inclination, arg_of_perigee, true_anomaly)

Calculates the quaternion representing the rotation from the inertial frame to the 
Local-Vertical/Local-Horizontal (LVLH) frame.

The LVLH frame is defined as:
- Z-axis: Points from the satellite to the Earth's center (nadir).
- Y-axis: Opposite to the orbital angular momentum vector.
- X-axis: Completes the right-handed system (generally in the velocity direction).

The resulting quaternion `q` transforms a vector `v_inertial` to `v_lvlh`.

# Arguments
- `raan::Float64`: Right Ascension of the Ascending Node (Ω) in radians.
- `inclination::Float64`: Inclination (i) in radians.
- `arg_of_perigee::Float64`: Argument of Perigee (ω) in radians.
- `true_anomaly::Float64`: True Anomaly (ν) in radians.

# Returns
- `Vector{Float64}`: A 4-element quaternion `[qx, qy, qz, qw]`.
"""
function orbital_elements_to_lvlh_quaternion(
    raan::Float64,
    inclination::Float64,
    arg_of_perigee::Float64,
    true_anomaly::Float64
)::SVector{4, Float64}

    # 1. Calculate Argument of Latitude (u)
    u = arg_of_perigee + true_anomaly

    # Pre-compute sines and cosines for clarity and efficiency
    su = sin(u)
    cu = cos(u)
    si = sin(inclination)
    ci = cos(inclination)
    sO = sin(raan)
    cO = cos(raan)

    # 2. Find LVLH axes in the inertial frame
    
    # Position unit vector (r_hat) in inertial frame
    r_hat_inertial = SVector{3, Float64}([
        cO * cu - sO * su * ci;
        sO * cu + cO * su * ci;
        su * si
    ])

    # Angular momentum unit vector (h_hat) in inertial frame
    h_hat_inertial = SVector{3, Float64}([
        sO * si;
        -cO * si;
        ci
    ])

    # Define LVLH axes based on r_hat and h_hat
    z_lvlh = -r_hat_inertial
    y_lvlh = -h_hat_inertial
    x_lvlh = cross(y_lvlh, z_lvlh)
    
    # Normalize to ensure perfect orthogonality due to potential floating point errors
    x_lvlh = normalize(x_lvlh)

    # 3. Construct the Direction Cosine Matrix (DCM) from Inertial to LVLH
    # The rows of the DCM are the basis vectors of the new frame (LVLH)
    # expressed in the old frame's (Inertial) coordinates.
    C = SMatrix{3, 3, Float64}([
        x_lvlh';
        y_lvlh';
        z_lvlh'
    ])

    # 4. Convert DCM to a quaternion (numerically stable method)
    trace = tr(C) # Trace of the matrix
    
    if trace > 0
        S = sqrt(trace + 1.0) * 2
        qw = 0.25 * S
        qx = (C[3, 2] - C[2, 3]) / S
        qy = (C[1, 3] - C[3, 1]) / S
        qz = (C[2, 1] - C[1, 2]) / S
    elseif (C[1, 1] > C[2, 2]) && (C[1, 1] > C[3, 3])
        S = sqrt(1.0 + C[1, 1] - C[2, 2] - C[3, 3]) * 2
        qw = (C[3, 2] - C[2, 3]) / S
        qx = 0.25 * S
        qy = (C[1, 2] + C[2, 1]) / S
        qz = (C[1, 3] + C[3, 1]) / S
    elseif C[2, 2] > C[3, 3]
        S = sqrt(1.0 + C[2, 2] - C[1, 1] - C[3, 3]) * 2
        qw = (C[1, 3] - C[3, 1]) / S
        qx = (C[1, 2] + C[2, 1]) / S
        qy = 0.25 * S
        qz = (C[2, 3] + C[3, 2]) / S
    else
        S = sqrt(1.0 + C[3, 3] - C[1, 1] - C[2, 2]) * 2
        qw = (C[2, 1] - C[1, 2]) / S
        qx = (C[1, 3] + C[3, 1]) / S
        qy = (C[2, 3] + C[3, 2]) / S
        qz = 0.25 * S
    end

    quaternion = SVector{4, Float64}([qx, qy, qz, qw])
    return normalize(quaternion) # Final normalization for safety
end

"""
    rotate_vector_by_quaternion(v, q)

Rotates a 3D vector `v` by a quaternion `q`.
Quaternion format is [qx, qy, qz, qw].
"""
function rotate_vector_by_quaternion(v::Vector{Float64}, q::Vector{Float64})
    q_vec = q[1:3]
    q_scalar = q[4]
    
    # Using the formula: v' = v + 2 * q_vec x (q_vec x v + q_scalar * v)
    t = 2 * cross(q_vec, v)
    v_rotated = v + q_scalar * t + cross(q_vec, t)
    
    return v_rotated
end

"""
    ned_to_ecef(v_ned::AbstractVector, date::DateTime, lat::Number, lon::Number)

Converts a vector `v_ned` from the local North-East-Down (NED) frame to the
Earth-Centered, Earth-Fixed (ECEF) frame.

# Args

- `v_ned`: A 3-element vector in the NED frame `[North, East, Down]`.
- `date`: The `DateTime` at which the conversion is to be performed. This is
          crucial for determining the Earth's orientation.
- `lat`: The geodetic latitude of the observer [radians].
- `lon`: The longitude of the observer [radians].

# Returns

- A 3-element `SVector` representing the vector in the GCRF frame.
"""
# function ned_to_ecef(v_ned::AbstractVector, date::DateTime, lat::Float64, lon::Float64, alt_m::Float64)
#     # Ensure the input vector is a 3-element SVector for performance.
#     v_ned_svector = SVector{3, Float64}(v_ned)

#     # == Step 1: Get the rotation from NED to ECEF (Earth-Centered, Earth-Fixed)
#     # This rotation depends only on the observer's position (latitude, longitude).
#     # It transforms the local frame to the global Earth-fixed frame.
#     R_NED_to_ECEF = ned_to_ecef(v_ned_svector, lat, lon, alt_m)

#     # == Step 2: Get the rotation from ECEF to GCRF (Inertial)
#     # This rotation depends on the time and accounts for Earth's rotation.
#     # SatelliteToolbox.jl handles the complex calculations of Earth's orientation
#     # (precession, nutation, polar motion) automatically.
#     R_ECEF_to_GCRF = r_ecef_to_eci(DCM, ITRF(), GCRF(), date_to_jd(date), eop_iau2000a)

#     # == Step 3: Combine the rotations and apply to the vector
#     # The order is crucial: first apply the NED->ECEF rotation, then ECEF->GCRF.
#     # v_gcrf = R_ECEF_to_GCRF * R_NED_to_ECEF * v_ned
#     v_gcrf = R_ECEF_to_GCRF * R_NED_to_ECEF

#     return v_gcrf
# end