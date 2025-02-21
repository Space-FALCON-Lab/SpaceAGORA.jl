using SPICE

import .config

function gravity_n_bodies(et, pos_ii, p, n_body, el_time)

    primary_body_name = p.name
    n_body_name = n_body.name

    if cmp(lowercase(primary_body_name), "mars") == 0 || cmp(lowercase(primary_body_name), "jupiter") == 0 || cmp(lowercase(primary_body_name), "saturn") == 0 || cmp(lowercase(primary_body_name), "uranus") == 0 || cmp(lowercase(primary_body_name), "neptune") == 0
        primary_body_name *= "_barycenter"
    end

    if cmp(lowercase(n_body_name), "mars") == 0 || cmp(lowercase(n_body_name), "jupiter") == 0 || cmp(lowercase(n_body_name), "saturn") == 0 || cmp(lowercase(n_body_name), "uranus") == 0 || cmp(lowercase(n_body_name), "neptune") == 0
        n_body_name *= "_barycenter"
    end

    pos_primary_k = spkpos(n_body_name, et, "ECLIPJ2000", "none", primary_body_name)[1] * 1e3 # pxform("ECLIPJ2000", "IAU_"*uppercase(p.name), -3.265e7)*
    # println("pos_primary_k: ", pos_primary_k)
    pos_spacecraft_k = pos_primary_k - pos_ii
    pos_spacecraft_k_mag = norm(pos_spacecraft_k)

    g = n_body.μ * ((pos_spacecraft_k / pos_spacecraft_k_mag^3) - (pos_primary_k / norm(pos_primary_k)^3))

    return g
end

function eclipse_area_calc(r_sat::Vector{Float64}, r_sun::Vector{Float64}, A::Float64, rp::Float64)
    """
    Calculate the exposed area of the satellite. Translated from Python to Julia.

    Parameters 
    ----------
    r_sat : Vector{Float64}
        Position vector of the satellite relative to the planet.
    r_sun : Vectro{Float64}
        Position vector of the Sun relative to the planet.
    A : Float64
        Area of the satellite.
    rp : Float64
        Radius of the planet.


    Returns
    -------
    A_exp : Float64
        Exposed area of the satellite.
    
    """
    shadow = "none"
    rs = 6.9634e8 # Radius of the Sun in meters 
    RP = norm(r_sun) # Distance from Sun to the planet 
    alpha_umb = asin((rs - rp) / RP)  # Umbra angle
    alpha_pen = asin((rs + rp) / RP)  # Penumbra angle

    if dot(r_sun, r_sat) < 0 # if the angle is greater than 90 degrees, satellite is potentially in an eclipse
        # Compute the satellite's horizontal and vertical distances
        sigma = acos(dot(-r_sat, r_sun) / (norm(r_sat) * norm(r_sun)))
        sat_horiz = norm(r_sat) * cos(sigma)
        sat_vert = norm(r_sat) * sin(sigma)
        # Determine the eclipse conditions
        x = rp / sin(alpha_pen)
        pen_vert = tan(alpha_pen) * (x + sat_horiz)

        if sat_vert <= pen_vert # if true, the satellite is in partial shadow(penumbra)
            shadow = "penumbra"
            y = rp / sin(alpha_umb)
            umb_vert = tan(alpha_umb) * (y - sat_horiz)
            if sat_vert <= umb_vert
                shadow = "umbra"
            end
        end
    end

    if shadow == "none"
        A_exp = A
    elseif shadow == "penumbra"
        A_exp = A * (1 - (1 - (sat_vert / pen_vert))^2)
    else 
        A_exp = 0
    end
    
    return A_exp
end

function srp(p, p_srp_unscaled::Float64, cR::Float64, A_sat::Float64, m::Float64, r_sat::Vector{Float64}, time_et::Float64)
    """
    Calculate acceleration due to solar radiation pressure. 

    Parameters 
    ----------
    p : Struct
        Contains planetary parameters, including equatorial radius.
    p_srp_unscaled : Float64
        Solar radiation pressure at 1 AU.
    cR : Float64 
        Coefficient of reflectivity.
    A_sat : Float64
        Area of the satellite
    m : Float64
        Mass of the satellite
    r_sat : Vector{Float64}
        Position vector of the satellite relative to the planet.
    time_et : String
        Ephemeris time of the simulation, used for SPICE calculations.
    
    Returns
    -------
    srp_accel : Vector{Float64}
        Acceleration due to solar radiation pressure.
    """
    rp = p.Rp_e # Equatorial radius of the planet
    r_sun, ltime = spkpos("SUN", time_et, "ECLIPJ2000", "NONE", uppercase(p.name))
    r_sun = r_sun .* 1e3 # Convert from km to m

    A_exp = eclipse_area_calc(r_sat, r_sun, A_sat, rp)
    p_srp = p_srp_unscaled * norm(r_sun - r_sat) / (1.496e11) # Scale SRP for distance from Sun
    r_sat_to_sun = r_sun - r_sat
    r_sat_to_sun_mag = norm(r_sat_to_sun)

    srp_accel = -p_srp * cR * A_exp / m * (r_sat_to_sun / r_sat_to_sun_mag)
    return srp_accel
end