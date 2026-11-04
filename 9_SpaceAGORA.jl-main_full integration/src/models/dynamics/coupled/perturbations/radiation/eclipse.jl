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
    @inline _clamp_unit(x::Float64) = clamp(x, -1.0, 1.0)

    r_sun_norm = norm(r_sun)
    r_sat_norm = norm(r_sat)
    if !isfinite(r_sun_norm) || !isfinite(r_sat_norm) || r_sun_norm <= 0.0 || r_sat_norm <= 0.0
        return 1.0
    end

    # Precompute dot product and reciprocals used multiple times below.
    dot_rs_rsat = dot(r_sun, r_sat)
    if dot_rs_rsat >= 0 # check the cos of the angle between the satellite and the Sun. If positive (angle less than 90 degrees), the satellite is not in eclipse
        return 1.0 # If the satellite is not in eclipse, return 1.0
    end
    inv_r_sun_norm = 1.0 / r_sun_norm
    inv_r_sat_norm = 1.0 / r_sat_norm
    # Eclipse conditions
    f1 = asin(_clamp_unit((rs + rp) * inv_r_sun_norm)) # Penumbra angle
    f2 = asin(_clamp_unit((rs - rp) * inv_r_sun_norm)) # Umbra angle
    s0 = -dot_rs_rsat * inv_r_sun_norm # Plane-axis intersection and planet center distance
    c1 = s0 + rp * sin(f1) # Distance from fundamental plane to cone vertex V1
    c2 = s0 - rp * sin(f2) # Distance from fundamental plane to cone vertex V2
    l1 = c1*tan(f1) # Radius of penumbra cone in fundamental plane
    l2 = c2*tan(f2) # Radius of umbra cone in fundamental plane
    l = √(max(r_sat_norm^2 - s0^2, 0.0)) # Distance from fundamental plane to satellite

    # Apparent radii of sun, planet, and apparent separation of sun and planet, respectively
    a = asin(_clamp_unit(rs * inv_r_sun_norm)) # Apparent radius of the Sun
    b = asin(_clamp_unit(rp * inv_r_sat_norm)) # Apparent radius of the planet
    c = acos(_clamp_unit(-dot_rs_rsat * inv_r_sun_norm * inv_r_sat_norm)) # Apparent separation of the Sun and planet
    if c < b - a # Total eclipse condition
        return 0.0 # If the satellite is in total eclipse, return 0.0
    elseif c < a - b # Annular eclipse condition
        return 1.0 - b^2 / a^2 # Exposed sun fraction
    elseif c < a + b # Partial eclipse condition
        x = (c^2 + a^2 - b^2) / (2 * c)
        y = √(max(a^2 - x^2, 0.0))
        A = a^2 * acos(_clamp_unit(x / a)) + b^2 * acos(_clamp_unit((c - x) / b)) - c * y
        return 1 - A / (π * a^2) # Shadow fraction
    else # No eclipse condition
        return 1.0 # If the satellite is not in eclipse, return 1.0
    end
end
