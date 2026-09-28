module Geodesy

using StaticArrays

# Shared inversion for altitude-based initial conditions and apsis targeting.
# Directions are unit vectors in planet-fixed coordinates; callers own validation.
@inline function geodetic_altitude(radius::Float64, u_pp::SVector{3, Float64}, planet)::Float64
    x = radius * u_pp[1]
    y = radius * u_pp[2]
    z = radius * u_pp[3]

    f = (planet.Rp_e - planet.Rp_p) / planet.Rp_e
    e2 = 1.0 - (1.0 - f)^2
    ep2 = e2 / (1.0 - e2)
    p_xy = sqrt(x^2 + y^2)
    θ = atan(z * planet.Rp_e, p_xy * planet.Rp_p)
    lat = atan(z + ep2 * planet.Rp_p * sin(θ)^3, p_xy - e2 * planet.Rp_e * cos(θ)^3)
    N = planet.Rp_e / sqrt(1.0 - e2 * sin(lat)^2)
    return p_xy * cos(lat) + (z + e2 * N * sin(lat)) * sin(lat) - N
end

@inline function ellipsoid_surface_radius(u_pp::SVector{3, Float64}, planet)::Float64
    return inv(sqrt((u_pp[1]^2 + u_pp[2]^2) / planet.Rp_e^2 + u_pp[3]^2 / planet.Rp_p^2))
end

function radius_for_geodetic_altitude(target_altitude_m::Float64, u_pp::SVector{3, Float64}, planet)::Float64
    target_altitude_m >= 0.0 || return NaN
    lo = ellipsoid_surface_radius(u_pp, planet)
    hi = lo + target_altitude_m + abs(planet.Rp_e - planet.Rp_p) + 1.0
    while geodetic_altitude(hi, u_pp, planet) < target_altitude_m
        hi += max(target_altitude_m, abs(planet.Rp_e - planet.Rp_p), 1.0)
    end
    for _ in 1:80
        mid = 0.5 * (lo + hi)
        if geodetic_altitude(mid, u_pp, planet) < target_altitude_m
            lo = mid
        else
            hi = mid
        end
    end
    return 0.5 * (lo + hi)
end

end # module Geodesy
