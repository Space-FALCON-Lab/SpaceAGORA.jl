module Geodesy

using StaticArrays

# Shared scalar Bowring expressions. Callers retain their input construction,
# flattening precision, prime-vertical radius, and altitude/topography policy.
@inline function _bowring_latitude(z, p_xy, e2, ep2, planet)
    θ = atan(z * planet.Rp_e, p_xy * planet.Rp_p)
    return atan(z + ep2 * planet.Rp_p * sin(θ)^3, p_xy - e2 * planet.Rp_e * cos(θ)^3)
end

@inline function _bowring_altitude(p_xy, z, e2, N, lat)
    return p_xy * cos(lat) + (z + e2 * N * sin(lat)) * sin(lat) - N
end

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
    lat = _bowring_latitude(z, p_xy, e2, ep2, planet)
    N = planet.Rp_e / sqrt(1.0 - e2 * sin(lat)^2)
    return _bowring_altitude(p_xy, z, e2, N, lat)
end

# The planet-fixed position (m) at geodetic latitude `lat_rad`, longitude
# `lon_rad` and altitude `alt_m` above the planet's ellipsoid (`Rp_e`, `Rp_p`):
# the inverse of `rtolatlong`'s conversion, with its flattening.
@inline function geodetic_to_pp(lat_rad::Float64, lon_rad::Float64, alt_m::Float64, planet)::SVector{3, Float64}
    f = (planet.Rp_e - planet.Rp_p) / planet.Rp_e
    e2 = 1.0 - (1.0 - f)^2
    N = planet.Rp_e / sqrt(1.0 - e2 * sin(lat_rad)^2)
    return SVector{3, Float64}((N + alt_m) * cos(lat_rad) * cos(lon_rad), (N + alt_m) * cos(lat_rad) * sin(lon_rad),
        (N * (1.0 - e2) + alt_m) * sin(lat_rad))
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
